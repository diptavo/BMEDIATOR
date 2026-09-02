#include "../plink_ld.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using bmediator::PlinkData;

struct ProteinInfo {
    std::string protein;
    std::string gene;
    std::string sumstat_file;
    int chr = 0;
    int tss = 0;
};

struct CandidateHit {
    std::string rsid;
    std::string a1;
    std::string a2;
    double freq = 0.0;
    double beta = 0.0;
    double se = 0.0;
    double p = 1.0;
    double n = 0.0;
    int chr = 0;
    int bp = 0;
    int bim_idx = -1;
};

struct ComboConfig {
    double p_thresh = 1.0;
    double r2_thresh = 0.0;
    std::string tag;
};

struct ComboWriter {
    ComboConfig cfg;
    std::ofstream pqtl_out;
    std::ofstream info_out;
    std::ofstream summary_out;
    int proteins_seen = 0;
    int proteins_with_hits = 0;
    long long candidates_total = 0;
    long long kept_total = 0;
};

bool parse_int(const std::string& s, int& out) {
    try {
        out = std::stoi(s);
        return true;
    } catch (...) {
        return false;
    }
}

bool parse_double(const std::string& s, double& out) {
    try {
        out = std::stod(s);
        return true;
    } catch (...) {
        return false;
    }
}

std::string trim(const std::string& s) {
    const size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    const size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, '\t')) out.push_back(item);
    return out;
}

std::string uppercase(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

bool derive_other_allele(const std::string& ref, const std::string& alt,
                         const std::string& a1, std::string& a2) {
    if (a1 == alt) {
        a2 = ref;
        return true;
    }
    if (a1 == ref) {
        a2 = alt;
        return true;
    }
    return false;
}

std::vector<ProteinInfo> read_resource(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("cannot open resource file: " + path);

    std::string line;
    if (!std::getline(fin, line)) throw std::runtime_error("empty resource file: " + path);

    auto header = split_tab(line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;

    const std::vector<std::string> required = {"protein", "gene", "sumstat_file", "chr", "tss"};
    for (const auto& key : required) {
        if (!col.count(key)) throw std::runtime_error("missing resource column: " + key);
    }

    std::vector<ProteinInfo> out;
    while (std::getline(fin, line)) {
        if (trim(line).empty()) continue;
        auto fields = split_tab(line);
        if (static_cast<int>(fields.size()) <= col["tss"]) continue;

        ProteinInfo p;
        p.protein = fields[col["protein"]];
        p.gene = fields[col["gene"]];
        p.sumstat_file = fields[col["sumstat_file"]];
        if (!parse_int(fields[col["chr"]], p.chr)) continue;
        if (!parse_int(fields[col["tss"]], p.tss)) continue;
        if (p.protein.empty() || p.sumstat_file.empty()) continue;
        out.push_back(p);
    }
    return out;
}

std::vector<CandidateHit> read_candidates_union(const ProteinInfo& protein,
                                                const PlinkData& plink,
                                                double max_p_thresh,
                                                int flank_kb) {
    std::ifstream fin(protein.sumstat_file);
    if (!fin.is_open()) throw std::runtime_error("cannot open ARIC sumstat: " + protein.sumstat_file);

    std::string line;
    if (!std::getline(fin, line)) return {};

    auto header = split_tab(line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;

    const std::vector<std::string> required = {
        "#CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "P", "ERRCODE"
    };
    for (const auto& key : required) {
        if (!col.count(key)) throw std::runtime_error("missing ARIC column " + key + " in " + protein.sumstat_file);
    }

    const int flank_bp = flank_kb * 1000;
    std::vector<CandidateHit> out;
    while (std::getline(fin, line)) {
        if (trim(line).empty()) continue;
        auto fields = split_tab(line);
        if (static_cast<int>(fields.size()) <= col["ERRCODE"]) continue;
        if (fields[col["TEST"]] != "ADD") continue;
        if (fields[col["ERRCODE"]] != ".") continue;

        CandidateHit hit;
        if (!parse_int(fields[col["#CHROM"]], hit.chr)) continue;
        if (hit.chr != protein.chr) continue;
        if (!parse_int(fields[col["POS"]], hit.bp)) continue;
        if (hit.bp < protein.tss - flank_bp || hit.bp > protein.tss + flank_bp) continue;
        if (!parse_double(fields[col["P"]], hit.p)) continue;
        if (hit.p >= max_p_thresh) continue;
        if (!parse_double(fields[col["BETA"]], hit.beta)) continue;
        if (!parse_double(fields[col["SE"]], hit.se)) continue;
        if (!parse_double(fields[col["A1_FREQ"]], hit.freq)) continue;
        if (!parse_double(fields[col["OBS_CT"]], hit.n)) continue;

        hit.rsid = fields[col["ID"]];
        if (hit.rsid.empty() || hit.rsid == ".") continue;
        hit.a1 = uppercase(fields[col["A1"]]);
        if (!derive_other_allele(uppercase(fields[col["REF"]]), uppercase(fields[col["ALT"]]), hit.a1, hit.a2)) {
            continue;
        }

        auto it = plink.rsid_to_idx.find(hit.rsid);
        if (it == plink.rsid_to_idx.end()) continue;
        hit.bim_idx = it->second;
        out.push_back(hit);
    }

    return out;
}

std::vector<CandidateHit> filter_candidates(const std::vector<CandidateHit>& hits, double p_thresh) {
    std::vector<CandidateHit> out;
    out.reserve(hits.size());
    for (const auto& hit : hits) {
        if (hit.p < p_thresh) out.push_back(hit);
    }
    return out;
}

std::vector<CandidateHit> ld_clump_hits(const PlinkData& plink,
                                        const std::vector<CandidateHit>& hits,
                                        double r2_thresh,
                                        int window_kb) {
    if (hits.empty()) return {};

    std::vector<int> order(hits.size());
    for (int i = 0; i < static_cast<int>(hits.size()); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) { return hits[a].p < hits[b].p; });

    auto geno = plink.extract_genotypes_double([&]() {
        std::vector<int> indices;
        indices.reserve(hits.size());
        for (const auto& hit : hits) indices.push_back(hit.bim_idx);
        return indices;
    }());

    std::vector<double> means(hits.size(), 0.0), sds(hits.size(), 0.0);
    for (int i = 0; i < static_cast<int>(hits.size()); ++i) {
        double s = 0.0, s2 = 0.0;
        for (double g : geno[i]) {
            s += g;
            s2 += g * g;
        }
        means[i] = s / geno[i].size();
        const double var = s2 / geno[i].size() - means[i] * means[i];
        sds[i] = (var > 1e-10) ? std::sqrt(var) : 0.0;
    }

    const int window_bp = window_kb * 1000;
    std::vector<bool> removed(hits.size(), false);
    std::vector<CandidateHit> kept;
    kept.reserve(hits.size());

    for (int oi = 0; oi < static_cast<int>(order.size()); ++oi) {
        const int idx = order[oi];
        if (removed[idx]) continue;
        kept.push_back(hits[idx]);
        for (int oj = oi + 1; oj < static_cast<int>(order.size()); ++oj) {
            const int jdx = order[oj];
            if (removed[jdx]) continue;
            if (hits[idx].chr != hits[jdx].chr) continue;
            if (std::abs(hits[idx].bp - hits[jdx].bp) > window_bp) continue;
            if (sds[idx] < 1e-10 || sds[jdx] < 1e-10) continue;
            double cov = 0.0;
            for (int s = 0; s < plink.n_samples; ++s) {
                cov += (geno[idx][s] - means[idx]) * (geno[jdx][s] - means[jdx]);
            }
            cov /= plink.n_samples;
            const double r = cov / (sds[idx] * sds[jdx]);
            if (r * r > r2_thresh) removed[jdx] = true;
        }
    }

    return kept;
}

std::vector<ComboConfig> combos() {
    return {
        {5e-6, 0.01, "aric_seqid_cis_tss_pm6_r2_0.01"},
        {5e-6, 0.1, "aric_seqid_cis_tss_pm6_r2_0.1"},
        {5e-8, 0.01, "aric_seqid_cis_tss_pm8_r2_0.01"},
        {5e-8, 0.1, "aric_seqid_cis_tss_pm8_r2_0.1"},
    };
}

void open_writer(ComboWriter& writer, const std::string& out_dir, const std::string& chunk_label) {
    const std::string base = out_dir + "/" + writer.cfg.tag + "." + chunk_label;
    writer.pqtl_out.open(base + ".tsv");
    writer.info_out.open(base + ".protein_info.tsv");
    writer.summary_out.open(base + ".build_summary.tsv");
    if (!writer.pqtl_out.is_open() || !writer.info_out.is_open() || !writer.summary_out.is_open()) {
        throw std::runtime_error("failed to open output files for " + base);
    }
    writer.pqtl_out << "PROTEIN\tSNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tN\tCHR\tBP\n";
    writer.info_out << "PROTEIN\tGENE\tCHR\tSTART\tEND\n";
    writer.summary_out << "PROTEIN\tGENE\tSTATUS\tN_CANDIDATES\tN_KEPT\tSOURCE\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: build_aric_pqtl_instruments_multi <resource_subset.tsv> <plink_prefix> <out_dir> <chunk_label>\n";
        return 1;
    }

    const std::string resource_path = argv[1];
    const std::string plink_prefix = argv[2];
    const std::string out_dir = argv[3];
    const std::string chunk_label = argv[4];

    auto proteins = read_resource(resource_path);
    PlinkData plink;
    if (!plink.load(plink_prefix)) {
        std::cerr << "failed to load PLINK reference: " << plink_prefix << "\n";
        return 1;
    }

    std::vector<ComboWriter> writers;
    for (const auto& cfg : combos()) {
        ComboWriter w;
        w.cfg = cfg;
        open_writer(w, out_dir, chunk_label);
        writers.push_back(std::move(w));
    }

    for (int i = 0; i < static_cast<int>(proteins.size()); ++i) {
        const auto& protein = proteins[i];
        std::vector<CandidateHit> union_hits;
        try {
            union_hits = read_candidates_union(protein, plink, 5e-6, 500);
        } catch (const std::exception& e) {
            for (auto& writer : writers) {
                ++writer.proteins_seen;
                writer.summary_out << protein.protein << '\t' << protein.gene << "\tread_error\t0\t0\t"
                                   << protein.sumstat_file << '\n';
            }
            std::cerr << chunk_label << ": " << e.what() << "\n";
            continue;
        }

        for (auto& writer : writers) {
            ++writer.proteins_seen;
            auto filtered = filter_candidates(union_hits, writer.cfg.p_thresh);
            writer.candidates_total += static_cast<long long>(filtered.size());
            auto kept = ld_clump_hits(plink, filtered, writer.cfg.r2_thresh, 10000);
            if (kept.empty()) {
                writer.summary_out << protein.protein << '\t' << protein.gene << "\tno_hits\t"
                                   << filtered.size() << "\t0\t" << protein.sumstat_file << '\n';
                continue;
            }
            ++writer.proteins_with_hits;
            writer.kept_total += static_cast<long long>(kept.size());
            writer.info_out << protein.protein << '\t' << protein.gene << '\t' << protein.chr
                            << '\t' << protein.tss << '\t' << protein.tss << '\n';
            writer.summary_out << protein.protein << '\t' << protein.gene << "\tok\t"
                               << filtered.size() << '\t' << kept.size() << '\t'
                               << protein.sumstat_file << '\n';
            for (const auto& hit : kept) {
                writer.pqtl_out << protein.protein << '\t' << hit.rsid << '\t' << hit.a1 << '\t'
                                << hit.a2 << '\t' << std::setprecision(10) << hit.freq << '\t'
                                << hit.beta << '\t' << hit.se << '\t' << hit.p << '\t'
                                << hit.n << '\t' << hit.chr << '\t' << hit.bp << '\n';
            }
        }

        if ((i + 1) % 100 == 0) {
            std::cerr << chunk_label << ": processed " << (i + 1) << " proteins from subset\n";
        }
    }

    for (const auto& writer : writers) {
        std::cerr << chunk_label << '\t' << writer.cfg.tag
                  << "\tproteins_scanned=" << writer.proteins_seen
                  << "\tproteins_with_hits=" << writer.proteins_with_hits
                  << "\tcandidates=" << writer.candidates_total
                  << "\tkept=" << writer.kept_total << '\n';
    }

    return 0;
}
