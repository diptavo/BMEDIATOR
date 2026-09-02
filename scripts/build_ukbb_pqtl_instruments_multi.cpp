#include "../plink_ld.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using bmediator::PlinkData;

struct ProteinInfo {
    std::string protein;
    std::string gene;
    std::string sumstat_dir;
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
    while (std::getline(ss, item, '\t')) out.push_back(trim(item));
    return out;
}

std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream iss(line);
    std::string item;
    while (iss >> item) out.push_back(item);
    return out;
}

std::string uppercase(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::string chr_token(int chr) {
    if (chr >= 1 && chr <= 22) return std::to_string(chr);
    if (chr == 23) return "X";
    if (chr == 24) return "Y";
    return "XY";
}

std::vector<ProteinInfo> read_resource(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("cannot open resource file: " + path);
    std::string line;
    if (!std::getline(fin, line)) throw std::runtime_error("empty resource file: " + path);
    auto header = split_tab(line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;
    const std::vector<std::string> required = {"protein", "gene", "sumstat_dir", "chr", "tss"};
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
        p.sumstat_dir = fields[col["sumstat_dir"]];
        if (!parse_int(fields[col["chr"]], p.chr)) continue;
        if (!parse_int(fields[col["tss"]], p.tss)) continue;
        if (p.protein.empty() || p.sumstat_dir.empty()) continue;
        out.push_back(p);
    }
    return out;
}

std::unordered_map<std::string, std::string> load_variant_map_chr(const std::string& variant_dir, int chr) {
    std::unordered_map<std::string, std::string> out;
    const std::string path = variant_dir + "/olink_rsid_map_mac5_info03_b0_7_chr" + chr_token(chr) + "_patched_v2.tsv.gz";
    const std::string cmd = "gzip -dc '" + path + "'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) throw std::runtime_error("cannot open variant map: " + path);

    char buffer[1 << 16];
    if (fgets(buffer, sizeof(buffer), pipe) == nullptr) {
        pclose(pipe);
        return out;
    }

    std::string line;
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        line = buffer;
        if (!line.empty() && line.back() == '\n') line.pop_back();
        if (trim(line).empty()) continue;
        auto fields = split_tab(line);
        if (fields.size() < 4) continue;
        out.emplace(fields[0], fields[3]);
    }
    pclose(pipe);
    return out;
}

std::string find_sumstat_file(const std::string& dir, int chr) {
    const std::string token = "_chr" + chr_token(chr) + "_";
    const std::string cmd = "find '" + dir + "' -maxdepth 1 -type f | grep '" + token + "' | head -n 1";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) throw std::runtime_error("cannot inspect UKBB directory: " + dir);
    char buffer[1 << 16];
    std::string line;
    if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        line = buffer;
        if (!line.empty() && line.back() == '\n') line.pop_back();
    }
    pclose(pipe);
    if (line.empty()) throw std::runtime_error("no chromosome file found for chr" + chr_token(chr) + " in " + dir);
    return line;
}

std::vector<CandidateHit> read_candidates_union(const ProteinInfo& protein,
                                                const PlinkData& plink,
                                                const std::string& variant_dir,
                                                std::unordered_map<int, std::unordered_map<std::string, std::string>>& variant_maps,
                                                double max_p_thresh,
                                                int flank_kb) {
    if (!variant_maps.count(protein.chr)) {
        variant_maps.emplace(protein.chr, load_variant_map_chr(variant_dir, protein.chr));
    }
    const auto& variant_map = variant_maps.at(protein.chr);

    const std::string sumstat_file = find_sumstat_file(protein.sumstat_dir, protein.chr);
    const std::string cmd = "gzip -dc '" + sumstat_file + "'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) throw std::runtime_error("cannot open UKBB sumstat: " + sumstat_file);

    char buffer[1 << 16];
    if (fgets(buffer, sizeof(buffer), pipe) == nullptr) {
        pclose(pipe);
        return {};
    }
    std::string line = buffer;
    if (!line.empty() && line.back() == '\n') line.pop_back();
    auto header = split_ws(line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;

    const std::vector<std::string> required = {
        "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "LOG10P"
    };
    for (const auto& key : required) {
        if (!col.count(key)) throw std::runtime_error("missing UKBB column " + key + " in " + sumstat_file);
    }

    const int flank_bp = flank_kb * 1000;
    std::vector<CandidateHit> out;
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        line = buffer;
        if (!line.empty() && line.back() == '\n') line.pop_back();
        if (trim(line).empty()) continue;
        auto fields = split_ws(line);
        if (static_cast<int>(fields.size()) <= col["LOG10P"]) continue;
        if (fields[col["TEST"]] != "ADD") continue;

        CandidateHit hit;
        if (!parse_int(fields[col["CHROM"]], hit.chr)) continue;
        if (hit.chr != protein.chr) continue;
        if (!parse_int(fields[col["GENPOS"]], hit.bp)) continue;
        if (hit.bp < protein.tss - flank_bp || hit.bp > protein.tss + flank_bp) continue;
        double log10p = 0.0;
        if (!parse_double(fields[col["LOG10P"]], log10p)) continue;
        hit.p = std::pow(10.0, -log10p);
        if (hit.p >= max_p_thresh) continue;
        if (!parse_double(fields[col["BETA"]], hit.beta)) continue;
        if (!parse_double(fields[col["SE"]], hit.se)) continue;
        if (!parse_double(fields[col["A1FREQ"]], hit.freq)) continue;
        if (!parse_double(fields[col["N"]], hit.n)) continue;

        hit.a1 = uppercase(fields[col["ALLELE1"]]);
        hit.a2 = uppercase(fields[col["ALLELE0"]]);
        const std::string raw_id = fields[col["ID"]];
        auto map_it = variant_map.find(raw_id);
        if (map_it == variant_map.end()) continue;
        hit.rsid = map_it->second;
        auto it = plink.rsid_to_idx.find(hit.rsid);
        if (it == plink.rsid_to_idx.end()) continue;
        hit.bim_idx = it->second;
        out.push_back(hit);
    }

    pclose(pipe);
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

std::vector<ComboConfig> combos(const std::string& prefix) {
    return {
        {5e-6, 0.01, prefix + "_cis_tss_pm6_r2_0.01"},
        {5e-6, 0.1, prefix + "_cis_tss_pm6_r2_0.1"},
        {5e-8, 0.01, prefix + "_cis_tss_pm8_r2_0.01"},
        {5e-8, 0.1, prefix + "_cis_tss_pm8_r2_0.1"},
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
    if (argc != 6) {
        std::cerr << "Usage: build_ukbb_pqtl_instruments_multi <resource_subset.tsv> <variant_dir> <plink_prefix> <out_dir> <chunk_label>\n";
        return 1;
    }

    const std::string resource_path = argv[1];
    const std::string variant_dir = argv[2];
    const std::string plink_prefix = argv[3];
    const std::string out_dir = argv[4];
    const std::string chunk_label = argv[5];

    std::string prefix = "ukbb";
    if (resource_path.find("EUR") != std::string::npos || resource_path.find("eur") != std::string::npos) {
        prefix = "ukbb_eur_oid";
    } else if (resource_path.find("COMBINED") != std::string::npos || resource_path.find("combined") != std::string::npos) {
        prefix = "ukbb_combined_oid";
    }

    auto proteins = read_resource(resource_path);
    PlinkData plink;
    if (!plink.load(plink_prefix)) {
        std::cerr << "failed to load PLINK reference: " << plink_prefix << "\n";
        return 1;
    }

    std::unordered_map<int, std::unordered_map<std::string, std::string>> variant_maps;
    std::vector<ComboWriter> writers;
    for (const auto& cfg : combos(prefix)) {
        ComboWriter w;
        w.cfg = cfg;
        open_writer(w, out_dir, chunk_label);
        writers.push_back(std::move(w));
    }

    for (int i = 0; i < static_cast<int>(proteins.size()); ++i) {
        const auto& protein = proteins[i];
        std::vector<CandidateHit> union_hits;
        try {
            union_hits = read_candidates_union(protein, plink, variant_dir, variant_maps, 5e-6, 500);
        } catch (const std::exception& e) {
            for (auto& writer : writers) {
                ++writer.proteins_seen;
                writer.summary_out << protein.protein << '\t' << protein.gene << "\tread_error\t0\t0\t"
                                   << protein.sumstat_dir << '\n';
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
                                   << filtered.size() << "\t0\t" << protein.sumstat_dir << '\n';
                continue;
            }
            ++writer.proteins_with_hits;
            writer.kept_total += static_cast<long long>(kept.size());
            writer.info_out << protein.protein << '\t' << protein.gene << '\t' << protein.chr
                            << '\t' << protein.tss << '\t' << protein.tss << '\n';
            writer.summary_out << protein.protein << '\t' << protein.gene << "\tok\t"
                               << filtered.size() << '\t' << kept.size() << '\t'
                               << protein.sumstat_dir << '\n';
            for (const auto& hit : kept) {
                writer.pqtl_out << protein.protein << '\t' << hit.rsid << '\t' << hit.a1 << '\t'
                                << hit.a2 << '\t' << std::setprecision(10) << hit.freq << '\t'
                                << hit.beta << '\t' << hit.se << '\t' << hit.p << '\t'
                                << hit.n << '\t' << hit.chr << '\t' << hit.bp << '\n';
            }
        }

        if ((i + 1) % 25 == 0) {
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
