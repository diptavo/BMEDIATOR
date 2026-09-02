#include "../plink_ld.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
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
    std::string seqid;
    std::string gene;
    std::string sumstat_file;
    int chr = 0;
    int gene_start = 0;
    int gene_end = 0;
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
    int chr = 0;
    int bp = 0;
    int bim_idx = -1;
};

struct ComboConfig {
    std::string cis_mode;
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
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, '\t')) {
        item = trim(item);
        if (item.size() >= 2 && item.front() == '"' && item.back() == '"') {
            item = item.substr(1, item.size() - 2);
        }
        parts.push_back(item);
    }
    return parts;
}

std::string uppercase(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::string> split_rsids(const std::string& raw) {
    std::vector<std::string> out;
    std::string token;
    auto flush = [&]() {
        token = trim(token);
        if (!token.empty() && token != "NA" && token != ".") out.push_back(token);
        token.clear();
    };
    for (char c : raw) {
        if (c == ';' || c == ',' || c == '|' || std::isspace(static_cast<unsigned char>(c))) {
            flush();
        } else {
            token.push_back(c);
        }
    }
    flush();
    return out;
}

int parse_chr(std::string s) {
    s = trim(s);
    if (s.rfind("chr", 0) == 0 || s.rfind("CHR", 0) == 0) s = s.substr(3);
    if (s == "X") return 23;
    if (s == "Y") return 24;
    if (s == "M" || s == "MT") return 25;
    int chr = 0;
    if (!parse_int(s, chr)) return 0;
    return chr;
}

std::vector<ProteinInfo> read_protein_info(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("cannot open protein resource file: " + path);
    }

    std::string line;
    if (!std::getline(fin, line)) {
        throw std::runtime_error("protein resource file is empty: " + path);
    }

    auto header = split_tab(line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;

    const std::vector<std::string> required = {
        "seqid", "symbol", "sumstat_file", "chromosome_hg38", "gene_start_hg38",
        "gene_end_hg38", "tss_hg38"
    };
    for (const auto& key : required) {
        if (!col.count(key)) {
            throw std::runtime_error("missing column in protein resource file: " + key);
        }
    }

    std::vector<ProteinInfo> out;
    while (std::getline(fin, line)) {
        if (trim(line).empty()) continue;
        auto fields = split_tab(line);
        if (static_cast<int>(fields.size()) <= col["tss_hg38"]) continue;

        ProteinInfo info;
        info.seqid = fields[col["seqid"]];
        info.gene = fields[col["symbol"]];
        info.sumstat_file = fields[col["sumstat_file"]];

        if (info.seqid.empty() || info.sumstat_file.empty()) continue;
        if (!parse_int(fields[col["chromosome_hg38"]], info.chr)) continue;
        if (!parse_int(fields[col["gene_start_hg38"]], info.gene_start)) continue;
        if (!parse_int(fields[col["gene_end_hg38"]], info.gene_end)) continue;
        if (!parse_int(fields[col["tss_hg38"]], info.tss)) continue;

        out.push_back(info);
    }
    return out;
}

bool in_cis_window(const ProteinInfo& protein, int variant_chr, int variant_bp,
                   const std::string& cis_mode, int flank_kb) {
    if (variant_chr != protein.chr) return false;
    const int flank_bp = flank_kb * 1000;
    if (cis_mode == "tss") {
        return variant_bp >= protein.tss - flank_bp && variant_bp <= protein.tss + flank_bp;
    }
    if (cis_mode == "genebody") {
        return variant_bp >= protein.gene_start - flank_bp &&
               variant_bp <= protein.gene_end + flank_bp;
    }
    throw std::runtime_error("unsupported cis mode: " + cis_mode);
}

bool in_union_window(const ProteinInfo& protein, int variant_chr, int variant_bp, int flank_kb) {
    if (variant_chr != protein.chr) return false;
    const int flank_bp = flank_kb * 1000;
    const int lo = std::min(protein.tss, protein.gene_start) - flank_bp;
    const int hi = std::max(protein.tss, protein.gene_end) + flank_bp;
    return variant_bp >= lo && variant_bp <= hi;
}

std::vector<CandidateHit> read_candidates_union(const ProteinInfo& protein,
                                                const PlinkData& plink,
                                                double max_p_thresh,
                                                int flank_kb) {
    const std::string cmd = "gzip -dc '" + protein.sumstat_file + "'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) {
        throw std::runtime_error("cannot open gzip stream for " + protein.sumstat_file);
    }

    char buffer[1 << 16];
    std::string line;
    if (fgets(buffer, sizeof(buffer), pipe) == nullptr) {
        pclose(pipe);
        return {};
    }
    line = buffer;
    if (!line.empty() && line.back() == '\n') line.pop_back();

    auto header = split_tab(line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;

    const std::vector<std::string> required = {
        "Chrom", "Pos", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "SE", "ImpMAF"
    };
    for (const auto& key : required) {
        if (!col.count(key)) {
            pclose(pipe);
            throw std::runtime_error("missing column in deCODE sumstat file " + protein.sumstat_file +
                                     ": " + key);
        }
    }

    std::vector<CandidateHit> hits;
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        line = buffer;
        if (!line.empty() && line.back() == '\n') line.pop_back();
        if (trim(line).empty()) continue;

        auto fields = split_tab(line);
        if (static_cast<int>(fields.size()) <= col["ImpMAF"]) continue;

        CandidateHit hit;
        if (!parse_double(fields[col["Pval"]], hit.p)) continue;
        if (hit.p >= max_p_thresh) continue;
        if (!parse_double(fields[col["Beta"]], hit.beta)) continue;
        if (!parse_double(fields[col["SE"]], hit.se)) continue;
        if (!parse_double(fields[col["ImpMAF"]], hit.freq)) continue;
        if (!parse_int(fields[col["Pos"]], hit.bp)) continue;

        hit.chr = parse_chr(fields[col["Chrom"]]);
        if (hit.chr == 0) continue;
        if (!in_union_window(protein, hit.chr, hit.bp, flank_kb)) continue;

        hit.a1 = uppercase(fields[col["effectAllele"]]);
        hit.a2 = uppercase(fields[col["otherAllele"]]);
        if (hit.a1.empty() || hit.a2.empty()) continue;

        const auto rsid_tokens = split_rsids(fields[col["rsids"]]);
        int chosen_idx = -1;
        std::string chosen_rsid;
        for (const auto& rsid : rsid_tokens) {
            auto it = plink.rsid_to_idx.find(rsid);
            if (it != plink.rsid_to_idx.end()) {
                chosen_rsid = rsid;
                chosen_idx = it->second;
                break;
            }
        }
        if (chosen_idx < 0) continue;

        hit.rsid = chosen_rsid;
        hit.bim_idx = chosen_idx;
        hits.push_back(hit);
    }

    pclose(pipe);
    return hits;
}

std::vector<CandidateHit> filter_candidates(const std::vector<CandidateHit>& union_hits,
                                            const ProteinInfo& protein,
                                            const std::string& cis_mode,
                                            double p_thresh,
                                            int flank_kb) {
    std::vector<CandidateHit> out;
    out.reserve(union_hits.size());
    for (const auto& hit : union_hits) {
        if (hit.p >= p_thresh) continue;
        if (!in_cis_window(protein, hit.chr, hit.bp, cis_mode, flank_kb)) continue;
        out.push_back(hit);
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
        double s = 0.0;
        double s2 = 0.0;
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

std::vector<ComboConfig> build_combos() {
    return {
        {"tss", 5e-6, 0.01, "decode_seqid_cis_tss_pm6_r2_0.01"},
        {"tss", 5e-6, 0.1, "decode_seqid_cis_tss_pm6_r2_0.1"},
        {"tss", 5e-8, 0.01, "decode_seqid_cis_tss_pm8_r2_0.01"},
        {"tss", 5e-8, 0.1, "decode_seqid_cis_tss_pm8_r2_0.1"},
        {"genebody", 5e-6, 0.01, "decode_seqid_cis_genebody_pm6_r2_0.01"},
        {"genebody", 5e-6, 0.1, "decode_seqid_cis_genebody_pm6_r2_0.1"},
        {"genebody", 5e-8, 0.01, "decode_seqid_cis_genebody_pm8_r2_0.01"},
        {"genebody", 5e-8, 0.1, "decode_seqid_cis_genebody_pm8_r2_0.1"},
    };
}

void open_writer(ComboWriter& writer, const std::string& out_dir, const std::string& chunk_label) {
    const std::string base = out_dir + "/" + writer.cfg.tag + "." + chunk_label;
    writer.pqtl_out.open(base + ".tsv");
    writer.info_out.open(base + ".protein_info.tsv");
    writer.summary_out.open(base + ".build_summary.tsv");
    if (!writer.pqtl_out.is_open() || !writer.info_out.is_open() || !writer.summary_out.is_open()) {
        throw std::runtime_error("failed to open output bundle for " + base);
    }
    writer.pqtl_out << "PROTEIN\tSNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tCHR\tBP\n";
    writer.info_out << "PROTEIN\tGENE\tCHR\tSTART\tEND\n";
    writer.summary_out << "PROTEIN\tGENE\tSTATUS\tN_CANDIDATES\tN_KEPT\tSOURCE\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: build_decode_pqtl_instruments_multi <decode_resource_subset.tsv> "
                     "<plink_prefix> <out_dir> <chunk_label>\n";
        return 1;
    }

    const std::string resource_path = argv[1];
    const std::string plink_prefix = argv[2];
    const std::string out_dir = argv[3];
    const std::string chunk_label = argv[4];
    const int flank_kb = 500;
    const double max_p_thresh = 5e-6;

    const auto proteins = read_protein_info(resource_path);

    PlinkData plink;
    if (!plink.load(plink_prefix)) {
        std::cerr << "failed to load PLINK reference: " << plink_prefix << "\n";
        return 1;
    }

    std::vector<ComboWriter> writers;
    for (const auto& cfg : build_combos()) {
        ComboWriter writer;
        writer.cfg = cfg;
        open_writer(writer, out_dir, chunk_label);
        writers.push_back(std::move(writer));
    }

    int proteins_seen = 0;
    for (const auto& protein : proteins) {
        ++proteins_seen;
        std::vector<CandidateHit> union_hits;
        try {
            union_hits = read_candidates_union(protein, plink, max_p_thresh, flank_kb);
        } catch (const std::exception& e) {
            for (auto& writer : writers) {
                ++writer.proteins_seen;
                writer.summary_out << protein.seqid << '\t' << protein.gene << "\tread_error\t0\t0\t"
                                   << protein.sumstat_file << '\n';
            }
            std::cerr << chunk_label << ": " << e.what() << "\n";
            continue;
        }

        for (auto& writer : writers) {
            ++writer.proteins_seen;
            const auto filtered = filter_candidates(
                union_hits, protein, writer.cfg.cis_mode, writer.cfg.p_thresh, flank_kb);
            writer.candidates_total += static_cast<long long>(filtered.size());

            const auto kept = ld_clump_hits(plink, filtered, writer.cfg.r2_thresh, 10000);
            if (kept.empty()) {
                writer.summary_out << protein.seqid << '\t' << protein.gene << "\tno_hits\t"
                                   << filtered.size() << "\t0\t" << protein.sumstat_file << '\n';
                continue;
            }

            ++writer.proteins_with_hits;
            writer.kept_total += static_cast<long long>(kept.size());
            writer.info_out << protein.seqid << '\t' << protein.gene << '\t' << protein.chr << '\t'
                            << protein.gene_start << '\t' << protein.gene_end << '\n';
            writer.summary_out << protein.seqid << '\t' << protein.gene << "\tok\t"
                               << filtered.size() << '\t' << kept.size() << '\t'
                               << protein.sumstat_file << '\n';
            for (const auto& hit : kept) {
                writer.pqtl_out << protein.seqid << '\t' << hit.rsid << '\t' << hit.a1 << '\t'
                                << hit.a2 << '\t' << std::setprecision(10) << hit.freq << '\t'
                                << hit.beta << '\t' << hit.se << '\t' << hit.p << '\t'
                                << hit.chr << '\t' << hit.bp << '\n';
            }
        }

        if (proteins_seen % 100 == 0) {
            std::cerr << chunk_label << ": processed " << proteins_seen
                      << " proteins from subset\n";
        }
    }

    for (const auto& writer : writers) {
        std::cerr << chunk_label << '\t' << writer.cfg.tag
                  << "\tproteins_scanned=" << writer.proteins_seen
                  << "\tproteins_with_hits=" << writer.proteins_with_hits
                  << "\tcandidates=" << writer.candidates_total
                  << "\tkept=" << writer.kept_total << "\n";
    }

    return 0;
}
