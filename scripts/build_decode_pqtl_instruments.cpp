#include "../plink_ld.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
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
    if (s.rfind("chr", 0) == 0 || s.rfind("CHR", 0) == 0) {
        s = s.substr(3);
    }
    if (s == "X") return 23;
    if (s == "Y") return 24;
    if (s == "M" || s == "MT") return 25;
    int chr = 0;
    if (!parse_int(s, chr)) return 0;
    return chr;
}

std::map<std::string, ProteinInfo> read_protein_info(const std::string& path) {
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

    std::map<std::string, ProteinInfo> out;
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

        out[info.seqid] = info;
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

std::vector<CandidateHit> read_candidates(const ProteinInfo& protein,
                                          const PlinkData& plink,
                                          double p_thresh,
                                          const std::string& cis_mode,
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
        if (hit.p >= p_thresh) continue;
        if (!parse_double(fields[col["Beta"]], hit.beta)) continue;
        if (!parse_double(fields[col["SE"]], hit.se)) continue;
        if (!parse_double(fields[col["ImpMAF"]], hit.freq)) continue;
        if (!parse_int(fields[col["Pos"]], hit.bp)) continue;

        hit.chr = parse_chr(fields[col["Chrom"]]);
        if (hit.chr == 0) continue;
        if (!in_cis_window(protein, hit.chr, hit.bp, cis_mode, flank_kb)) continue;

        hit.a1 = uppercase(fields[col["effectAllele"]]);
        hit.a2 = uppercase(fields[col["otherAllele"]]);
        if (hit.a1.empty() || hit.a2.empty()) continue;

        auto rsid_tokens = split_rsids(fields[col["rsids"]]);
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

std::vector<CandidateHit> ld_clump_hits(const PlinkData& plink,
                                        const std::vector<CandidateHit>& hits,
                                        double r2_thresh,
                                        int window_kb) {
    if (hits.empty()) return {};

    std::vector<int> order(hits.size());
    for (int i = 0; i < static_cast<int>(hits.size()); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return hits[a].p < hits[b].p;
    });

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
        double var = s2 / geno[i].size() - means[i] * means[i];
        sds[i] = (var > 1e-10) ? std::sqrt(var) : 0.0;
    }

    const int window_bp = window_kb * 1000;
    std::vector<bool> removed(hits.size(), false);
    std::vector<CandidateHit> kept;
    kept.reserve(hits.size());

    for (int oi = 0; oi < static_cast<int>(order.size()); ++oi) {
        int idx = order[oi];
        if (removed[idx]) continue;
        kept.push_back(hits[idx]);

        for (int oj = oi + 1; oj < static_cast<int>(order.size()); ++oj) {
            int jdx = order[oj];
            if (removed[jdx]) continue;
            if (hits[idx].chr != hits[jdx].chr) continue;
            if (std::abs(hits[idx].bp - hits[jdx].bp) > window_bp) continue;
            if (sds[idx] < 1e-10 || sds[jdx] < 1e-10) continue;

            double cov = 0.0;
            for (int s = 0; s < plink.n_samples; ++s) {
                cov += (geno[idx][s] - means[idx]) * (geno[jdx][s] - means[jdx]);
            }
            cov /= plink.n_samples;
            double r = cov / (sds[idx] * sds[jdx]);
            if (r * r > r2_thresh) {
                removed[jdx] = true;
            }
        }
    }

    return kept;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 9) {
        std::cerr << "Usage: build_decode_pqtl_instruments <decode_resource.tsv> <plink_prefix> "
                     "<cis_mode:tss|genebody> <out_pqtl.tsv> <out_protein_info.tsv> "
                     "<out_summary.tsv> <p_thresh> <r2_thresh>\n";
        return 1;
    }

    const std::string resource_path = argv[1];
    const std::string plink_prefix = argv[2];
    const std::string cis_mode = argv[3];
    const std::string out_pqtl = argv[4];
    const std::string out_info = argv[5];
    const std::string out_summary = argv[6];
    const double p_thresh = std::stod(argv[7]);
    const double r2_thresh = std::stod(argv[8]);
    const int flank_kb = 500;

    if (cis_mode != "tss" && cis_mode != "genebody") {
        std::cerr << "cis_mode must be one of: tss, genebody\n";
        return 1;
    }

    auto proteins = read_protein_info(resource_path);

    PlinkData plink;
    if (!plink.load(plink_prefix)) {
        std::cerr << "failed to load PLINK reference: " << plink_prefix << "\n";
        return 1;
    }

    std::ofstream pqtl_out(out_pqtl);
    std::ofstream info_out(out_info);
    std::ofstream summary_out(out_summary);
    if (!pqtl_out.is_open() || !info_out.is_open() || !summary_out.is_open()) {
        std::cerr << "failed to open output files\n";
        return 1;
    }

    pqtl_out << "PROTEIN\tSNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tCHR\tBP\n";
    info_out << "PROTEIN\tGENE\tCHR\tSTART\tEND\n";
    summary_out << "PROTEIN\tGENE\tSTATUS\tN_CANDIDATES\tN_KEPT\tSOURCE\n";

    int proteins_seen = 0;
    int proteins_with_hits = 0;
    long long candidates_total = 0;
    long long kept_total = 0;

    for (const auto& kv : proteins) {
        const auto& protein = kv.second;
        ++proteins_seen;

        std::vector<CandidateHit> candidates;
        try {
            candidates = read_candidates(protein, plink, p_thresh, cis_mode, flank_kb);
        } catch (const std::exception& e) {
            summary_out << protein.seqid << '\t' << protein.gene << "\tread_error\t0\t0\t"
                        << protein.sumstat_file << '\n';
            std::cerr << e.what() << "\n";
            continue;
        }

        candidates_total += static_cast<long long>(candidates.size());
        auto kept = ld_clump_hits(plink, candidates, r2_thresh, 10000);
        if (kept.empty()) {
            summary_out << protein.seqid << '\t' << protein.gene << "\tno_hits\t"
                        << candidates.size() << "\t0\t" << protein.sumstat_file << '\n';
            continue;
        }

        ++proteins_with_hits;
        kept_total += static_cast<long long>(kept.size());
        info_out << protein.seqid << '\t' << protein.gene << '\t' << protein.chr << '\t'
                 << protein.tss << '\t' << protein.tss << '\n';
        summary_out << protein.seqid << '\t' << protein.gene << "\tok\t"
                    << candidates.size() << '\t' << kept.size() << '\t'
                    << protein.sumstat_file << '\n';

        for (const auto& hit : kept) {
            pqtl_out << protein.seqid << '\t' << hit.rsid << '\t' << hit.a1 << '\t' << hit.a2
                     << '\t' << std::setprecision(10) << hit.freq << '\t'
                     << hit.beta << '\t' << hit.se << '\t' << hit.p << '\t'
                     << hit.chr << '\t' << hit.bp << '\n';
        }

        if (proteins_seen % 250 == 0) {
            std::cerr << "processed " << proteins_seen << " proteins, kept "
                      << kept_total << " instruments so far\n";
        }
    }

    std::cerr << "proteins scanned: " << proteins_seen << "\n";
    std::cerr << "proteins with instruments: " << proteins_with_hits << "\n";
    std::cerr << "candidate hits: " << candidates_total << "\n";
    std::cerr << "kept instruments after clumping: " << kept_total << "\n";

    return 0;
}
