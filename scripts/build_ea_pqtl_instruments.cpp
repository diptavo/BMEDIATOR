#include "../plink_ld.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

using bmediator::PlinkData;

struct SeqInfo {
    std::string seqid;
    std::string gene;
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

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, '\t')) {
        parts.push_back(item);
    }
    return parts;
}

std::map<std::string, SeqInfo> read_seqinfo(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("cannot open seqid mapping: " + path);
    }

    std::string line;
    std::getline(fin, line);
    auto header = split_tab(line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) {
        col[header[i]] = i;
    }

    const std::vector<std::string> required = {
        "seqid_in_sample", "entrezgenesymbol", "chromosome_name", "transcription_start_site"
    };
    for (const auto& key : required) {
        if (!col.count(key)) {
            throw std::runtime_error("missing column in seqid mapping: " + key);
        }
    }

    std::map<std::string, SeqInfo> out;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        auto fields = split_tab(line);
        if (static_cast<int>(fields.size()) <= col["transcription_start_site"]) continue;

        SeqInfo info;
        info.seqid = fields[col["seqid_in_sample"]];
        info.gene = fields[col["entrezgenesymbol"]];
        if (!parse_int(fields[col["chromosome_name"]], info.chr)) continue;
        if (!parse_int(fields[col["transcription_start_site"]], info.tss)) continue;
        out[info.seqid] = info;
    }
    return out;
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

std::vector<CandidateHit> read_candidates(const std::string& path,
                                          const PlinkData& plink,
                                          double p_thresh) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("cannot open cis summary file: " + path);
    }

    std::string line;
    std::getline(fin, line);
    auto header = split_tab(line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) {
        col[header[i]] = i;
    }

    const std::vector<std::string> required = {
        "#CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ", "TEST", "BETA", "SE", "P"
    };
    for (const auto& key : required) {
        if (!col.count(key)) {
            throw std::runtime_error("missing column in cis summary file: " + key);
        }
    }

    std::vector<CandidateHit> hits;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        auto fields = split_tab(line);
        if (static_cast<int>(fields.size()) <= col["P"]) continue;
        if (fields[col["TEST"]] != "ADD") continue;
        if (fields[col["ID"]] == "." || fields[col["ERRCODE"]] != ".") continue;

        CandidateHit hit;
        hit.rsid = fields[col["ID"]];
        hit.a1 = fields[col["A1"]];
        if (!derive_other_allele(fields[col["REF"]], fields[col["ALT"]], hit.a1, hit.a2)) {
            continue;
        }
        if (!parse_double(fields[col["A1_FREQ"]], hit.freq)) continue;
        if (!parse_double(fields[col["BETA"]], hit.beta)) continue;
        if (!parse_double(fields[col["SE"]], hit.se)) continue;
        if (!parse_double(fields[col["P"]], hit.p)) continue;
        if (!parse_int(fields[col["#CHROM"]], hit.chr)) continue;
        if (!parse_int(fields[col["POS"]], hit.bp)) continue;
        if (hit.p >= p_thresh) continue;

        auto it = plink.rsid_to_idx.find(hit.rsid);
        if (it == plink.rsid_to_idx.end()) continue;
        hit.bim_idx = it->second;
        hits.push_back(hit);
    }
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
    if (argc != 8) {
        std::cerr << "Usage: build_ea_pqtl_instruments <seqid.txt> <cis_dir> <plink_prefix> "
                     "<out_pqtl.tsv> <out_protein_info.tsv> <p_thresh> <r2_thresh>\n";
        return 1;
    }

    const std::string seq_path = argv[1];
    const std::string cis_dir = argv[2];
    const std::string plink_prefix = argv[3];
    const std::string out_pqtl = argv[4];
    const std::string out_info = argv[5];
    const double p_thresh = std::stod(argv[6]);
    const double r2_thresh = std::stod(argv[7]);

    auto seqinfo = read_seqinfo(seq_path);

    PlinkData plink;
    if (!plink.load(plink_prefix)) {
        std::cerr << "failed to load PLINK reference: " << plink_prefix << "\n";
        return 1;
    }

    std::ofstream pqtl_out(out_pqtl);
    std::ofstream info_out(out_info);
    if (!pqtl_out.is_open() || !info_out.is_open()) {
        std::cerr << "failed to open output files\n";
        return 1;
    }

    pqtl_out << "PROTEIN\tSNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tCHR\tBP\n";
    info_out << "PROTEIN\tGENE\tCHR\tSTART\tEND\n";

    int proteins_seen = 0;
    int proteins_with_hits = 0;
    long long candidates_total = 0;
    long long kept_total = 0;

    for (const auto& kv : seqinfo) {
        const auto& seqid = kv.first;
        const auto& meta = kv.second;
        ++proteins_seen;

        const std::string cis_file = cis_dir + "/" + seqid + ".PHENO1.glm.linear";
        std::ifstream test(cis_file);
        if (!test.is_open()) continue;
        test.close();

        auto candidates = read_candidates(cis_file, plink, p_thresh);
        candidates_total += static_cast<long long>(candidates.size());
        auto kept = ld_clump_hits(plink, candidates, r2_thresh, 10000);
        if (kept.empty()) continue;

        ++proteins_with_hits;
        kept_total += static_cast<long long>(kept.size());
        info_out << seqid << '\t' << meta.gene << '\t' << meta.chr << '\t'
                 << meta.tss << '\t' << meta.tss << '\n';

        for (const auto& hit : kept) {
            pqtl_out << seqid << '\t' << hit.rsid << '\t' << hit.a1 << '\t' << hit.a2 << '\t'
                     << std::setprecision(10) << hit.freq << '\t'
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
