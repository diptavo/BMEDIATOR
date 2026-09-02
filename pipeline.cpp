#include "bmediator.h"
#include "plink_ld.h"
#include "gsmr_qc.h"

#include <cassert>

namespace bmediator {

namespace {

struct PqtlData {
    std::map<std::string, std::map<std::string, SumStat>> ss_by_protein;
    std::map<std::string, std::map<std::string, double>> pval_by_protein;
};

struct ProteinFileSpec {
    std::string protein_id;
    std::string path;
};

struct RfInstrumentRef {
    std::string rsid;
    int chr = 0;
    int bp = 0;
    int bim_idx = -1;
};

struct ProxyMatch {
    bool found = false;
    std::string rsid;
    int bim_idx = -1;
    double r = 0.0;
    double r2 = 0.0;
};

struct CisVariantRef {
    std::string rsid;
    int chr = 0;
    int bp = 0;
    int bim_idx = -1;
    double pval = 1.0;
};

static std::unordered_set<std::string> collect_needed_outcome_rsids(
        const std::map<std::string, double>& rf_pval,
        const PqtlData& pqtl,
        double rf_p_thresh) {
    std::unordered_set<std::string> keep_rsids;
    for (const auto& kv : rf_pval) {
        if (kv.second < rf_p_thresh) {
            keep_rsids.insert(kv.first);
        }
    }
    for (const auto& prot_kv : pqtl.ss_by_protein) {
        for (const auto& snp_kv : prot_kv.second) {
            keep_rsids.insert(snp_kv.first);
        }
    }
    return keep_rsids;
}

static void read_needed_outcome_sumstats(const std::string& file,
                                         const std::unordered_set<std::string>& keep_rsids,
                                         std::map<std::string, SumStat>& cancer_ss,
                                         std::map<std::string, double>& cancer_pval,
                                         const Options& opts) {
    std::cout << "Reading outcome summary statistics for "
              << keep_rsids.size() << " needed SNPs only...\n";
    read_sumstats_with_pval_subset(file, cancer_ss, cancer_pval, keep_rsids,
                                   true, 0.01, opts.verbose);
}

static bool harmonize_snp(SumStat& s, const BimEntry& ref) {
    HarmonizeResult hr = harmonize_alleles(s.a1, s.a2, ref.a1, ref.a2, s.beta);
    if (hr.action == 0) return false;
    s.beta = hr.beta_adj;
    if (hr.action == -1) {
        std::swap(s.a1, s.a2);
        s.freq = 1.0 - s.freq;
    }
    return true;
}

static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static std::string basename_no_ext(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
    const std::string suffix = ".PHENO1.glm.linear";
    if (name.size() > suffix.size() &&
        name.substr(name.size() - suffix.size()) == suffix) {
        return name.substr(0, name.size() - suffix.size());
    }
    size_t dot = name.find('.');
    return (dot == std::string::npos) ? name : name.substr(0, dot);
}

static std::vector<ProteinFileSpec> read_protein_manifest(const std::string& file) {
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open protein GWAS manifest " << file << "\n";
        std::exit(1);
    }

    std::vector<ProteinFileSpec> specs;
    std::string line;
    while (std::getline(fin, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) tokens.push_back(token);
        if (tokens.empty()) continue;
        ProteinFileSpec spec;
        if (tokens.size() == 1) {
            spec.path = tokens[0];
            spec.protein_id = basename_no_ext(spec.path);
        } else {
            spec.protein_id = tokens[0];
            spec.path = tokens[1];
        }
        specs.push_back(spec);
    }
    return specs;
}

static bool derive_other_allele(const std::string& ref, const std::string& alt,
                                const std::string& a1, std::string& a2) {
    if (a1 == alt) { a2 = ref; return true; }
    if (a1 == ref) { a2 = alt; return true; }
    return false;
}

static void read_single_protein_gwas_file(const ProteinFileSpec& spec,
                                          PqtlData& out,
                                          const Options& opts) {
    std::ifstream fin(spec.path);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open protein GWAS file " << spec.path << "\n";
        std::exit(1);
    }

    std::string line;
    if (!std::getline(fin, line)) return;
    std::vector<std::string> header;
    {
        std::istringstream iss(line);
        std::string field;
        while (std::getline(iss, field, '\t')) header.push_back(field);
    }
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;

    const std::vector<std::string> required = {
        "#CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ", "TEST", "BETA", "SE", "P"
    };
    for (const auto& key : required) {
        if (!col.count(key)) {
            std::cerr << "Error: missing column '" << key << "' in " << spec.path << "\n";
            std::exit(1);
        }
    }

    int n = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields;
        {
            std::istringstream iss(line);
            std::string field;
            while (std::getline(iss, field, '\t')) fields.push_back(field);
        }
        if (static_cast<int>(fields.size()) <= col["P"]) continue;
        if (fields[col["TEST"]] != "ADD") continue;
        if (col.count("ERRCODE") && fields[col["ERRCODE"]] != ".") continue;
        if (fields[col["ID"]] == ".") continue;

        SumStat s;
        double pval;
        s.rsid = fields[col["ID"]];
        s.a1 = fields[col["A1"]];
        if (!derive_other_allele(fields[col["REF"]], fields[col["ALT"]], s.a1, s.a2)) {
            continue;
        }
        try {
            s.freq = std::stod(fields[col["A1_FREQ"]]);
            s.beta = std::stod(fields[col["BETA"]]);
            s.se = std::stod(fields[col["SE"]]);
            pval = std::stod(fields[col["P"]]);
            s.chr = std::stoi(fields[col["#CHROM"]]);
            s.bp = std::stoi(fields[col["POS"]]);
        } catch (...) {
            continue;
        }
        for (auto& c : s.a1) c = toupper(c);
        for (auto& c : s.a2) c = toupper(c);
        out.ss_by_protein[spec.protein_id][s.rsid] = s;
        out.pval_by_protein[spec.protein_id][s.rsid] = pval;
        n++;
    }
    write_log("  Loaded " + std::to_string(n) + " variants for protein " + spec.protein_id, opts);
}

static PqtlData read_protein_gwas_manifest(const std::string& manifest_file,
                                           const Options& opts) {
    PqtlData data;
    auto specs = read_protein_manifest(manifest_file);
    std::cout << "  Reading protein GWAS manifest with " << specs.size() << " entries\n";
    int done = 0;
    for (const auto& spec : specs) {
        read_single_protein_gwas_file(spec, data, opts);
        done++;
        if (opts.verbose && done % 100 == 0) {
            std::cout << "    processed " << done << " / " << specs.size()
                      << " protein GWAS files\n";
        }
    }
    return data;
}

static void harmonize_map_against_reference(std::map<std::string, SumStat>& ss,
                                            const PlinkData& plink,
                                            const std::string& label) {
    int n_harmonized = 0, n_excluded = 0;
    for (auto it = ss.begin(); it != ss.end(); ) {
        auto ref_it = plink.rsid_to_idx.find(it->first);
        if (ref_it == plink.rsid_to_idx.end()) {
            it = ss.erase(it);
            continue;
        }
        if (!harmonize_snp(it->second, plink.bim[ref_it->second])) {
            n_excluded++;
            it = ss.erase(it);
            continue;
        }
        n_harmonized++;
        ++it;
    }
    std::cout << "  " << label << ": " << n_harmonized
              << " harmonized, " << n_excluded << " excluded\n";
}

static void harmonize_pqtl_data(PqtlData& pqtl, const PlinkData& plink, const Options& opts) {
    int proteins_done = 0;
    for (auto& prot_kv : pqtl.ss_by_protein) {
        for (auto it = prot_kv.second.begin(); it != prot_kv.second.end(); ) {
            auto ref_it = plink.rsid_to_idx.find(it->first);
            if (ref_it == plink.rsid_to_idx.end()) {
                it = prot_kv.second.erase(it);
                continue;
            }
            if (!harmonize_snp(it->second, plink.bim[ref_it->second])) {
                it = prot_kv.second.erase(it);
                continue;
            }
            ++it;
        }
        proteins_done++;
        if (opts.verbose && proteins_done % 250 == 0) {
            std::cout << "  Harmonized pQTL data for " << proteins_done << " proteins\n";
        }
    }
}

static QCParams build_qc_params(const Options& opts) {
    QCParams qc;
    qc.maf_min = opts.maf_thresh;
    qc.freq_diff_thresh = opts.freq_diff_thresh;
    qc.clump_window_kb = opts.clump_window_kb;
    qc.ld_r2_thresh = opts.r2_thresh;
    qc.heidi_thresh = opts.heidi_thresh;
    qc.heidi_flag = opts.heidi_flag;
    qc.heidi_global = opts.heidi_global;
    qc.steiger_flag = opts.steiger_flag;
    qc.remove_palindromic = opts.remove_palindromic;
    qc.ld_fdr_thresh = opts.ld_fdr_thresh;
    qc.min_instruments = opts.min_instruments;
    return qc;
}

static std::map<int, std::vector<RfInstrumentRef>> build_rf_index(
    const std::set<std::string>& rf_instruments,
    const std::map<std::string, SumStat>& rf_ss,
    const PlinkData& plink) {
    std::map<int, std::vector<RfInstrumentRef>> by_chr;
    for (const auto& rsid : rf_instruments) {
        auto rf_it = rf_ss.find(rsid);
        auto ref_it = plink.rsid_to_idx.find(rsid);
        if (rf_it == rf_ss.end() || ref_it == plink.rsid_to_idx.end()) continue;
        by_chr[rf_it->second.chr].push_back({rsid, rf_it->second.chr, rf_it->second.bp, ref_it->second});
    }
    for (auto& kv : by_chr) {
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const RfInstrumentRef& a, const RfInstrumentRef& b) { return a.bp < b.bp; });
    }
    return by_chr;
}

static std::vector<RfInstrumentRef> nearby_rf_instruments(
    const std::map<int, std::vector<RfInstrumentRef>>& rf_by_chr,
    int chr, int bp, int window_bp) {
    std::vector<RfInstrumentRef> out;
    auto it = rf_by_chr.find(chr);
    if (it == rf_by_chr.end()) return out;
    for (const auto& rf : it->second) {
        if (rf.bp < bp - window_bp) continue;
        if (rf.bp > bp + window_bp) break;
        out.push_back(rf);
    }
    return out;
}

static double ld_r2_with_cache(const PlinkData& plink,
                               std::map<std::pair<int, int>, double>& cache,
                               int a, int b) {
    if (a > b) std::swap(a, b);
    std::pair<int, int> key(a, b);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    double r2 = plink.compute_ld_r2(a, b);
    cache[key] = r2;
    return r2;
}

static double ld_r_with_cache(const PlinkData& plink,
                              std::map<std::pair<int, int>, double>& cache,
                              int a, int b) {
    if (a > b) std::swap(a, b);
    std::pair<int, int> key(a, b);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    double r = plink.compute_ld_r(a, b);
    cache[key] = r;
    return r;
}

static ProxyMatch best_rf_proxy_for_cis(const PlinkData& plink,
                                        std::map<std::pair<int, int>, double>& r_cache,
                                        const std::map<int, std::vector<RfInstrumentRef>>& rf_by_chr,
                                        const SumStat& s,
                                        double r2_thresh,
                                        int window_bp) {
    ProxyMatch best;
    auto self_it = plink.rsid_to_idx.find(s.rsid);
    if (self_it == plink.rsid_to_idx.end()) return best;
    auto neighbors = nearby_rf_instruments(rf_by_chr, s.chr, s.bp, window_bp);
    for (const auto& rf : neighbors) {
        double r = ld_r_with_cache(plink, r_cache, self_it->second, rf.bim_idx);
        double r2 = r * r;
        if (r2 < r2_thresh) continue;
        if (!best.found || r2 > best.r2) {
            best.found = true;
            best.r = r;
            best.r2 = r2;
            best.rsid = rf.rsid;
            best.bim_idx = rf.bim_idx;
        }
    }
    return best;
}

static ProxyMatch best_cis_proxy_for_rf(const PlinkData& plink,
                                        std::map<std::pair<int, int>, double>& r_cache,
                                        const std::vector<CisVariantRef>& cis_variants,
                                        const RfInstrumentRef& rf,
                                        double r2_thresh,
                                        int window_bp) {
    ProxyMatch best;
    for (const auto& cis : cis_variants) {
        if (cis.chr != rf.chr) continue;
        if (std::abs(cis.bp - rf.bp) > window_bp) continue;
        double r = ld_r_with_cache(plink, r_cache, cis.bim_idx, rf.bim_idx);
        double r2 = r * r;
        if (r2 < r2_thresh) continue;
        if (!best.found || r2 > best.r2) {
            best.found = true;
            best.r = r;
            best.r2 = r2;
            best.rsid = cis.rsid;
            best.bim_idx = cis.bim_idx;
        }
    }
    return best;
}

static std::vector<double> invert_matrix_diag(std::vector<std::vector<double>> a) {
    int n = static_cast<int>(a.size());
    std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) inv[i][i] = 1.0;
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int r = i + 1; r < n; ++r) {
            if (std::fabs(a[r][i]) > std::fabs(a[pivot][i])) pivot = r;
        }
        if (std::fabs(a[pivot][i]) < 1e-8) {
            a[i][i] += 1e-4;
            pivot = i;
        }
        if (pivot != i) {
            std::swap(a[pivot], a[i]);
            std::swap(inv[pivot], inv[i]);
        }
        double diag = a[i][i];
        if (std::fabs(diag) < 1e-8) diag = (diag >= 0.0 ? 1e-8 : -1e-8);
        for (int c = 0; c < n; ++c) {
            a[i][c] /= diag;
            inv[i][c] /= diag;
        }
        for (int r = 0; r < n; ++r) {
            if (r == i) continue;
            double factor = a[r][i];
            if (std::fabs(factor) < 1e-12) continue;
            for (int c = 0; c < n; ++c) {
                a[r][c] -= factor * a[i][c];
                inv[r][c] -= factor * inv[i][c];
            }
        }
    }
    std::vector<double> diag(n, 1.0);
    for (int i = 0; i < n; ++i) diag[i] = std::max(inv[i][i], 1e-6);
    return diag;
}

static void assign_ld_weights(ProteinData& prot,
                              const PlinkData& plink,
                              std::map<std::pair<int, int>, double>& cache,
                              int block_max_size = 64) {
    int nAC = prot.nA() + prot.nC();
    int nU = nAC + prot.nB();
    prot.ld_weight_alpha_ac.assign(nAC, 1.0);
    prot.ld_weight_cancer_union.assign(nU, 1.0);

    struct UnionSnp { int chr; int bp; int bim_idx; };
    std::vector<UnionSnp> ac_snps;
    std::vector<UnionSnp> union_snps;
    ac_snps.reserve(nAC);
    union_snps.reserve(nU);

    for (const auto& rsid : prot.setA_rsid) {
        int bim = plink.rsid_to_idx.at(rsid);
        const auto& b = plink.bim[bim];
        ac_snps.push_back({b.chr, b.bp, bim});
        union_snps.push_back({b.chr, b.bp, bim});
    }
    for (const auto& rsid : prot.setC_rsid) {
        int bim = plink.rsid_to_idx.at(rsid);
        const auto& b = plink.bim[bim];
        ac_snps.push_back({b.chr, b.bp, bim});
        union_snps.push_back({b.chr, b.bp, bim});
    }
    for (const auto& rsid : prot.setB_rsid) {
        int bim = plink.rsid_to_idx.at(rsid);
        const auto& b = plink.bim[bim];
        union_snps.push_back({b.chr, b.bp, bim});
    }

    std::vector<double> alpha_sum(nAC, 0.0), union_sum(nU, 0.0);
    const int window_bp = 1000 * 1000;

    for (int i = 0; i < nAC; ++i) {
        for (int j = i + 1; j < nAC; ++j) {
            if (ac_snps[i].chr != ac_snps[j].chr) continue;
            if (std::abs(ac_snps[i].bp - ac_snps[j].bp) > window_bp) continue;
            double r2 = ld_r2_with_cache(plink, cache, ac_snps[i].bim_idx, ac_snps[j].bim_idx);
            alpha_sum[i] += r2;
            alpha_sum[j] += r2;
        }
    }
    for (int i = 0; i < nU; ++i) {
        for (int j = i + 1; j < nU; ++j) {
            if (union_snps[i].chr != union_snps[j].chr) continue;
            if (std::abs(union_snps[i].bp - union_snps[j].bp) > window_bp) continue;
            double r2 = ld_r2_with_cache(plink, cache, union_snps[i].bim_idx, union_snps[j].bim_idx);
            union_sum[i] += r2;
            union_sum[j] += r2;
        }
    }
    auto rescale_weights = [](const std::vector<double>& corr_sum) {
        std::vector<double> w(corr_sum.size(), 1.0);
        if (corr_sum.empty()) return w;
        double raw_sum = 0.0;
        double pair_sum = 0.0;
        for (double x : corr_sum) pair_sum += x;
        pair_sum *= 0.5;
        for (size_t i = 0; i < corr_sum.size(); ++i) {
            // Stronger downweighting than 1/(1+sum r^2), but preserve total effective information.
            w[i] = 1.0 / std::sqrt(1.0 + corr_sum[i]);
            raw_sum += w[i];
        }
        double n = static_cast<double>(corr_sum.size());
        double n_eff = (n * n) / std::max(n + 2.0 * pair_sum, 1.0);
        double scale = (raw_sum > 0.0) ? (n_eff / raw_sum) : 1.0;
        for (double& wi : w) wi *= scale;
        return w;
    };
    prot.ld_weight_alpha_ac = rescale_weights(alpha_sum);
    prot.ld_weight_cancer_union = rescale_weights(union_sum);

    auto refine_with_local_inverse = [&](const std::vector<UnionSnp>& snps,
                                         std::vector<double>& weights) {
        if (snps.empty()) return;
        int start = 0;
        while (start < (int)snps.size()) {
            int end = start + 1;
            while (end < (int)snps.size() &&
                   snps[end].chr == snps[end - 1].chr &&
                   std::abs(snps[end].bp - snps[end - 1].bp) <= window_bp &&
                   (end - start) < block_max_size) {
                end++;
            }
            int block_n = end - start;
            if (block_n >= 2) {
                std::vector<int> block_idx;
                block_idx.reserve(block_n);
                for (int i = start; i < end; ++i) block_idx.push_back(snps[i].bim_idx);
                auto ld = plink.compute_ld_matrix(block_idx);
                for (int i = 0; i < block_n; ++i) ld[i][i] += 0.05;
                auto inv_diag = invert_matrix_diag(ld);
                double raw_sum = 0.0;
                for (double x : inv_diag) raw_sum += x;
                if (raw_sum > 0.0) {
                    double target = 0.0;
                    for (int i = start; i < end; ++i) target += weights[i];
                    double scale = target / raw_sum;
                    for (int i = start; i < end; ++i) {
                        weights[i] = std::max(1e-4, inv_diag[i - start] * scale);
                    }
                }
            }
            start = end;
        }
    };
    refine_with_local_inverse(ac_snps, prot.ld_weight_alpha_ac);
    refine_with_local_inverse(union_snps, prot.ld_weight_cancer_union);
}

static std::set<std::string> select_rf_instruments(const std::map<std::string, SumStat>& rf_ss,
                                                   const std::map<std::string, double>& rf_pval,
                                                   const std::map<std::string, SumStat>& cancer_ss,
                                                   const std::map<std::string, double>& cancer_pval,
                                                   const PlinkData& plink,
                                                   const QCParams& qc,
                                                   const Options& opts) {
    std::cout << "\n--- Pre-MR Quality Control ---\n";
    std::cout << "QC parameters:\n";
    std::cout << "  MAF filter:               >= " << qc.maf_min << "\n";
    std::cout << "  Freq concordance:         |diff| < " << qc.freq_diff_thresh << "\n";
    std::cout << "  Palindromic SNP removal:  " << (qc.remove_palindromic ? "YES" : "NO") << "\n";
    std::cout << "  HEIDI-outlier:            " << (qc.heidi_flag ? (qc.heidi_global ? "global (GSMR2)" : "single-SNP") : "OFF") << "\n";
    std::cout << "  HEIDI p threshold:        " << qc.heidi_thresh << "\n";
    std::cout << "  Steiger filter:           " << (qc.steiger_flag ? "YES" : "NO") << "\n";
    std::cout << "  LD clumping r2:           " << qc.ld_r2_thresh << "\n";
    std::cout << "  LD FDR (chance LD):       " << qc.ld_fdr_thresh << "\n";
    std::cout << "  Min instruments:          " << qc.min_instruments << "\n\n";

    QCReport rf_qc_report;
    std::cout << "Running SNP QC on RF GWAS...\n";
    std::vector<std::string> rf_qc_snps = run_snp_qc(rf_ss, rf_pval, plink, qc, rf_qc_report);
    std::cout << "  Input: " << rf_qc_report.n_input << " -> after MAF: " << rf_qc_report.n_after_maf
              << " -> after freq: " << rf_qc_report.n_after_freq_check
              << " -> after palindromic: " << rf_qc_report.n_after_palindromic
              << " -> after SE/beta: " << rf_qc_report.n_after_se_check
              << " -> QC pass: " << rf_qc_report.n_final << "\n";

    QCReport ca_qc_report;
    std::cout << "Running SNP QC on cancer GWAS...\n";
    std::vector<std::string> ca_qc_snps = run_snp_qc(cancer_ss, cancer_pval, plink, qc, ca_qc_report);
    std::cout << "  QC pass: " << ca_qc_report.n_final << "\n";

    std::set<std::string> ca_qc_set(ca_qc_snps.begin(), ca_qc_snps.end());
    std::vector<int> rf_candidate_bim;
    std::vector<double> rf_candidate_pval;
    std::vector<std::string> rf_candidate_rsid;
    for (const auto& rsid : rf_qc_snps) {
        auto p_it = rf_pval.find(rsid);
        auto ref_it = plink.rsid_to_idx.find(rsid);
        if (p_it == rf_pval.end() || ref_it == plink.rsid_to_idx.end()) continue;
        if (p_it->second >= opts.p_thresh_rf) continue;
        if (!ca_qc_set.count(rsid)) continue;
        rf_candidate_bim.push_back(ref_it->second);
        rf_candidate_pval.push_back(p_it->second);
        rf_candidate_rsid.push_back(rsid);
    }

    std::cout << "\nSelecting RF instruments (p < " << opts.p_thresh_rf
              << ", LD r2 < " << qc.ld_r2_thresh << ")...\n";
    std::cout << "  " << rf_candidate_bim.size()
              << " RF SNPs pass p-value threshold + QC\n";

    ClumpResult rf_clump = ld_clump(plink, rf_candidate_bim, rf_candidate_pval,
                                    qc.ld_r2_thresh, qc.clump_window_kb);

    std::vector<std::string> rf_clumped_rsid;
    std::vector<int> rf_clumped_bim;
    for (int idx : rf_clump.index_snps) {
        rf_clumped_rsid.push_back(rf_candidate_rsid[idx]);
        rf_clumped_bim.push_back(rf_candidate_bim[idx]);
    }
    std::cout << "  " << rf_clumped_rsid.size()
              << " independent RF instruments after LD clumping\n";

    if (rf_clumped_bim.size() > 1) {
        auto ld_mat = plink.compute_ld_matrix(rf_clumped_bim);
        auto keep_ld = filter_chance_ld(ld_mat, plink.n_samples, qc.ld_fdr_thresh);
        std::vector<std::string> rf_after_ld;
        int n_removed = 0;
        for (size_t i = 0; i < rf_clumped_rsid.size(); ++i) {
            if (keep_ld[i]) rf_after_ld.push_back(rf_clumped_rsid[i]);
            else n_removed++;
        }
        if (n_removed > 0) {
            std::cout << "  " << n_removed << " instruments removed by chance LD filter\n";
        }
        rf_clumped_rsid.swap(rf_after_ld);
    }

    std::set<std::string> rf_instruments(rf_clumped_rsid.begin(), rf_clumped_rsid.end());
    std::cout << "  " << rf_instruments.size() << " RF instruments after all QC\n";
    return rf_instruments;
}

static void copy_mediation_outputs(const Options& opts,
                                   std::vector<ProteinData>& proteins) {
    std::cout << "\nRunning Bayesian Mediation MR...\n";
    Hyperparams hyp;
    hyp.p0 = opts.prior_p0;
    hyp.p1 = opts.prior_p1;
    hyp.p2 = opts.prior_p2;
    hyp.p3 = opts.prior_p3;
    hyp.p4 = opts.prior_p4;
    hyp.p5 = opts.prior_p5;
    hyp.sigma2_beta1 = opts.prior_sigma2_beta1;
    hyp.sigma2_beta2 = opts.prior_sigma2_beta2;
    hyp.sigma2_beta3 = opts.prior_sigma2_beta3;

    std::vector<ProteinResult> results;
    run_empirical_bayes(proteins, hyp, results, opts);

    std::sort(results.begin(), results.end(),
              [](const ProteinResult& a, const ProteinResult& b) {
                  if (a.selection_probability != b.selection_probability) {
                      return a.selection_probability > b.selection_probability;
                  }
                  return a.prob_M1 > b.prob_M1;
              });

    write_results(results, hyp, opts);
}

static void write_instrument_log(const std::vector<ProteinData>& proteins,
                                 const Options& opts) {
    std::string inst_file = opts.out_prefix + ".instruments";
    std::ofstream iout(inst_file);
    iout << "Protein\tSet\tSNP\tGamma\tAlpha\tGamma_cancer\n";
    for (const auto& prot : proteins) {
        for (int k = 0; k < prot.nA(); ++k) {
            iout << prot.protein_id << "\tA\t" << prot.setA_rsid[k] << "\t"
                 << prot.setA_gamma[k] << "\t" << prot.setA_alpha[k] << "\t"
                 << prot.setA_Gamma[k] << "\n";
        }
        for (int l = 0; l < prot.nB(); ++l) {
            iout << prot.protein_id << "\tB\t" << prot.setB_rsid[l] << "\tNA\t"
                 << prot.setB_alpha_cis[l] << "\t" << prot.setB_Gamma_cis[l] << "\n";
        }
        for (int c = 0; c < prot.nC(); ++c) {
            iout << prot.protein_id << "\tC\t" << prot.setC_rsid[c] << "\t"
                 << prot.setC_gamma[c] << "\t" << prot.setC_alpha[c] << "\t"
                 << prot.setC_Gamma[c] << "\n";
        }
    }
    iout.close();
    std::cout << "Instrument details saved to " << inst_file << "\n";
}

static void apply_set_qc(ProteinData& prot,
                         const PlinkData& plink,
                         int n_ref,
                         const QCParams& qc,
                         double n_pqtl,
                         double n_cancer,
                         int& total_heidi_removed,
                         const Options& opts) {
    if (qc.heidi_flag && prot.nA() >= qc.heidi_min_snps) {
        double bxy_total = ivw_estimate(prot.setA_gamma, prot.setA_Gamma, prot.setA_se_Gamma);
        std::vector<int> setA_bim;
        for (const auto& rsid : prot.setA_rsid) setA_bim.push_back(plink.rsid_to_idx.at(rsid));
        auto ld_A = plink.compute_ld_matrix(setA_bim);
        auto keep_heidi = heidi_outlier_test(prot.setA_gamma, prot.setA_se_gamma,
                                             prot.setA_Gamma, prot.setA_se_Gamma,
                                             ld_A, n_ref, bxy_total,
                                             qc.heidi_thresh, qc.heidi_global);
        int n_before = prot.nA();
        std::vector<std::string> new_rsid;
        std::vector<double> new_gamma, new_se_gamma, new_alpha, new_se_alpha, new_alpha_rel, new_Gamma, new_se_Gamma;
        std::vector<bool> new_alpha_observed;
        for (int k = 0; k < n_before; ++k) {
            if (!keep_heidi[k]) continue;
            new_rsid.push_back(prot.setA_rsid[k]);
            new_gamma.push_back(prot.setA_gamma[k]);
            new_se_gamma.push_back(prot.setA_se_gamma[k]);
            new_alpha.push_back(prot.setA_alpha[k]);
            new_se_alpha.push_back(prot.setA_se_alpha[k]);
            new_alpha_rel.push_back(prot.setA_alpha_reliability[k]);
            new_alpha_observed.push_back(prot.setA_alpha_observed[k]);
            new_Gamma.push_back(prot.setA_Gamma[k]);
            new_se_Gamma.push_back(prot.setA_se_Gamma[k]);
        }
        total_heidi_removed += n_before - static_cast<int>(new_rsid.size());
        prot.setA_rsid = new_rsid;
        prot.setA_gamma = new_gamma;
        prot.setA_se_gamma = new_se_gamma;
        prot.setA_alpha = new_alpha;
        prot.setA_se_alpha = new_se_alpha;
        prot.setA_alpha_reliability = new_alpha_rel;
        prot.setA_alpha_observed = new_alpha_observed;
        prot.setA_Gamma = new_Gamma;
        prot.setA_se_Gamma = new_se_Gamma;
    }

    if (qc.heidi_flag && prot.nB() >= qc.heidi_min_snps) {
        double bxy_cis = ivw_estimate(prot.setB_alpha_cis, prot.setB_Gamma_cis, prot.setB_se_Gamma_cis);
        std::vector<int> setB_bim;
        for (const auto& rsid : prot.setB_rsid) setB_bim.push_back(plink.rsid_to_idx.at(rsid));
        auto ld_B = plink.compute_ld_matrix(setB_bim);
        auto keep_heidi_B = heidi_outlier_test(prot.setB_alpha_cis, prot.setB_se_alpha_cis,
                                               prot.setB_Gamma_cis, prot.setB_se_Gamma_cis,
                                               ld_B, n_ref, bxy_cis,
                                               qc.heidi_thresh, qc.heidi_global);
        int n_before = prot.nB();
        std::vector<std::string> new_rsid;
        std::vector<double> new_alpha, new_se_alpha, new_Gamma, new_se_Gamma;
        for (int l = 0; l < n_before; ++l) {
            if (!keep_heidi_B[l]) continue;
            new_rsid.push_back(prot.setB_rsid[l]);
            new_alpha.push_back(prot.setB_alpha_cis[l]);
            new_se_alpha.push_back(prot.setB_se_alpha_cis[l]);
            new_Gamma.push_back(prot.setB_Gamma_cis[l]);
            new_se_Gamma.push_back(prot.setB_se_Gamma_cis[l]);
        }
        total_heidi_removed += n_before - static_cast<int>(new_rsid.size());
        prot.setB_rsid = new_rsid;
        prot.setB_alpha_cis = new_alpha;
        prot.setB_se_alpha_cis = new_se_alpha;
        prot.setB_Gamma_cis = new_Gamma;
        prot.setB_se_Gamma_cis = new_se_Gamma;
    }

    if (qc.steiger_flag && prot.nB() >= 3) {
        bool correct_dir = steiger_direction_test(prot.setB_alpha_cis, prot.setB_se_alpha_cis, n_pqtl,
                                                  prot.setB_Gamma_cis, prot.setB_se_Gamma_cis, n_cancer);
        if (!correct_dir) {
            write_log("    Steiger cleared Set B for " + prot.protein_id, opts);
            prot.setB_rsid.clear();
            prot.setB_alpha_cis.clear();
            prot.setB_se_alpha_cis.clear();
            prot.setB_Gamma_cis.clear();
            prot.setB_se_Gamma_cis.clear();
        }
    }

    prot.ld_weight_alpha_ac.clear();
    prot.ld_weight_cancer_union.clear();
    std::map<std::pair<int, int>, double> qc_ld_cache;
    assign_ld_weights(prot, plink, qc_ld_cache, opts.ld_block_max_size);
}

static void build_legacy_sets(std::vector<ProteinData>& proteins,
                              const std::map<std::string, SumStat>& rf_ss,
                              const std::map<std::string, SumStat>& cancer_ss,
                              const PqtlData& pqtl,
                              const std::set<std::string>& rf_instruments,
                              const PlinkData& plink,
                              const Options& opts) {
    std::cout << "\nBuilding instrument sets...\n";
    int n_with = 0;
    int logged = 0;
    int cis_window_bp = opts.cis_window_kb * 1000;
    int overlap_window_bp = opts.clump_window_kb * 1000;
    auto rf_by_chr = build_rf_index(rf_instruments, rf_ss, plink);
    std::map<std::pair<int, int>, double> ld_r2_cache;
    std::map<std::pair<int, int>, double> ld_r_cache;
    for (auto& prot : proteins) {
        prot.ld_reference_used = true;
        auto pq_it = pqtl.ss_by_protein.find(prot.protein_id);
        auto pp_it = pqtl.pval_by_protein.find(prot.protein_id);
        if (pq_it == pqtl.ss_by_protein.end() || pp_it == pqtl.pval_by_protein.end()) continue;

        std::set<std::string> cis_instruments;
        std::vector<CisVariantRef> cis_all;
        for (const auto& kv : pp_it->second) {
            auto sit = pq_it->second.find(kv.first);
            if (sit == pq_it->second.end()) continue;
            const SumStat& s = sit->second;
            if (s.chr == prot.gene_chr &&
                s.bp >= prot.gene_start - cis_window_bp &&
                s.bp <= prot.gene_end + cis_window_bp) {
                auto ref_it = plink.rsid_to_idx.find(kv.first);
                if (ref_it == plink.rsid_to_idx.end()) continue;
                cis_all.push_back({kv.first, s.chr, s.bp, ref_it->second, kv.second});
                if (kv.second < opts.p_thresh_cis) cis_instruments.insert(kv.first);
            }
        }

        std::set<std::string> used_cis;
        std::set<std::string> consumed_rf;
        for (const auto& rsid : cis_instruments) {
            if (!rf_instruments.count(rsid)) continue;
            auto rf_it = rf_ss.find(rsid);
            auto pq_snp = pq_it->second.find(rsid);
            auto ca_it = cancer_ss.find(rsid);
            if (rf_it == rf_ss.end() || pq_snp == pq_it->second.end() || ca_it == cancer_ss.end()) continue;
            prot.setC_rsid.push_back(rsid);
            prot.setC_gamma.push_back(rf_it->second.beta);
            prot.setC_se_gamma.push_back(rf_it->second.se);
            prot.setC_alpha.push_back(pq_snp->second.beta);
            prot.setC_se_alpha.push_back(pq_snp->second.se);
            prot.setC_alpha_reliability.push_back(1.0);
            prot.setC_Gamma.push_back(ca_it->second.beta);
            prot.setC_se_Gamma.push_back(ca_it->second.se);
            prot.nC_exact++;
            used_cis.insert(rsid);
            consumed_rf.insert(rsid);
        }

        for (const auto& rsid : cis_instruments) {
            auto pq_snp = pq_it->second.find(rsid);
            auto ca_it = cancer_ss.find(rsid);
            if (used_cis.count(rsid) || pq_snp == pq_it->second.end() || ca_it == cancer_ss.end()) continue;
            ProxyMatch proxy = best_rf_proxy_for_cis(plink, ld_r_cache, rf_by_chr, pq_snp->second,
                                                     opts.proxy_r2_thresh, overlap_window_bp);
            if (!proxy.found || consumed_rf.count(proxy.rsid)) continue;
            auto rf_it = rf_ss.find(proxy.rsid);
            if (rf_it == rf_ss.end()) continue;
            prot.setC_rsid.push_back(rsid);
            prot.setC_gamma.push_back(proxy.r * rf_it->second.beta);
            prot.setC_se_gamma.push_back(rf_it->second.se);
            prot.setC_alpha.push_back(pq_snp->second.beta);
            prot.setC_se_alpha.push_back(pq_snp->second.se);
            prot.setC_alpha_reliability.push_back(std::fabs(proxy.r));
            prot.setC_Gamma.push_back(ca_it->second.beta);
            prot.setC_se_Gamma.push_back(ca_it->second.se);
            prot.nC_proxy++;
            used_cis.insert(rsid);
            consumed_rf.insert(proxy.rsid);
        }

        for (const auto& rsid : rf_instruments) {
            if (consumed_rf.count(rsid)) continue;
            auto rf_it = rf_ss.find(rsid);
            auto ca_it = cancer_ss.find(rsid);
            if (rf_it == rf_ss.end() || ca_it == cancer_ss.end()) continue;
            auto pq_snp = pq_it->second.find(rsid);
            double alpha_val = std::numeric_limits<double>::quiet_NaN();
            double se_alpha_val = std::numeric_limits<double>::quiet_NaN();
            bool alpha_observed = false;
            double alpha_rel = 0.0;
            if (pq_snp != pq_it->second.end()) {
                alpha_val = pq_snp->second.beta;
                se_alpha_val = pq_snp->second.se;
                alpha_observed = true;
                alpha_rel = 1.0;
            } else {
                auto ref_it = plink.rsid_to_idx.find(rsid);
                if (ref_it != plink.rsid_to_idx.end()) {
                    ProxyMatch proxy = best_cis_proxy_for_rf(plink, ld_r_cache, cis_all,
                                                             {rsid, rf_it->second.chr, rf_it->second.bp, ref_it->second},
                                                             opts.proxy_r2_thresh, overlap_window_bp);
                    if (proxy.found) {
                        auto proxy_it = pq_it->second.find(proxy.rsid);
                        if (proxy_it != pq_it->second.end()) {
                            alpha_val = proxy.r * proxy_it->second.beta;
                            se_alpha_val = proxy_it->second.se;
                            alpha_observed = true;
                            alpha_rel = std::fabs(proxy.r);
                            prot.nA_proxy++;
                        }
                    }
                }
            }
            prot.setA_rsid.push_back(rsid);
            prot.setA_gamma.push_back(rf_it->second.beta);
            prot.setA_se_gamma.push_back(rf_it->second.se);
            prot.setA_alpha.push_back(alpha_val);
            prot.setA_se_alpha.push_back(se_alpha_val);
            prot.setA_alpha_observed.push_back(alpha_observed);
            prot.setA_alpha_reliability.push_back(alpha_rel);
            prot.setA_Gamma.push_back(ca_it->second.beta);
            prot.setA_se_Gamma.push_back(ca_it->second.se);
        }

        for (const auto& rsid : cis_instruments) {
            if (used_cis.count(rsid)) continue;
            auto pq_snp = pq_it->second.find(rsid);
            auto ca_it = cancer_ss.find(rsid);
            if (pq_snp == pq_it->second.end() || ca_it == cancer_ss.end()) continue;
            prot.setB_rsid.push_back(rsid);
            prot.setB_alpha_cis.push_back(pq_snp->second.beta);
            prot.setB_se_alpha_cis.push_back(pq_snp->second.se);
            prot.setB_Gamma_cis.push_back(ca_it->second.beta);
            prot.setB_se_Gamma_cis.push_back(ca_it->second.se);
        }

        assign_ld_weights(prot, plink, ld_r2_cache, opts.ld_block_max_size);

        if (prot.nTotal() > 0) n_with++;
        if (opts.verbose && logged < 15) {
            std::cout << "  Protein " << prot.protein_id << ": nA=" << prot.nA()
                      << ", nB=" << prot.nB() << ", nC=" << prot.nC() << "\n";
            logged++;
        }
    }
    std::cout << "  " << n_with << " / " << proteins.size() << " proteins have instruments.\n";
}

static void build_full_sets(std::vector<ProteinData>& proteins,
                            const std::map<std::string, SumStat>& rf_ss,
                            const std::map<std::string, SumStat>& cancer_ss,
                            const PqtlData& pqtl,
                            const std::set<std::string>& rf_instruments,
                            const PlinkData& plink,
                            const QCParams& qc,
                            const Options& opts) {
    std::cout << "\nBuilding instrument sets for " << proteins.size()
              << " proteins (raw protein GWAS + cis clumping)...\n";
    int cis_window_bp = opts.cis_window_kb * 1000;
    int n_with_instruments = 0;
    int total_heidi_removed = 0;
    int processed = 0;
    int overlap_window_bp = opts.clump_window_kb * 1000;
    auto rf_by_chr = build_rf_index(rf_instruments, rf_ss, plink);
    std::map<std::pair<int, int>, double> ld_r2_cache;
    std::map<std::pair<int, int>, double> ld_r_cache;

    for (auto& prot : proteins) {
        prot.ld_reference_used = true;
        auto pq_it = pqtl.ss_by_protein.find(prot.protein_id);
        auto pp_it = pqtl.pval_by_protein.find(prot.protein_id);
        if (pq_it == pqtl.ss_by_protein.end() || pp_it == pqtl.pval_by_protein.end()) continue;

        std::vector<int> cis_candidate_bim;
        std::vector<double> cis_candidate_pval;
        std::vector<std::string> cis_candidate_rsid;
        std::vector<CisVariantRef> cis_all;
        for (const auto& snp_kv : pq_it->second) {
            const std::string& rsid = snp_kv.first;
            const SumStat& s = snp_kv.second;
            if (s.chr != prot.gene_chr) continue;
            if (s.bp < prot.gene_start - cis_window_bp || s.bp > prot.gene_end + cis_window_bp) continue;
            auto pval_it = pp_it->second.find(rsid);
            auto ref_it = plink.rsid_to_idx.find(rsid);
            if (pval_it == pp_it->second.end() || ref_it == plink.rsid_to_idx.end()) continue;
            cis_all.push_back({rsid, s.chr, s.bp, ref_it->second, pval_it->second});
            auto outcome_it = cancer_ss.find(rsid);
            if (outcome_it != cancer_ss.end() &&
                std::isfinite(s.beta) && std::isfinite(s.se) && s.se > 0.0 &&
                std::isfinite(outcome_it->second.beta) &&
                std::isfinite(outcome_it->second.se) && outcome_it->second.se > 0.0) {
                prot.regional_cis_rsid.push_back(rsid);
                prot.regional_pp_beta.push_back(s.beta);
                prot.regional_pp_se.push_back(s.se);
                prot.regional_outcome_beta.push_back(outcome_it->second.beta);
                prot.regional_outcome_se.push_back(outcome_it->second.se);
            }
            if (pval_it->second >= opts.p_thresh_cis) continue;
            cis_candidate_bim.push_back(ref_it->second);
            cis_candidate_pval.push_back(pval_it->second);
            cis_candidate_rsid.push_back(rsid);
        }

        std::set<std::string> cis_instruments;
        if (!cis_candidate_bim.empty()) {
            ClumpResult cis_clump = ld_clump(plink, cis_candidate_bim, cis_candidate_pval,
                                             opts.r2_thresh, cis_window_bp / 1000);
            for (int idx : cis_clump.index_snps) cis_instruments.insert(cis_candidate_rsid[idx]);
        }

        std::set<std::string> used_cis;
        std::set<std::string> consumed_rf;
        for (const auto& rsid : cis_instruments) {
            if (!rf_instruments.count(rsid)) continue;
            auto rf_it = rf_ss.find(rsid);
            auto pq_snp = pq_it->second.find(rsid);
            auto ca_it = cancer_ss.find(rsid);
            if (rf_it == rf_ss.end() || pq_snp == pq_it->second.end() || ca_it == cancer_ss.end()) continue;
            prot.setC_rsid.push_back(rsid);
            prot.setC_gamma.push_back(rf_it->second.beta);
            prot.setC_se_gamma.push_back(rf_it->second.se);
            prot.setC_alpha.push_back(pq_snp->second.beta);
            prot.setC_se_alpha.push_back(pq_snp->second.se);
            prot.setC_alpha_reliability.push_back(1.0);
            prot.setC_Gamma.push_back(ca_it->second.beta);
            prot.setC_se_Gamma.push_back(ca_it->second.se);
            prot.nC_exact++;
            used_cis.insert(rsid);
            consumed_rf.insert(rsid);
        }

        for (const auto& rsid : cis_instruments) {
            auto pq_snp = pq_it->second.find(rsid);
            auto ca_it = cancer_ss.find(rsid);
            if (used_cis.count(rsid) || pq_snp == pq_it->second.end() || ca_it == cancer_ss.end()) continue;
            ProxyMatch proxy = best_rf_proxy_for_cis(plink, ld_r_cache, rf_by_chr, pq_snp->second,
                                                     opts.proxy_r2_thresh, overlap_window_bp);
            if (!proxy.found || consumed_rf.count(proxy.rsid)) continue;
            auto rf_it = rf_ss.find(proxy.rsid);
            if (rf_it == rf_ss.end()) continue;
            prot.setC_rsid.push_back(rsid);
            prot.setC_gamma.push_back(proxy.r * rf_it->second.beta);
            prot.setC_se_gamma.push_back(rf_it->second.se);
            prot.setC_alpha.push_back(pq_snp->second.beta);
            prot.setC_se_alpha.push_back(pq_snp->second.se);
            prot.setC_alpha_reliability.push_back(std::fabs(proxy.r));
            prot.setC_Gamma.push_back(ca_it->second.beta);
            prot.setC_se_Gamma.push_back(ca_it->second.se);
            prot.nC_proxy++;
            used_cis.insert(rsid);
            consumed_rf.insert(proxy.rsid);
        }

        for (const auto& rsid : rf_instruments) {
            if (consumed_rf.count(rsid)) continue;
            auto rf_it = rf_ss.find(rsid);
            auto ca_it = cancer_ss.find(rsid);
            if (rf_it == rf_ss.end() || ca_it == cancer_ss.end()) continue;
            auto pq_snp = pq_it->second.find(rsid);
            double alpha_val = std::numeric_limits<double>::quiet_NaN();
            double se_alpha_val = std::numeric_limits<double>::quiet_NaN();
            bool alpha_observed = false;
            double alpha_rel = 0.0;
            if (pq_snp != pq_it->second.end()) {
                alpha_val = pq_snp->second.beta;
                se_alpha_val = pq_snp->second.se;
                alpha_observed = true;
                alpha_rel = 1.0;
            } else {
                auto ref_it = plink.rsid_to_idx.find(rsid);
                if (ref_it != plink.rsid_to_idx.end()) {
                    ProxyMatch proxy = best_cis_proxy_for_rf(plink, ld_r_cache, cis_all,
                                                             {rsid, rf_it->second.chr, rf_it->second.bp, ref_it->second},
                                                             opts.proxy_r2_thresh, overlap_window_bp);
                    if (proxy.found) {
                        auto proxy_it = pq_it->second.find(proxy.rsid);
                        if (proxy_it != pq_it->second.end()) {
                            alpha_val = proxy.r * proxy_it->second.beta;
                            se_alpha_val = proxy_it->second.se;
                            alpha_observed = true;
                            alpha_rel = std::fabs(proxy.r);
                            prot.nA_proxy++;
                        }
                    }
                }
            }
            prot.setA_rsid.push_back(rsid);
            prot.setA_gamma.push_back(rf_it->second.beta);
            prot.setA_se_gamma.push_back(rf_it->second.se);
            prot.setA_alpha.push_back(alpha_val);
            prot.setA_se_alpha.push_back(se_alpha_val);
            prot.setA_alpha_observed.push_back(alpha_observed);
            prot.setA_alpha_reliability.push_back(alpha_rel);
            prot.setA_Gamma.push_back(ca_it->second.beta);
            prot.setA_se_Gamma.push_back(ca_it->second.se);
        }

        for (const auto& rsid : cis_instruments) {
            if (used_cis.count(rsid)) continue;
            auto pq_snp = pq_it->second.find(rsid);
            auto ca_it = cancer_ss.find(rsid);
            if (pq_snp == pq_it->second.end() || ca_it == cancer_ss.end()) continue;
            prot.setB_rsid.push_back(rsid);
            prot.setB_alpha_cis.push_back(pq_snp->second.beta);
            prot.setB_se_alpha_cis.push_back(pq_snp->second.se);
            prot.setB_Gamma_cis.push_back(ca_it->second.beta);
            prot.setB_se_Gamma_cis.push_back(ca_it->second.se);
        }

        apply_set_qc(prot, plink, plink.n_samples, qc, opts.n_pqtl, opts.n_cancer,
                     total_heidi_removed, opts);
        assign_ld_weights(prot, plink, ld_r2_cache, opts.ld_block_max_size);
        prot.regional_data_complete = prot.regional_cis_rsid.size() >= 2;

        if (prot.nTotal() > 0) n_with_instruments++;
        processed++;
        if (opts.verbose && (processed <= 15 || processed % 250 == 0)) {
            std::cout << "  Protein " << prot.protein_id << ": nA=" << prot.nA()
                      << ", nB=" << prot.nB() << ", nC=" << prot.nC() << "\n";
        }
    }

    std::cout << "  " << n_with_instruments << " / " << proteins.size()
              << " proteins have instruments\n";
    if (qc.heidi_flag) {
        std::cout << "  " << total_heidi_removed
                  << " total instruments removed by HEIDI-outlier across all proteins\n";
    }
}

} // namespace

void run_legacy_pipeline(const Options& opts) {
    std::cout << "Running legacy mode: pre-clumped stacked pQTL input with LD-clumped RF instruments\n";
    if (opts.bfile_prefix.empty()) {
        std::cout << "  No --bfile provided; falling back to legacy behavior without RF clumping\n";
        std::cout << "Reading summary statistics...\n";
        std::map<std::string, SumStat> rf_ss;
        std::map<std::string, double> rf_pval;
        read_sumstats_with_pval(opts.rf_sumstat_file, rf_ss, rf_pval);
        PqtlData pqtl;
        read_pqtl_sumstats(opts.pqtl_sumstat_file, pqtl.ss_by_protein, pqtl.pval_by_protein);
        std::vector<ProteinData> proteins;
        read_protein_info(opts.protein_info_file, proteins);
        auto needed_outcome_rsids = collect_needed_outcome_rsids(rf_pval, pqtl, opts.p_thresh_rf);
        std::map<std::string, SumStat> cancer_ss;
        std::map<std::string, double> cancer_pval;
        read_needed_outcome_sumstats(opts.cancer_sumstat_file, needed_outcome_rsids,
                                     cancer_ss, cancer_pval, opts);
        std::cout << "\nBuilding instrument sets...\n";
        int n_with = 0;
        for (auto& prot : proteins) {
            auto pq_it = pqtl.ss_by_protein.find(prot.protein_id);
            auto pp_it = pqtl.pval_by_protein.find(prot.protein_id);
            if (pq_it != pqtl.ss_by_protein.end() && pp_it != pqtl.pval_by_protein.end()) {
                build_instrument_sets(rf_ss, cancer_ss, pq_it->second, pp_it->second, prot, opts);
            }
            if (prot.nTotal() > 0) n_with++;
        }
        std::cout << "  " << n_with << " / " << proteins.size() << " proteins have instruments.\n";
        copy_mediation_outputs(opts, proteins);
        write_instrument_log(proteins, opts);
        return;
    }

    std::cout << "Loading LD reference panel...\n";
    PlinkData plink;
    if (!plink.load(opts.bfile_prefix)) {
        std::cerr << "Error: failed to load PLINK files\n";
        std::exit(1);
    }

    std::cout << "\nReading summary statistics...\n";
    std::map<std::string, SumStat> rf_ss;
    std::map<std::string, double> rf_pval;
    read_sumstats_with_pval(opts.rf_sumstat_file, rf_ss, rf_pval);

    PqtlData pqtl;
    read_pqtl_sumstats(opts.pqtl_sumstat_file, pqtl.ss_by_protein, pqtl.pval_by_protein);

    std::vector<ProteinData> proteins;
    read_protein_info(opts.protein_info_file, proteins);

    auto needed_outcome_rsids = collect_needed_outcome_rsids(rf_pval, pqtl, opts.p_thresh_rf);
    std::map<std::string, SumStat> cancer_ss;
    std::map<std::string, double> cancer_pval;
    read_needed_outcome_sumstats(opts.cancer_sumstat_file, needed_outcome_rsids,
                                 cancer_ss, cancer_pval, opts);

    std::cout << "\nHarmonizing alleles against reference panel...\n";
    harmonize_map_against_reference(rf_ss, plink, "RF GWAS");
    harmonize_map_against_reference(cancer_ss, plink, "Cancer GWAS");
    harmonize_pqtl_data(pqtl, plink, opts);

    auto qc = build_qc_params(opts);
    auto rf_instruments = select_rf_instruments(rf_ss, rf_pval, cancer_ss, cancer_pval,
                                                plink, qc, opts);
    build_legacy_sets(proteins, rf_ss, cancer_ss, pqtl, rf_instruments, plink, opts);
    copy_mediation_outputs(opts, proteins);
    write_instrument_log(proteins, opts);
}

void run_full_pipeline(const Options& opts) {
    std::cout << "Running full mode: per-protein GWAS manifest with LD-clumped RF and cis instruments\n";
    if (opts.bfile_prefix.empty()) {
        std::cerr << "Error: --bfile is required in full mode\n";
        std::exit(1);
    }

    std::cout << "Loading LD reference panel...\n";
    PlinkData plink;
    if (!plink.load(opts.bfile_prefix)) {
        std::cerr << "Error: failed to load PLINK files\n";
        std::exit(1);
    }

    std::cout << "\nReading summary statistics...\n";
    std::map<std::string, SumStat> rf_ss;
    std::map<std::string, double> rf_pval;
    read_sumstats_with_pval(opts.rf_sumstat_file, rf_ss, rf_pval);

    PqtlData pqtl = read_protein_gwas_manifest(opts.protein_gwas_list_file, opts);

    std::vector<ProteinData> proteins;
    read_protein_info(opts.protein_info_file, proteins);

    auto needed_outcome_rsids = collect_needed_outcome_rsids(rf_pval, pqtl, opts.p_thresh_rf);
    std::map<std::string, SumStat> cancer_ss;
    std::map<std::string, double> cancer_pval;
    read_needed_outcome_sumstats(opts.cancer_sumstat_file, needed_outcome_rsids,
                                 cancer_ss, cancer_pval, opts);

    std::cout << "\nHarmonizing alleles against reference panel...\n";
    harmonize_map_against_reference(rf_ss, plink, "RF GWAS");
    harmonize_map_against_reference(cancer_ss, plink, "Cancer GWAS");
    harmonize_pqtl_data(pqtl, plink, opts);

    auto qc = build_qc_params(opts);
    auto rf_instruments = select_rf_instruments(rf_ss, rf_pval, cancer_ss, cancer_pval,
                                                plink, qc, opts);
    build_full_sets(proteins, rf_ss, cancer_ss, pqtl, rf_instruments, plink, qc, opts);
    copy_mediation_outputs(opts, proteins);
    write_instrument_log(proteins, opts);
}

} // namespace bmediator
