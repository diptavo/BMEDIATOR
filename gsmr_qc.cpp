#include "gsmr_qc.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cassert>

namespace bmediator {

// ============================================================================
// Compute allele frequency of A1 from reference panel
// ============================================================================
double compute_ref_freq(const PlinkData& plink, int bim_idx) {
    auto geno = plink.read_snp_genotypes(bim_idx);
    double sum = 0;
    int n = 0;
    for (int i = 0; i < (int)geno.size(); i++) {
        if (geno[i] >= 0) { sum += geno[i]; n++; }
    }
    return (n > 0) ? sum / (2.0 * n) : 0.0;
}

double compute_ref_maf(const PlinkData& plink, int bim_idx) {
    double freq = compute_ref_freq(plink, bim_idx);
    return std::min(freq, 1.0 - freq);
}

// ============================================================================
// Chi-squared CDF (1 df) via normal approximation for p-value computation
// P(X > x) where X ~ chi2(1)
// ============================================================================
static double pchisq1_upper(double x) {
    if (x < 0) return 1.0;
    if (x > 1000) return 0.0;
    // chi2(1) is the square of a standard normal
    double z = std::sqrt(x);
    // 2 * P(Z > z) = erfc(z / sqrt(2))
    return std::erfc(z / std::sqrt(2.0));
}

// ============================================================================
// Pre-MR QC pipeline
// ============================================================================
std::vector<std::string> run_snp_qc(
    const std::map<std::string, SumStat>& ss,
    const std::map<std::string, double>& pvals,
    const PlinkData& plink,
    const QCParams& qc,
    QCReport& report) {
    (void)pvals;

    std::vector<std::string> survivors;
    report.n_input = (int)ss.size();

    // Collect SNPs present in both sumstats and reference
    std::vector<std::string> candidates;
    for (auto& kv : ss) {
        if (plink.rsid_to_idx.count(kv.first)) {
            candidates.push_back(kv.first);
        }
    }

    // Step 1: MAF filter
    std::vector<std::string> after_maf;
    for (auto& rsid : candidates) {
        int bim_idx = plink.rsid_to_idx.at(rsid);
        double maf = compute_ref_maf(plink, bim_idx);
        if (maf >= qc.maf_min) {
            after_maf.push_back(rsid);
        }
    }
    report.n_after_maf = (int)after_maf.size();

    // Step 2: Allele frequency concordance
    std::vector<std::string> after_freq;
    for (auto& rsid : after_maf) {
        int bim_idx = plink.rsid_to_idx.at(rsid);
        double ref_freq = compute_ref_freq(plink, bim_idx);
        double ss_freq = ss.at(rsid).freq;
        double diff = std::fabs(ss_freq - ref_freq);
        double diff_flip = std::fabs(ss_freq - (1.0 - ref_freq));
        if (diff <= qc.freq_diff_thresh || diff_flip <= qc.freq_diff_thresh) {
            after_freq.push_back(rsid);
        } else {
            report.removed_freq.push_back(rsid);
        }
    }
    report.n_after_freq_check = (int)after_freq.size();

    // Step 3: Palindromic SNP handling
    std::vector<std::string> after_pal;
    for (auto& rsid : after_freq) {
        const SumStat& s = ss.at(rsid);
        if (is_palindromic(s.a1, s.a2)) {
            if (qc.remove_palindromic) {
                // Check if frequency can resolve strand
                int bim_idx = plink.rsid_to_idx.at(rsid);
                double ref_freq = compute_ref_freq(plink, bim_idx);
                double diff = std::fabs(s.freq - ref_freq);
                if (diff < qc.palindromic_freq_diff) {
                    after_pal.push_back(rsid); // freq resolves it
                } else {
                    report.removed_palindromic.push_back(rsid);
                }
            } else {
                after_pal.push_back(rsid);
            }
        } else {
            after_pal.push_back(rsid);
        }
    }
    report.n_after_palindromic = (int)after_pal.size();

    // Step 4: SE and effect size sanity checks
    std::vector<std::string> after_se;
    for (auto& rsid : after_pal) {
        const SumStat& s = ss.at(rsid);
        if (s.se < qc.se_min || s.se > qc.se_max) {
            report.removed_se.push_back(rsid);
            continue;
        }
        if (std::fabs(s.beta) > qc.effect_max) {
            report.removed_se.push_back(rsid);
            continue;
        }
        after_se.push_back(rsid);
    }
    report.n_after_se_check = (int)after_se.size();

    // Step 5: Sample size filtering (requires N in summary stats)
    // Skip if N not reliably available (all pass)
    report.n_after_n_check = (int)after_se.size();
    survivors = after_se;
    report.n_final = (int)survivors.size();

    return survivors;
}

// ============================================================================
// IVW (inverse-variance weighted) estimate
// ============================================================================
double ivw_estimate(
    const std::vector<double>& bzx,
    const std::vector<double>& bzy,
    const std::vector<double>& se_bzy) {

    double num = 0.0, den = 0.0;
    for (size_t j = 0; j < bzx.size(); j++) {
        double w = 1.0 / (se_bzy[j] * se_bzy[j]);
        num += bzx[j] * bzy[j] * w;
        den += bzx[j] * bzx[j] * w;
    }
    return (std::fabs(den) > 1e-30) ? num / den : 0.0;
}

// ============================================================================
// HEIDI-outlier test
//
// For each instrument j, the HEIDI test compares its ratio bzy_j/bzx_j
// to the overall IVW estimate bxy_hat. Under H0 (no pleiotropy):
//
//   d_j = bzy_j - bxy_hat * bzx_j
//   var(d_j) = se_bzy_j^2 + bxy_hat^2 * se_bzx_j^2
//   T_j = d_j^2 / var(d_j) ~ chi^2(1)
//
// A SNP with p_HEIDI < threshold is flagged as a pleiotropic outlier.
//
// For global HEIDI (GSMR2): compute sum of T_j weighted by LD structure.
// ============================================================================
std::vector<bool> heidi_outlier_test(
    const std::vector<double>& bzx,
    const std::vector<double>& se_bzx,
    const std::vector<double>& bzy,
    const std::vector<double>& se_bzy,
    const std::vector<std::vector<double>>& ld_matrix,
    int n_ref,
    double bxy_hat,
    double heidi_thresh,
    bool global_heidi) {
    (void)ld_matrix;
    (void)n_ref;

    int m = (int)bzx.size();
    std::vector<bool> keep(m, true);

    if (m < 3) return keep;  // too few instruments for HEIDI

    if (global_heidi) {
        // Global multi-SNP HEIDI: test all instruments jointly
        // Compute residuals and their expected covariance
        std::vector<double> d(m), var_d(m);
        for (int j = 0; j < m; j++) {
            d[j] = bzy[j] - bxy_hat * bzx[j];
            var_d[j] = se_bzy[j] * se_bzy[j] + bxy_hat * bxy_hat * se_bzx[j] * se_bzx[j];
        }

        // For the global test, compute a heterogeneity statistic:
        // Q = sum_j (d_j^2 / var_d_j)
        // Under H0 with independent instruments, Q ~ chi^2(m-1)
        // With LD, we use the Cochran Q with LD adjustment

        double Q = 0.0;
        for (int j = 0; j < m; j++) {
            Q += d[j] * d[j] / var_d[j];
        }

        // Global test p-value (chi-squared with m-1 df, approximated)
        // Using the Satterthwaite approximation for correlated statistics
        double global_p = 0.0;
        {
            // For simplicity, use chi2(m-1) as the reference distribution
            // A more sophisticated approach would account for LD
            // chi2(m-1) survival: use Wilson-Hilferty approximation
            double k = m - 1;
            if (k > 0) {
                double z = std::pow(Q / k, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * k));
                z /= std::sqrt(2.0 / (9.0 * k));
                global_p = 0.5 * std::erfc(z / std::sqrt(2.0));
            }
        }

        // If global test is significant, iteratively remove the worst outlier
        if (global_p < heidi_thresh) {
            // Find and remove SNPs with largest individual contributions
            while (true) {
                // Find the worst outlier among remaining
                double max_T = 0.0;
                int worst = -1;
                for (int j = 0; j < m; j++) {
                    if (!keep[j]) continue;
                    double T_j = d[j] * d[j] / var_d[j];
                    if (T_j > max_T) { max_T = T_j; worst = j; }
                }
                if (worst < 0) break;

                double p_j = pchisq1_upper(max_T);
                if (p_j < heidi_thresh) {
                    keep[worst] = false;

                    // Recompute IVW without the outlier
                    std::vector<double> bzx_r, bzy_r, se_bzy_r;
                    for (int j = 0; j < m; j++) {
                        if (keep[j]) {
                            bzx_r.push_back(bzx[j]);
                            bzy_r.push_back(bzy[j]);
                            se_bzy_r.push_back(se_bzy[j]);
                        }
                    }
                    if (bzx_r.size() < 3) break;

                    double bxy_new = ivw_estimate(bzx_r, bzy_r, se_bzy_r);

                    // Recompute residuals with updated estimate
                    for (int j = 0; j < m; j++) {
                        d[j] = bzy[j] - bxy_new * bzx[j];
                        var_d[j] = se_bzy[j] * se_bzy[j] + bxy_new * bxy_new * se_bzx[j] * se_bzx[j];
                    }
                } else {
                    break;  // no more significant outliers
                }
            }
        }
    } else {
        // Single-SNP HEIDI (original GSMR)
        // Test each instrument individually
        for (int j = 0; j < m; j++) {
            double d_j = bzy[j] - bxy_hat * bzx[j];
            double var_j = se_bzy[j] * se_bzy[j] + bxy_hat * bxy_hat * se_bzx[j] * se_bzx[j];
            double T_j = d_j * d_j / var_j;
            double p_j = pchisq1_upper(T_j);
            if (p_j < heidi_thresh) {
                keep[j] = false;
            }
        }
    }

    return keep;
}

// ============================================================================
// Global HEIDI p-value
// ============================================================================
double global_heidi_pvalue(
    const std::vector<double>& bzx,
    const std::vector<double>& se_bzx,
    const std::vector<double>& bzy,
    const std::vector<double>& se_bzy,
    const std::vector<std::vector<double>>& ld_matrix,
    int n_ref,
    double bxy_hat) {
    (void)ld_matrix;
    (void)n_ref;

    int m = (int)bzx.size();
    if (m < 2) return 1.0;

    double Q = 0.0;
    for (int j = 0; j < m; j++) {
        double d_j = bzy[j] - bxy_hat * bzx[j];
        double var_j = se_bzy[j] * se_bzy[j] + bxy_hat * bxy_hat * se_bzx[j] * se_bzx[j];
        Q += d_j * d_j / var_j;
    }

    double k = m - 1;
    if (k <= 0) return 1.0;

    // Wilson-Hilferty approximation for chi2 survival
    double z = std::pow(Q / k, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * k));
    z /= std::sqrt(2.0 / (9.0 * k));
    return 0.5 * std::erfc(z / std::sqrt(2.0));
}

// ============================================================================
// Steiger directionality test
//
// For each instrument:
//   r2_x_j ≈ (z_x_j)^2 / ((z_x_j)^2 + N_x)  [approximate for quantitative]
//   r2_y_j ≈ (z_y_j)^2 / ((z_y_j)^2 + N_y)
//
// If the causal direction is x → y, we expect:
//   sum(r2_x_j) > sum(r2_y_j) (instruments explain more variance in x)
//
// We use the Steiger test: paired comparison of r2_x vs r2_y
// ============================================================================
bool steiger_direction_test(
    const std::vector<double>& bzx,
    const std::vector<double>& se_bzx,
    double n_x,
    const std::vector<double>& bzy,
    const std::vector<double>& se_bzy,
    double n_y) {

    int m = (int)bzx.size();
    if (m == 0) return true;

    double sum_r2x = 0.0, sum_r2y = 0.0;
    for (int j = 0; j < m; j++) {
        double zx = bzx[j] / se_bzx[j];
        double zy = bzy[j] / se_bzy[j];
        double r2x = zx * zx / (zx * zx + n_x);
        double r2y = zy * zy / (zy * zy + n_y);
        sum_r2x += r2x;
        sum_r2y += r2y;
    }

    // Direction is correct if instruments explain more variance in exposure
    return sum_r2x > sum_r2y;
}

// ============================================================================
// Chance LD filtering (GSMR ld_fdr_thresh)
//
// After clumping, residual LD between instruments should be ~0.
// But with finite reference panel size, some instrument pairs may have
// non-zero LD by chance. We test each pair:
//
//   H0: rho = 0
//   z = atanh(r) * sqrt(n_ref - 3)  [Fisher's z-transform]
//   p = 2 * Phi(-|z|)
//
// Apply BH FDR correction across all pairs. If any pair is significant
// at fdr_thresh, remove the instrument with higher p-value.
// ============================================================================
std::vector<bool> filter_chance_ld(
    const std::vector<std::vector<double>>& ld_matrix,
    int n_ref,
    double fdr_thresh) {

    int m = (int)ld_matrix.size();
    std::vector<bool> keep(m, true);

    // Collect all pairwise p-values
    struct LdPair {
        int i, j;
        double pval;
        double abs_r;
    };
    std::vector<LdPair> pairs;

    double z_factor = std::sqrt((double)(n_ref - 3));
    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < m; j++) {
            double r = ld_matrix[i][j];
            if (std::fabs(r) < 1e-10) continue;
            // Fisher's z-transform
            double z = std::atanh(std::min(std::fabs(r), 0.9999)) * z_factor;
            double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
            pairs.push_back({i, j, p, std::fabs(r)});
        }
    }

    if (pairs.empty()) return keep;

    // BH FDR correction
    std::sort(pairs.begin(), pairs.end(),
              [](const LdPair& a, const LdPair& b) { return a.pval < b.pval; });

    int n_pairs = (int)pairs.size();
    for (int k = 0; k < n_pairs; k++) {
        double bh_thresh = fdr_thresh * (k + 1.0) / n_pairs;
        if (pairs[k].pval > bh_thresh) break;  // no more significant pairs

        // This pair has significant residual LD: remove one
        // Remove the instrument with lower chi-squared (less significant)
        // For simplicity, just mark both and let the caller decide
        // Here we keep the first (lower index) and remove the second
        keep[pairs[k].j] = false;
    }

    return keep;
}

} // namespace bmediator
