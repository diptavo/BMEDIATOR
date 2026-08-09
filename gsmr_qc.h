#ifndef GSMR_QC_H
#define GSMR_QC_H

#include "plink_ld.h"
#include "bmediator.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>

namespace bmediator {

// ============================================================================
// QC configuration (mirrors GSMR control parameters)
// ============================================================================
struct QCParams {
    // MAF filtering
    double maf_min;               // minimum MAF in reference panel (default 0.01)

    // Allele frequency concordance
    double freq_diff_thresh;      // max |freq_ss - freq_ref| (default 0.2)

    // Palindromic SNP handling
    bool remove_palindromic;      // remove A/T and C/G SNPs (default true)
    double palindromic_freq_diff; // for palindromic SNPs, require freq_diff < this (default 0.15)

    // SE / effect size sanity
    double se_min;                // minimum SE (remove SNPs with SE < this, default 1e-6)
    double se_max;                // maximum SE (remove outlier SEs, default 10.0)
    double effect_max;            // maximum |beta| (remove outlier effects, default 10.0)

    // Sample size filtering
    double n_min_frac;            // remove SNPs with N < n_min_frac * median(N) (default 0.5)

    // LD-based filtering
    double ld_r2_thresh;          // LD r2 threshold for clumping (default 0.05, GSMR uses 0.05)
    double ld_fdr_thresh;         // FDR threshold for removing chance LD (default 0.05)
    int    clump_window_kb;       // clumping window in kb (default 10000)

    // HEIDI outlier filtering
    bool   heidi_flag;            // run HEIDI-outlier test (default true)
    double heidi_thresh;          // p-value threshold for HEIDI-outlier (default 0.01)
    int    heidi_min_snps;        // minimum instruments for HEIDI (default 3)
    bool   heidi_global;          // use global multi-SNP HEIDI (GSMR2, default true)

    // Steiger directionality
    bool   steiger_flag;          // run Steiger filtering (default true)

    // Minimum instruments
    int    min_instruments;       // skip protein/trait if fewer instruments (default 5)

    QCParams() :
        maf_min(0.01),
        freq_diff_thresh(0.2),
        remove_palindromic(true),
        palindromic_freq_diff(0.15),
        se_min(1e-6), se_max(10.0), effect_max(10.0),
        n_min_frac(0.5),
        ld_r2_thresh(0.05),
        ld_fdr_thresh(0.05),
        clump_window_kb(10000),
        heidi_flag(true),
        heidi_thresh(0.01),
        heidi_min_snps(3),
        heidi_global(true),
        steiger_flag(true),
        min_instruments(5) {}
};

// ============================================================================
// QC report for a single trait or protein
// ============================================================================
struct QCReport {
    int n_input;
    int n_after_maf;
    int n_after_freq_check;
    int n_after_palindromic;
    int n_after_se_check;
    int n_after_n_check;
    int n_after_clump;
    int n_after_heidi;
    int n_after_steiger;
    int n_final;
    std::vector<std::string> removed_palindromic;
    std::vector<std::string> removed_freq;
    std::vector<std::string> removed_se;
    std::vector<std::string> removed_heidi;
    std::vector<std::string> removed_steiger;
};

// ============================================================================
// Check if SNP is palindromic (A/T or C/G)
// ============================================================================
inline bool is_palindromic(const std::string& a1, const std::string& a2) {
    if (a1.size() != 1 || a2.size() != 1) return false;
    char c1 = toupper(a1[0]), c2 = toupper(a2[0]);
    return (c1 == 'A' && c2 == 'T') || (c1 == 'T' && c2 == 'A') ||
           (c1 == 'C' && c2 == 'G') || (c1 == 'G' && c2 == 'C');
}

// ============================================================================
// Compute MAF from reference panel genotypes
// ============================================================================
double compute_ref_maf(const PlinkData& plink, int bim_idx);

// ============================================================================
// Compute allele frequency from reference panel
// Returns frequency of A1 (the allele coded in the BED file)
// ============================================================================
double compute_ref_freq(const PlinkData& plink, int bim_idx);

// ============================================================================
// Pre-MR QC pipeline for a set of SNPs
//
// Takes a set of candidate SNP rsids (present in both summary stats and ref),
// applies all QC filters, returns the surviving SNP set.
//
// Steps:
//   1. MAF filter (reference panel)
//   2. Allele frequency concordance (sumstat vs reference)
//   3. Palindromic SNP removal (or strict freq concordance)
//   4. SE and effect size sanity checks
//   5. Sample size filtering
// ============================================================================
std::vector<std::string> run_snp_qc(
    const std::map<std::string, SumStat>& ss,
    const std::map<std::string, double>& pvals,
    const PlinkData& plink,
    const QCParams& qc,
    QCReport& report);

// ============================================================================
// HEIDI-outlier test (GSMR style)
//
// Tests whether the ratio bzy/bzx is homogeneous across instruments.
// Under the causal model H0: bzy_j = bxy * bzx_j for all j.
// The HEIDI test statistic for SNP j vs index SNP i:
//
//   d_j = bzy_j - bxy_hat * bzx_j
//   var(d_j) = var(bzy_j) + bxy_hat^2 * var(bzx_j)
//              - 2*bxy_hat*rho_j*se_bzy_j*se_bzx_j  [if overlap exists]
//   T_j = d_j^2 / var(d_j) ~ chi^2(1)
//
// For the global HEIDI (GSMR2), we compute a joint test statistic.
//
// Inputs:
//   bzx, se_bzx: SNP effects on exposure (vector)
//   bzy, se_bzy: SNP effects on outcome (vector)
//   ld_matrix: correlation matrix (r) of instruments from reference panel
//   n_ref: reference panel sample size
//   bxy_hat: estimated causal effect (from IVW or similar)
//
// Returns: vector of booleans (true = keep, false = HEIDI outlier)
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
    bool global_heidi);

// ============================================================================
// Global HEIDI test (GSMR2 multi-SNP version)
//
// Joint test of homogeneity across all instruments simultaneously.
// Uses the full LD matrix to account for correlation between test statistics.
//
// Returns: p-value for the global HEIDI test
// ============================================================================
double global_heidi_pvalue(
    const std::vector<double>& bzx,
    const std::vector<double>& se_bzx,
    const std::vector<double>& bzy,
    const std::vector<double>& se_bzy,
    const std::vector<std::vector<double>>& ld_matrix,
    int n_ref,
    double bxy_hat);

// ============================================================================
// Steiger directionality test
//
// Tests whether the variance explained in the exposure is greater than
// the variance explained in the outcome, which should hold if the
// causal direction is exposure -> outcome (not reverse).
//
// For each instrument j:
//   r2_xj = bzx_j^2 / (bzx_j^2 + N_x * se_bzx_j^2)
//   r2_yj = bzy_j^2 / (bzy_j^2 + N_y * se_bzy_j^2)
//
// If median(r2_xj) < median(r2_yj), the direction may be wrong.
//
// Returns: true if direction is correct, false if reverse
// ============================================================================
bool steiger_direction_test(
    const std::vector<double>& bzx,
    const std::vector<double>& se_bzx,
    double n_x,
    const std::vector<double>& bzy,
    const std::vector<double>& se_bzy,
    double n_y);

// ============================================================================
// IVW estimate of causal effect (used as input to HEIDI)
//
// bxy_ivw = sum(bzx_j * bzy_j / se_bzy_j^2) / sum(bzx_j^2 / se_bzy_j^2)
// ============================================================================
double ivw_estimate(
    const std::vector<double>& bzx,
    const std::vector<double>& bzy,
    const std::vector<double>& se_bzy);

// ============================================================================
// Chance LD removal (GSMR ld_fdr_thresh)
//
// For each pair of instruments, test whether their LD in the reference
// panel is significantly different from zero using Fisher's z-transform.
// Remove one of a pair if significant LD remains after clumping
// (accounts for finite reference panel size).
// ============================================================================
std::vector<bool> filter_chance_ld(
    const std::vector<std::vector<double>>& ld_matrix,
    int n_ref,
    double fdr_thresh);

} // namespace bmediator

#endif // GSMR_QC_H
