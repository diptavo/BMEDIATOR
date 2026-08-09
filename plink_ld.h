#ifndef PLINK_LD_H
#define PLINK_LD_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cstdint>

namespace bmediator {

// ============================================================================
// PLINK BIM entry
// ============================================================================
struct BimEntry {
    int chr;
    std::string rsid;
    double cm;
    int bp;
    std::string a1;  // allele 1 (usually minor/effect allele in PLINK)
    std::string a2;  // allele 2
    int index;       // position in .bed file (0-based)
};

// ============================================================================
// PLINK FAM entry
// ============================================================================
struct FamEntry {
    std::string fid;
    std::string iid;
    std::string father;
    std::string mother;
    int sex;
    double pheno;
};

// ============================================================================
// PLINK binary dataset handle
// ============================================================================
class PlinkData {
public:
    std::string prefix;
    std::vector<BimEntry> bim;
    std::vector<FamEntry> fam;
    std::map<std::string, int> rsid_to_idx;  // rsid -> bim index
    int n_samples;
    int n_snps;

    bool load(const std::string& prefix);
    bool read_bim(const std::string& fname);
    bool read_fam(const std::string& fname);

    // Extract genotypes for a set of SNPs (by bim index)
    // Returns n_samples x n_snps matrix (0, 1, 2, or -9 for missing)
    // Genotype is coded as count of A1 allele
    std::vector<std::vector<int8_t>> extract_genotypes(
        const std::vector<int>& snp_indices) const;

    // Extract genotypes as doubles (with mean imputation for missing)
    std::vector<std::vector<double>> extract_genotypes_double(
        const std::vector<int>& snp_indices) const;

    // Compute LD r (correlation) between two SNPs
    double compute_ld_r(int snp_idx1, int snp_idx2) const;

    // Compute LD r^2 between two SNPs
    double compute_ld_r2(int snp_idx1, int snp_idx2) const;

    // Compute pairwise LD matrix for a set of SNPs
    // Returns n x n matrix of r values
    std::vector<std::vector<double>> compute_ld_matrix(
        const std::vector<int>& snp_indices) const;

    // Read genotypes for a single SNP from .bed file
    std::vector<int8_t> read_snp_genotypes(int snp_idx) const;

private:
    std::string bed_file;
    int bytes_per_snp;
};

// ============================================================================
// LD clumping (PLINK/GSMR style)
//
// Greedy clumping: sort SNPs by chi-squared (or p-value), then iteratively
// select the top SNP and remove all SNPs in LD (r^2 > threshold) with it
// within a specified window.
// ============================================================================
struct ClumpResult {
    std::vector<int> index_snps;      // indices (in the input list) of retained SNPs
    std::vector<int> removed_snps;    // indices of removed SNPs
};

ClumpResult ld_clump(
    const PlinkData& plink,
    const std::vector<int>& candidate_bim_indices,  // bim indices of candidates
    const std::vector<double>& pvalues,              // p-values (same order as candidates)
    double r2_thresh,                                 // LD r^2 threshold (e.g., 0.1)
    int window_kb                                     // clumping window in kb (e.g., 10000)
);

// ============================================================================
// Allele harmonization
//
// Given effect alleles from summary statistics and from the reference panel,
// check strand consistency and flip effects if needed.
// Returns: +1 if alleles match, -1 if flipped, 0 if incompatible
// ============================================================================
struct HarmonizeResult {
    int action;        // +1 = match, -1 = flip, 0 = exclude
    double beta_adj;   // adjusted beta (flipped if needed)
};

HarmonizeResult harmonize_alleles(
    const std::string& ss_a1, const std::string& ss_a2,  // summary stat alleles
    const std::string& ref_a1, const std::string& ref_a2, // reference panel alleles
    double beta                                            // effect size
);

// Complement base
inline char complement_base(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

inline std::string complement_allele(const std::string& a) {
    std::string result = a;
    for (auto& c : result) c = complement_base(c);
    return result;
}

// ============================================================================
// Frequency check between summary stats and reference panel
// ============================================================================
bool freq_check(double ss_freq, double ref_freq, double max_diff = 0.2);

} // namespace bmediator

#endif // PLINK_LD_H
