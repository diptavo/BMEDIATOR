#include "plink_ld.h"
#include <cassert>

namespace bmediator {

// ============================================================================
// Read .bim file
// Format: CHR  RSID  CM  BP  A1  A2
// ============================================================================
bool PlinkData::read_bim(const std::string& fname) {
    std::ifstream fin(fname);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open BIM file " << fname << "\n";
        return false;
    }
    std::string line;
    int idx = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        BimEntry b;
        if (!(iss >> b.chr >> b.rsid >> b.cm >> b.bp >> b.a1 >> b.a2))
            continue;
        b.index = idx;
        // Uppercase alleles
        for (auto& c : b.a1) c = toupper(c);
        for (auto& c : b.a2) c = toupper(c);
        rsid_to_idx[b.rsid] = idx;
        bim.push_back(b);
        idx++;
    }
    fin.close();
    n_snps = (int)bim.size();
    return true;
}

// ============================================================================
// Read .fam file
// Format: FID  IID  FATHER  MOTHER  SEX  PHENO
// ============================================================================
bool PlinkData::read_fam(const std::string& fname) {
    std::ifstream fin(fname);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open FAM file " << fname << "\n";
        return false;
    }
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        FamEntry f;
        if (!(iss >> f.fid >> f.iid >> f.father >> f.mother >> f.sex >> f.pheno))
            continue;
        fam.push_back(f);
    }
    fin.close();
    n_samples = (int)fam.size();
    return true;
}

// ============================================================================
// Load PLINK dataset (prefix.bed, prefix.bim, prefix.fam)
// ============================================================================
bool PlinkData::load(const std::string& pfx) {
    prefix = pfx;
    bed_file = pfx + ".bed";

    if (!read_bim(pfx + ".bim")) return false;
    if (!read_fam(pfx + ".fam")) return false;

    // Validate .bed file exists and has correct magic number
    std::ifstream bed(bed_file, std::ios::binary);
    if (!bed.is_open()) {
        std::cerr << "Error: cannot open BED file " << bed_file << "\n";
        return false;
    }
    // PLINK BED magic number: 0x6C 0x1B 0x01 (SNP-major mode)
    unsigned char magic[3];
    bed.read(reinterpret_cast<char*>(magic), 3);
    if (magic[0] != 0x6C || magic[1] != 0x1B) {
        std::cerr << "Error: " << bed_file << " is not a valid PLINK BED file\n";
        return false;
    }
    if (magic[2] != 0x01) {
        std::cerr << "Error: " << bed_file << " is in individual-major mode. "
                  << "Please convert to SNP-major mode.\n";
        return false;
    }
    bed.close();

    // Calculate bytes per SNP (each sample takes 2 bits, packed into bytes)
    bytes_per_snp = (n_samples + 3) / 4;

    std::cout << "  Loaded PLINK dataset: " << n_snps << " SNPs, "
              << n_samples << " samples from " << pfx << "\n";
    return true;
}

// ============================================================================
// Read genotypes for a single SNP from .bed file
// Returns vector of int8_t: 0, 1, 2 = dosage of A1; -9 = missing
//
// PLINK BED encoding (2 bits per sample):
//   00 = homozygous A1 (A1/A1) -> dosage 2
//   01 = missing                -> -9
//   10 = heterozygous (A1/A2)  -> dosage 1
//   11 = homozygous A2 (A2/A2) -> dosage 0
// ============================================================================
std::vector<int8_t> PlinkData::read_snp_genotypes(int snp_idx) const {
    std::vector<int8_t> geno(n_samples);
    std::ifstream bed(bed_file, std::ios::binary);

    // Seek to the SNP's data (3-byte header + snp_idx * bytes_per_snp)
    long offset = 3L + (long)snp_idx * bytes_per_snp;
    bed.seekg(offset, std::ios::beg);

    std::vector<unsigned char> buf(bytes_per_snp);
    bed.read(reinterpret_cast<char*>(buf.data()), bytes_per_snp);
    bed.close();

    // Decode
    // Lookup table for 2-bit genotype codes
    static const int8_t geno_lookup[4] = {2, -9, 1, 0}; // 00, 01, 10, 11

    int sample_i = 0;
    for (int byte_i = 0; byte_i < bytes_per_snp && sample_i < n_samples; byte_i++) {
        unsigned char byte = buf[byte_i];
        for (int bit_i = 0; bit_i < 4 && sample_i < n_samples; bit_i++) {
            int code = (byte >> (bit_i * 2)) & 0x03;
            geno[sample_i] = geno_lookup[code];
            sample_i++;
        }
    }

    return geno;
}

// ============================================================================
// Extract genotypes for multiple SNPs
// ============================================================================
std::vector<std::vector<int8_t>> PlinkData::extract_genotypes(
    const std::vector<int>& snp_indices) const {
    std::vector<std::vector<int8_t>> result(snp_indices.size());
    for (size_t i = 0; i < snp_indices.size(); i++) {
        result[i] = read_snp_genotypes(snp_indices[i]);
    }
    return result;
}

// ============================================================================
// Extract genotypes as doubles with mean imputation
// ============================================================================
std::vector<std::vector<double>> PlinkData::extract_genotypes_double(
    const std::vector<int>& snp_indices) const {
    std::vector<std::vector<double>> result(snp_indices.size());
    for (size_t i = 0; i < snp_indices.size(); i++) {
        auto geno = read_snp_genotypes(snp_indices[i]);
        result[i].resize(n_samples);
        // Compute mean for imputation
        double sum = 0.0;
        int count = 0;
        for (int j = 0; j < n_samples; j++) {
            if (geno[j] >= 0) { sum += geno[j]; count++; }
        }
        double mean = (count > 0) ? sum / count : 0.0;
        for (int j = 0; j < n_samples; j++) {
            result[i][j] = (geno[j] >= 0) ? (double)geno[j] : mean;
        }
    }
    return result;
}

// ============================================================================
// Compute LD r between two SNPs (Pearson correlation of genotypes)
// ============================================================================
double PlinkData::compute_ld_r(int snp_idx1, int snp_idx2) const {
    auto g1 = read_snp_genotypes(snp_idx1);
    auto g2 = read_snp_genotypes(snp_idx2);

    double sum1 = 0, sum2 = 0, sum12 = 0, sum11 = 0, sum22 = 0;
    int n = 0;
    for (int i = 0; i < n_samples; i++) {
        if (g1[i] < 0 || g2[i] < 0) continue;
        double x = g1[i], y = g2[i];
        sum1 += x; sum2 += y;
        sum12 += x * y;
        sum11 += x * x; sum22 += y * y;
        n++;
    }
    if (n < 10) return 0.0;

    double mean1 = sum1 / n, mean2 = sum2 / n;
    double var1 = sum11 / n - mean1 * mean1;
    double var2 = sum22 / n - mean2 * mean2;
    if (var1 < 1e-10 || var2 < 1e-10) return 0.0;

    double cov = sum12 / n - mean1 * mean2;
    return cov / std::sqrt(var1 * var2);
}

double PlinkData::compute_ld_r2(int snp_idx1, int snp_idx2) const {
    double r = compute_ld_r(snp_idx1, snp_idx2);
    return r * r;
}

// ============================================================================
// Compute LD matrix for a set of SNPs
// ============================================================================
std::vector<std::vector<double>> PlinkData::compute_ld_matrix(
    const std::vector<int>& snp_indices) const {
    int m = (int)snp_indices.size();

    // Read all genotypes first (much faster than reading each pair)
    auto geno = extract_genotypes_double(snp_indices);

    // Standardize
    std::vector<double> means(m), sds(m);
    for (int i = 0; i < m; i++) {
        double s = 0, s2 = 0;
        for (int j = 0; j < n_samples; j++) {
            s += geno[i][j]; s2 += geno[i][j] * geno[i][j];
        }
        means[i] = s / n_samples;
        double var = s2 / n_samples - means[i] * means[i];
        sds[i] = (var > 1e-10) ? std::sqrt(var) : 0.0;
    }

    // Compute correlation matrix
    std::vector<std::vector<double>> ld(m, std::vector<double>(m, 0.0));
    for (int i = 0; i < m; i++) {
        ld[i][i] = 1.0;
        if (sds[i] < 1e-10) continue;
        for (int j = i + 1; j < m; j++) {
            if (sds[j] < 1e-10) continue;
            double cov = 0;
            for (int k = 0; k < n_samples; k++) {
                cov += (geno[i][k] - means[i]) * (geno[j][k] - means[j]);
            }
            cov /= n_samples;
            double r = cov / (sds[i] * sds[j]);
            ld[i][j] = r;
            ld[j][i] = r;
        }
    }
    return ld;
}

// ============================================================================
// LD clumping (PLINK/GSMR greedy algorithm)
//
// 1. Sort candidate SNPs by p-value (ascending = most significant first)
// 2. Select top SNP as index SNP
// 3. Remove all candidates within window_kb that have r^2 > r2_thresh
// 4. Move to next unremoved SNP, repeat
// ============================================================================
ClumpResult ld_clump(
    const PlinkData& plink,
    const std::vector<int>& candidate_bim_indices,
    const std::vector<double>& pvalues,
    double r2_thresh,
    int window_kb) {

    int n = (int)candidate_bim_indices.size();
    assert(n == (int)pvalues.size());

    // Sort by p-value (index in the candidates list)
    std::vector<int> order(n);
    for (int i = 0; i < n; i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return pvalues[a] < pvalues[b];
    });

    std::vector<bool> removed(n, false);
    ClumpResult result;

    int window_bp = window_kb * 1000;

    for (int i = 0; i < n; i++) {
        int idx = order[i];
        if (removed[idx]) continue;

        // This SNP becomes an index SNP
        result.index_snps.push_back(idx);
        int bim_i = candidate_bim_indices[idx];

        // Remove all SNPs in LD with this one
        for (int j = i + 1; j < n; j++) {
            int jdx = order[j];
            if (removed[jdx]) continue;

            int bim_j = candidate_bim_indices[jdx];

            // Check window
            if (plink.bim[bim_i].chr != plink.bim[bim_j].chr) continue;
            if (std::abs(plink.bim[bim_i].bp - plink.bim[bim_j].bp) > window_bp) continue;

            // Check LD
            double r2 = plink.compute_ld_r2(bim_i, bim_j);
            if (r2 > r2_thresh) {
                removed[jdx] = true;
                result.removed_snps.push_back(jdx);
            }
        }
    }

    return result;
}

// ============================================================================
// Allele harmonization
//
// Handles: direct match, flip (A1<->A2), strand flip, strand+allele flip
// Returns action: +1 (match), -1 (flip beta sign), 0 (incompatible)
// ============================================================================
HarmonizeResult harmonize_alleles(
    const std::string& ss_a1, const std::string& ss_a2,
    const std::string& ref_a1, const std::string& ref_a2,
    double beta) {

    HarmonizeResult res;
    res.beta_adj = beta;
    res.action = 0;

    // Direct match: ss_a1=ref_a1, ss_a2=ref_a2
    if (ss_a1 == ref_a1 && ss_a2 == ref_a2) {
        res.action = 1;
        return res;
    }

    // Flipped alleles: ss_a1=ref_a2, ss_a2=ref_a1
    if (ss_a1 == ref_a2 && ss_a2 == ref_a1) {
        res.action = -1;
        res.beta_adj = -beta;
        return res;
    }

    // Strand flip
    std::string ss_a1c = complement_allele(ss_a1);
    std::string ss_a2c = complement_allele(ss_a2);

    // Strand flip + direct match
    if (ss_a1c == ref_a1 && ss_a2c == ref_a2) {
        res.action = 1;
        return res;
    }

    // Strand flip + allele flip
    if (ss_a1c == ref_a2 && ss_a2c == ref_a1) {
        res.action = -1;
        res.beta_adj = -beta;
        return res;
    }

    // Incompatible (ambiguous A/T or C/G SNPs could be problematic)
    // For A/T and C/G SNPs, we cannot resolve strand ambiguity
    res.action = 0;
    return res;
}

// ============================================================================
// Frequency check
// ============================================================================
bool freq_check(double ss_freq, double ref_freq, double max_diff) {
    // Check if frequencies are concordant (within max_diff)
    // Also check flipped: if |ss_freq - (1-ref_freq)| < max_diff, alleles may be flipped
    double diff = std::fabs(ss_freq - ref_freq);
    return diff <= max_diff;
}

} // namespace bmediator
