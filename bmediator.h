#ifndef BMEDIATOR_H
#define BMEDIATOR_H

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <limits>
#include <unordered_set>

// ============================================================================
// BMEDIATOR: Bayesian Mediation MR
//
// A tool for identifying mediating plasma proteins between risk factors
// and disease outcomes using a Bayesian framework with summary statistics.
//
// Model: RF -> PP -> Cancer  (mediation)
//        RF -------> Cancer  (direct effect)
//
// Six protein-pathway scenarios:
//   M=0: Null
//   M=1: Partial mediation (RF->PP->Cancer, residual RF->Cancer allowed)
//   M=2: RF->PP only
//   M=3: RF->Cancer direct/residual only
//   M=4: PP->Cancer only (protein target, not RF mediator)
//   M=5: Correlated/shared pleiotropy masquerading as mediation
//
// Reference: Dey et al. (2026)
// ============================================================================

namespace bmediator {

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
constexpr double LOG2PI = 1.8378770664093453;
constexpr double TINY   = 1e-300;
constexpr double EPS    = 1e-10;
constexpr int    MAX_CAVI_ITER = 200;
constexpr double ELBO_TOL = 1e-6;

// ---------------------------------------------------------------------------
// Input data structures
// ---------------------------------------------------------------------------

// Summary statistics for a single SNP-trait association
struct SumStat {
    std::string rsid;
    std::string a1, a2;    // effect allele, other allele
    double beta;           // effect estimate
    double se;             // standard error
    double freq;           // allele frequency (optional, for diagnostics)
    int    chr;
    int    bp;
};

struct RegionalSignalPairResult {
    int protein_signal = 0;
    int outcome_signal = 0;
    std::string protein_lead;
    std::string outcome_lead;
    double pp_h0 = 0.0;
    double pp_h1 = 0.0;
    double pp_h2 = 0.0;
    double pp_h3 = 0.0;
    double pp_h4 = 0.0;
    double shared_given_both = 0.0;
    double lead_pair_r2 = 0.0;
    double max_credible_set_pair_r2 = 0.0;
    std::string interpretation;
};

// Per-protein data: three instrument sets (A, B, C)
struct ProteinData {
    std::string protein_id;
    std::string gene_name;
    int gene_chr;
    int gene_start;
    int gene_end;

    // Set A: RF-only instruments (genome-wide sig for RF, NOT in cis region)
    // For each SNP k in A:
    //   gamma_k  = SNP effect on RF
    //   alpha_k  = SNP effect on PP (from pQTL study)
    //   Gamma_k  = SNP effect on Cancer
    std::vector<std::string> setA_rsid;
    std::vector<double> setA_gamma;      // RF effect
    std::vector<double> setA_se_gamma;
    std::vector<double> setA_alpha;      // PP effect
    std::vector<double> setA_se_alpha;
    std::vector<bool> setA_alpha_observed;
    std::vector<double> setA_alpha_reliability;
    std::vector<double> setA_Gamma;      // Cancer effect
    std::vector<double> setA_se_Gamma;

    // Set B: Cis-only instruments (cis-pQTL for PP, NOT sig for RF)
    std::vector<std::string> setB_rsid;
    std::vector<double> setB_alpha_cis;  // cis effect on PP
    std::vector<double> setB_se_alpha_cis;
    std::vector<double> setB_Gamma_cis;  // Cancer effect
    std::vector<double> setB_se_Gamma_cis;

    // Set C: Overlapping instruments (both RF-sig AND in cis region)
    std::vector<std::string> setC_rsid;
    std::vector<double> setC_gamma;
    std::vector<double> setC_se_gamma;
    std::vector<double> setC_alpha;      // total PP effect (cis + RF-mediated)
    std::vector<double> setC_se_alpha;
    std::vector<double> setC_alpha_reliability;
    std::vector<double> setC_Gamma;
    std::vector<double> setC_se_Gamma;

    // Unpruned cis-region statistics used for shared-vs-distinct causal evidence.
    std::vector<std::string> regional_cis_rsid;
    std::vector<double> regional_pp_beta;
    std::vector<double> regional_pp_se;
    std::vector<double> regional_outcome_beta;
    std::vector<double> regional_outcome_se;
    std::vector<int> regional_bim_index;
    bool regional_data_complete = false;
    bool ld_reference_used = false;
    bool regional_multisignal_evaluated = false;
    std::string regional_method = "single";
    int regional_protein_signals = 0;
    int regional_outcome_signals = 0;
    int regional_independent_shared_signals = 0;
    double regional_best_pp_shared = 0.0;
    double regional_best_pp_distinct = 0.0;
    double regional_best_shared_given_both = 0.0;
    double regional_best_cs_pair_r2 = 0.0;
    std::string regional_multisignal_interpretation;
    std::vector<RegionalSignalPairResult> regional_signal_pairs;

    int nC_exact = 0;
    int nC_proxy = 0;
    int nA_proxy = 0;

    // LD-aware weights used by the approximate likelihood
    std::vector<double> ld_weight_alpha_ac;    // weights for A/C alpha equations
    std::vector<double> ld_weight_cancer_union; // weights for A/C/B cancer equations

    // Signed LD correlation matrices for factorized structural inference.
    // Empty matrices imply identity LD (legacy input without a reference).
    std::vector<std::vector<double>> setA_ld;
    std::vector<std::vector<double>> setB_ld;
    double setAB_max_r2 = 0.0;

    int nA() const { return (int)setA_rsid.size(); }
    int nB() const { return (int)setB_rsid.size(); }
    int nC() const { return (int)setC_rsid.size(); }
    int nTotal() const { return nA() + nB() + nC(); }
};

// ---------------------------------------------------------------------------
// Variational parameters for a single protein under a given scenario
// ---------------------------------------------------------------------------
struct VarParams {
    // Causal effects
    double mu_beta1, s2_beta1;   // RF -> PP
    double mu_beta2, s2_beta2;   // PP -> Cancer
    double mu_beta3, s2_beta3;   // RF -> Cancer (direct)
    double prior_sigma2_beta1, prior_sigma2_beta2, prior_sigma2_beta3;

    // Per-SNP spike-and-slab: pleiotropy terms
    // Set A + C: delta_k (stage 1 pleiotropy)
    std::vector<double> mu_delta, s2_delta, omega_delta;
    // Set B: phi_l (cis pleiotropy)
    std::vector<double> mu_phi, s2_phi, omega_phi;
    // Set A + C: psi_k (cancer pleiotropy)
    std::vector<double> mu_psi, s2_psi, omega_psi;

    // Scenario M5: bivariate (delta_k, psi_k) posterior
    // Joint inclusion probability (replaces independent omega_delta, omega_psi under M5)
    std::vector<double> omega_joint;     // P(z_k = 1) for the bivariate pair
    std::vector<double> cov_delta_psi;   // posterior covariance between delta_k and psi_k
    double rho_delta_psi;                // correlation parameter (prior or estimated)

    double elbo;
};

// ---------------------------------------------------------------------------
// Hyperparameters (shared across all proteins, updated via EB)
// ---------------------------------------------------------------------------
struct Hyperparams {
    // Prior scenario probabilities
    double p0, p1, p2, p3, p4, p5;

    // Prior variances on causal effects (slab)
    double sigma2_beta1;
    double sigma2_beta2;
    double sigma2_beta3;

    // Spike-and-slab: pleiotropy proportions and variances
    double pi1;          // proportion of pleiotropic RF instruments (delta)
    double tau2_1;       // slab variance for delta
    double pi2_cis;      // proportion of pleiotropic cis instruments (phi)
    double tau2_2_cis;   // slab variance for phi
    double pi3;          // proportion of pleiotropic RF instruments (psi)
    double tau2_3;       // slab variance for psi

    // Scenario M5: correlated pleiotropy prior
    double rho_prior;    // prior on delta-psi correlation

    Hyperparams() :
        p0(0.85), p1(0.03), p2(0.05), p3(0.03), p4(0.02), p5(0.02),
        sigma2_beta1(0.1), sigma2_beta2(0.1), sigma2_beta3(0.1),
        pi1(0.1), tau2_1(0.01),
        pi2_cis(0.05), tau2_2_cis(0.01),
        pi3(0.1), tau2_3(0.01),
        rho_prior(0.5) {}
};

// ---------------------------------------------------------------------------
// Per-protein output
// ---------------------------------------------------------------------------
struct ProteinResult {
    std::string protein_id;
    std::string gene_name;
    int nA, nB, nC;
    int nC_exact, nC_proxy;
    int n_rf_to_pp_obs;
    bool rf_to_pp_identifiable;

    // Posterior scenario probabilities
    double prob_M0, prob_M1, prob_M2, prob_M3, prob_M4, prob_M5;
    double prob_mediator;
    double prob_protein_disease;
    double prob_rf_responsive;
    double prob_rf_direct;
    double prob_mediator_ld_resolved;
    // Deprecated compatibility alias for prob_mediator_ld_resolved.
    double prob_mediator_identified;

    // Regional shared-versus-distinct evidence. Full mode uses LD-aware
    // conditional multi-signal inference by default.
    int regional_n_variants;
    double regional_pp_shared;
    double regional_pp_distinct;
    double regional_shared_given_both;
    std::string regional_method;
    int regional_protein_signals;
    int regional_outcome_signals;
    int regional_signal_pair_count;
    int regional_independent_shared_signals;
    double regional_max_credible_set_pair_r2;
    std::vector<RegionalSignalPairResult> regional_signal_pairs;
    std::string mediation_identifiability;

    // Estimates under M=1 (true mediation)
    double beta1_est, beta1_se;
    double beta2_est, beta2_se;
    double beta3_est, beta3_se;
    double mediated_effect;     // beta1 * beta2
    double mediated_effect_se;  // delta method
    double ivw_rf_to_pp_beta;
    double ivw_rf_to_pp_se;
    double ivw_rf_to_pp_p;
    double ivw_pp_to_outcome_beta;
    double ivw_pp_to_outcome_se;
    double ivw_pp_to_outcome_p;
    double ivw_rf_to_outcome_beta;
    double ivw_rf_to_outcome_se;
    double ivw_rf_to_outcome_p;
    std::string indirect_direction;
    std::string rf_to_outcome_direction;
    std::string direction_consistent;
    double direction_consistency_prob;
    double proportion_mediated;
    double directional_mediator_prob;
    double selection_probability;
    double selection_local_fdr;
    double selection_cum_fdr;
    int selection_rank;
    double posterior_local_fdr;
    double target_local_fdr;
    double posterior_cum_fdr;
    double posterior_cum_fdr5;
    int mediation_rank;
    bool selected_fdr_10;
    bool selected_fdr_5;
    std::string evidence_tier;

    // Factorized effects are allowed to coexist. They are not normalized over
    // the legacy mutually exclusive M0-M5 scenario list.
    double factor_beta1, factor_beta1_se, factor_p_xm, factor_p_xm_strict,
           factor_log_bf_xm;
    double factor_beta1_ci_lower, factor_beta1_ci_upper;
    double factor_beta2, factor_beta2_se, factor_p_my, factor_p_my_strict,
           factor_log_bf_my;
    double factor_beta2_ci_lower, factor_beta2_ci_upper;
    double factor_beta3, factor_beta3_se, factor_p_xy, factor_p_xy_strict,
           factor_log_bf_xy;
    double factor_log_bf_heterogeneity_xm;
    double factor_log_bf_heterogeneity_my;
    double factor_log_bf_heterogeneity_xy;
    double factor_log_bf_directional_xm;
    double factor_log_bf_directional_my;
    double factor_log_bf_slope_only_xm;
    double factor_log_bf_directional_only_xm;
    double factor_log_bf_slope_directional_xm;
    double factor_log_bf_slope_only_my;
    double factor_log_bf_directional_only_my;
    double factor_log_bf_slope_directional_my;
    double factor_pp_directional_xm;
    double factor_pp_directional_my;
    double factor_directional_intercept_xm;
    double factor_directional_intercept_xm_se;
    double factor_directional_intercept_my;
    double factor_directional_intercept_my_se;
    double factor_directional_collinearity_xm;
    double factor_directional_collinearity_my;
    double factor_cross_set_max_r2;
    double factor_log_e_xm, factor_log_e_my, factor_log_e_xy;
    double factor_log_e_mediation, factor_e_q_ebh;
    // Experimental analytical calibration tracks. The balanced score assumes
    // mean-zero (InSIDE) pleiotropy. The adaptive e-value retains the oriented
    // intercept but treats the reported outcome covariance as known.
    double factor_p_xm_balanced, factor_p_my_balanced;
    double factor_balanced_conjunction_p, factor_balanced_conjunction_q_bh;
    double factor_balanced_conjunction_q_by;
    double factor_balanced_conjunction_q_adafilter;
    double factor_log_e_xm_balanced, factor_log_e_my_balanced;
    double factor_log_e_mediation_balanced, factor_e_q_balanced_ebh;
    double factor_log_e_p2e_balanced_mediation;
    double factor_e_q_p2e_balanced_ebh;
    double factor_log_e_xm_adaptive, factor_log_e_my_adaptive;
    double factor_log_e_mediation_adaptive, factor_e_q_adaptive_ebh;
    double factor_tau_xm, factor_tau_my, factor_tau_xy;
    double factor_indirect, factor_indirect_se;
    double factor_indirect_ci_lower, factor_indirect_ci_upper;
    double factor_conjunction_p, factor_conjunction_q_by;
    double factor_strict_conjunction_p, factor_strict_conjunction_q_by;
    double factor_min_log_bf;
    double factor_log_e_p2e_mediation, factor_e_q_p2e_ebh;
    double factor_pp_xm, factor_pp_my, factor_pp_two_stage;
    double factor_posterior_local_fdr, factor_posterior_cum_fdr;
    int factor_posterior_rank;
    double factor_pleiotropy_rho, factor_pleiotropy_p;
    int factor_nA, factor_nB;
    std::string factor_ld_source;
    std::string factor_selection_design;
    std::string factor_effect_estimator;
    std::string factor_pattern;
    std::string factor_two_stage_status;
    std::string factor_mediation_status;
    std::string factor_frequentist_status;
    std::string factor_strict_status;
    std::string factor_p2e_status;
    std::string factor_ebh_status;
    std::string factor_balanced_status;
    std::string factor_balanced_bh_status;
    std::string factor_adafilter_status;
    std::string factor_balanced_ebh_status;
    std::string factor_balanced_p2e_status;
    std::string factor_adaptive_ebh_status;
    std::string factor_posterior_status;

    // EB sufficient-stat approximations from the M1 fit
    double eb_beta1_second_moment;
    double eb_beta2_second_moment;
    double eb_beta3_second_moment;
    double eb_delta_pi;
    double eb_delta_second_moment;
    double eb_phi_pi;
    double eb_phi_second_moment;
    double eb_psi_pi;
    double eb_psi_second_moment;
    double eb_m3_resid_corr;

    // ELBO per scenario
    double elbo_M0, elbo_M1, elbo_M2, elbo_M3, elbo_M4, elbo_M5;

    // Flags
    bool converged;
};

// ---------------------------------------------------------------------------
// CLI options
// ---------------------------------------------------------------------------
struct Options {
    // Input files
    std::string rf_sumstat_file;       // RF GWAS summary stats
    std::string pqtl_sumstat_file;     // pQTL summary stats (stacked/pre-clumped, legacy mode)
    std::string protein_gwas_list_file; // manifest of per-protein GWAS files (full mode)
    std::string cancer_sumstat_file;   // Cancer GWAS summary stats
    std::string protein_info_file;     // Protein gene annotations
    std::string bfile_prefix;          // PLINK LD reference prefix

    // Instrument selection
    double p_thresh_rf;        // p-value threshold for RF instruments
    double p_thresh_cis;       // p-value threshold for cis-pQTL instruments
    double r2_thresh;          // LD r2 clumping threshold (instruments assumed pre-clumped)
    int    cis_window_kb;      // cis window in kb (default 1000 = +/- 1Mb)
    int    clump_window_kb;    // LD clumping window in kb
    double maf_thresh;         // MAF threshold for QC
    double freq_diff_thresh;   // Allele frequency concordance threshold
    double n_pqtl;             // pQTL sample size (for Steiger)
    double n_cancer;           // Cancer GWAS sample size (for Steiger)
    double heidi_thresh;       // HEIDI p-value threshold
    double ld_fdr_thresh;      // residual LD FDR threshold
    int    min_instruments;    // minimum instruments for analysis/QC
    bool   heidi_flag;         // run HEIDI
    bool   heidi_global;       // global or single-SNP HEIDI
    bool   steiger_flag;       // run Steiger
    bool   remove_palindromic; // remove palindromic SNPs

    // Inference
    int    max_cavi_iter;
    double elbo_tol;
    int    max_eb_iter;        // empirical Bayes outer iterations
    double eb_tol;
    bool   fixed_priors;       // keep input scenario priors fixed during EB
    bool   legacy_adaptive_priors; // opt in to data-reused soft priors/scales
    double eb_prior_strength;  // pseudo-count strength for EB scenario priors

    // Priors (can override defaults)
    double prior_p0, prior_p1, prior_p2, prior_p3, prior_p4, prior_p5;
    double prior_sigma2_beta1, prior_sigma2_beta2, prior_sigma2_beta3;
    double prior_sigma2_floor;
    double proxy_r2_thresh;
    int    ld_block_max_size;
    std::string direction_mode; // report, prioritize, soft, hard
    double direction_weight;    // soft-prior penalty strength
    double direction_min_prob;  // hard-filter consistency threshold
    int    m1_min_cis_only;     // minimum Set B instruments required for M1
    double m1_min_first_stage_z;  // minimum absolute RF->PP IVW z for M1
    double m1_min_second_stage_z; // minimum absolute PP->outcome IVW z for M1
    double m1_resid_corr_threshold; // residual corr threshold for M1 penalty
    double m1_resid_corr_penalty;   // residual corr penalty strength for M1
    double regional_prior_pp;       // per-variant protein association prior
    double regional_prior_outcome;  // per-variant outcome association prior
    double regional_prior_shared;   // per-variant shared-causal prior
    double regional_prior_var_pp;   // Wakefield effect prior variance
    double regional_prior_var_outcome;
    double regional_min_both;       // minimum P(shared or distinct)
    double regional_min_shared;     // minimum P(shared | both associated)
    std::string regional_method;     // ld-multisignal or single
    int regional_max_signals;
    double regional_signal_p;
    double regional_coverage;
    double regional_high_ld_r2;
    bool allow_unresolved_selection;
    std::string structural_method; // legacy-six-state or factorized
    double sampling_corr_rf_pqtl;
    double sampling_corr_rf_outcome;
    double sampling_corr_pqtl_outcome;
    int factor_min_set_a;
    int factor_min_set_b;
    double factor_alpha;
    double factor_bf_threshold;
    double factor_prior_xm;
    double factor_prior_my;
    double factor_prior_directional;
    double factor_directional_variance;
    int factor_quadrature_points;
    double factor_pleio_sd_xm;
    double factor_pleio_sd_my;
    double factor_pleio_sd_xy;
    int factor_pleio_quadrature_points;
    bool factor_independent_selection;

    // Output
    std::string out_prefix;
    int    threads;
    bool   verbose;

    Options() :
        p_thresh_rf(5e-6), p_thresh_cis(5e-6),
        r2_thresh(0.1), cis_window_kb(1000),
        clump_window_kb(10000), maf_thresh(0.01), freq_diff_thresh(0.2),
        n_pqtl(50000), n_cancer(100000),
        heidi_thresh(0.01), ld_fdr_thresh(0.05), min_instruments(5),
        heidi_flag(true), heidi_global(true), steiger_flag(true), remove_palindromic(true),
        max_cavi_iter(MAX_CAVI_ITER), elbo_tol(ELBO_TOL),
        max_eb_iter(20), eb_tol(1e-4), fixed_priors(true), legacy_adaptive_priors(false),
        eb_prior_strength(0.0),
        prior_p0(0.85), prior_p1(0.03), prior_p2(0.05),
        prior_p3(0.03), prior_p4(0.02), prior_p5(0.02),
        prior_sigma2_beta1(0.1), prior_sigma2_beta2(0.1), prior_sigma2_beta3(0.1),
        prior_sigma2_floor(1e-4), proxy_r2_thresh(0.8), ld_block_max_size(64),
        direction_mode("report"), direction_weight(1.0), direction_min_prob(0.80),
        m1_min_cis_only(0), m1_min_first_stage_z(0.0), m1_min_second_stage_z(0.0),
        m1_resid_corr_threshold(1.1), m1_resid_corr_penalty(0.0),
        regional_prior_pp(1e-4), regional_prior_outcome(1e-4),
        regional_prior_shared(1e-8), regional_prior_var_pp(0.04),
        regional_prior_var_outcome(0.04), regional_min_both(0.80),
        regional_min_shared(0.80), regional_method("ld-multisignal"),
        regional_max_signals(10), regional_signal_p(5e-6),
        regional_coverage(0.95), regional_high_ld_r2(0.80),
        allow_unresolved_selection(false),
        structural_method("legacy-six-state"),
        sampling_corr_rf_pqtl(0.0), sampling_corr_rf_outcome(0.0),
        sampling_corr_pqtl_outcome(0.0),
        factor_min_set_a(3), factor_min_set_b(3), factor_alpha(0.05),
        factor_bf_threshold(10.0), factor_prior_xm(0.5), factor_prior_my(0.25),
        factor_prior_directional(0.1), factor_directional_variance(0.01),
        factor_quadrature_points(161),
        factor_pleio_sd_xm(0.1), factor_pleio_sd_my(0.1),
        factor_pleio_sd_xy(0.1), factor_pleio_quadrature_points(41),
        factor_independent_selection(false),
        out_prefix("bmediator"), threads(1), verbose(false) {}
};

// ---------------------------------------------------------------------------
// Core functions (declared here, defined in .cpp files)
// ---------------------------------------------------------------------------

// IO
void parse_args(int argc, char* argv[], Options& opts);
void read_sumstats(const std::string& file, std::map<std::string, SumStat>& ss);
void read_sumstats_with_pval(const std::string& file,
                             std::map<std::string, SumStat>& ss,
                             std::map<std::string, double>& pvals,
                             bool use_selection_p = false);
void read_sumstats_with_pval_subset(const std::string& file,
                                    std::map<std::string, SumStat>& ss,
                                    std::map<std::string, double>& pvals,
                                    const std::unordered_set<std::string>& keep_rsids,
                                    bool snps_only,
                                    double min_freq,
                                    bool verbose);
void read_pqtl_sumstats(const std::string& file,
                        std::map<std::string, std::map<std::string, SumStat>>& pqtl_by_protein,
                        std::map<std::string, std::map<std::string, double>>& pqtl_pval,
                        bool use_selection_p = false);
void read_protein_info(const std::string& file,
                       std::vector<ProteinData>& proteins);
void build_instrument_sets(const std::map<std::string, SumStat>& rf_ss,
                           const std::map<std::string, SumStat>& cancer_ss,
                           const std::map<std::string, SumStat>& pqtl_ss,
                           const std::map<std::string, double>& pqtl_pval,
                           ProteinData& protein,
                           const Options& opts);

// CAVI engine
VarParams init_var_params(const ProteinData& prot, int scenario,
                          const Hyperparams& hyp);
double run_cavi(const ProteinData& prot, int scenario,
                VarParams& vp, const Hyperparams& hyp,
                const Options& opts);
void update_beta1(const ProteinData& prot, VarParams& vp,
                  const Hyperparams& hyp);
void update_beta2(const ProteinData& prot, VarParams& vp,
                  const Hyperparams& hyp);
void update_beta3(const ProteinData& prot, VarParams& vp,
                  const Hyperparams& hyp);
void update_spike_slab_delta(const ProteinData& prot, int k,
                             VarParams& vp, const Hyperparams& hyp);
void update_spike_slab_phi(const ProteinData& prot, int l,
                           VarParams& vp, const Hyperparams& hyp);
void update_spike_slab_psi(const ProteinData& prot, int k,
                           VarParams& vp, const Hyperparams& hyp);
// Scenario M5: joint bivariate update for (delta_k, psi_k) with correlation
void update_bivariate_delta_psi(const ProteinData& prot, int k,
                                VarParams& vp, const Hyperparams& hyp);
double compute_elbo(const ProteinData& prot, int scenario,
                    const VarParams& vp, const Hyperparams& hyp);

// Main analysis
ProteinResult analyze_protein(const ProteinData& prot,
                              const Hyperparams& hyp,
                              const Options& opts);
void run_factorized_inference(const ProteinData& prot,
                              ProteinResult& result,
                              const Options& opts);
void finalize_factorized_multiple_testing(std::vector<ProteinResult>& results,
                                          const Options& opts);
void run_empirical_bayes(std::vector<ProteinData>& proteins,
                         Hyperparams& hyp,
                         std::vector<ProteinResult>& results,
                         const Options& opts);

// Output
void write_results(const std::vector<ProteinResult>& results,
                   const Hyperparams& hyp,
                   const Options& opts);
void write_log(const std::string& msg, const Options& opts);

// Pipelines
void run_legacy_pipeline(const Options& opts);
void run_full_pipeline(const Options& opts);

// Utility
inline double sigmoid(double x) {
    if (x > 500.0) return 1.0;
    if (x < -500.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}

inline double log_sum_exp(double a, double b) {
    double mx = std::max(a, b);
    return mx + std::log(std::exp(a - mx) + std::exp(b - mx));
}

inline double log_sum_exp(const std::vector<double>& v) {
    double mx = *std::max_element(v.begin(), v.end());
    double s = 0.0;
    for (auto x : v) s += std::exp(x - mx);
    return mx + std::log(s);
}

inline double log_dnorm(double x, double mu, double sigma2) {
    return -0.5 * (LOG2PI + std::log(sigma2) + (x - mu) * (x - mu) / sigma2);
}

} // namespace bmediator

#endif // BMEDIATOR_H
