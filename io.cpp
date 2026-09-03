#include "bmediator.h"
#include <cstdlib>

namespace bmediator {

// ============================================================================
// Print banner (GCTA style)
// ============================================================================
static void print_banner() {
    std::cout << "*******************************************************************\n";
    std::cout << "* BMEDIATOR (Bayesian Mediation MR)\n";
    std::cout << "* Version 1.2.0-dev\n";
    std::cout << "* Bayesian framework for identifying mediating plasma proteins\n";
    std::cout << "* between risk factors and disease outcomes using summary statistics\n";
    std::cout << "*******************************************************************\n";
}

// ============================================================================
// Print usage
// ============================================================================
static void print_usage() {
    std::cout << "\nUsage: bmediator [options]\n\n";
    std::cout << "Required options:\n";
    std::cout << "  --rf-sumstat    <file>   RF GWAS summary statistics\n";
    std::cout << "  --cancer-sumstat <file>  Cancer GWAS summary statistics\n";
    std::cout << "  --protein-info  <file>   Protein gene annotation file\n";
    std::cout << "  --out           <prefix> Output file prefix\n";
    std::cout << "\nProtein input (choose one mode):\n";
    std::cout << "  --pqtl-sumstat  <file>   Pre-clumped stacked pQTL summary statistics (legacy mode)\n";
    std::cout << "  --protein-gwas-list <file> Manifest of per-protein GWAS files (full mode)\n";
    std::cout << "\nLD reference / QC:\n";
    std::cout << "  --bfile         <prefix> PLINK LD reference prefix (required for RF clumping)\n";
    std::cout << "\nInstrument selection:\n";
    std::cout << "  --p-thresh-rf   <val>    p-value threshold for RF instruments (default 5e-6)\n";
    std::cout << "  --p-thresh-cis  <val>    p-value threshold for cis-pQTL instruments (default 5e-6)\n";
    std::cout << "  --cis-window    <kb>     cis window in kb (default 1000, i.e. +/-1Mb)\n";
    std::cout << "  --clump-kb      <kb>     LD clumping window in kb (default 10000)\n";
    std::cout << "  --clump-r2      <val>    LD r-squared clumping threshold (default 0.1)\n";
    std::cout << "\nPrior specification:\n";
    std::cout << "  --prior-p0      <val>    Prior probability of null scenario (default 0.85)\n";
    std::cout << "  --prior-p1      <val>    Prior probability of mediation (default 0.03)\n";
    std::cout << "  --prior-p2      <val>    Prior probability of RF->PP only (default 0.05)\n";
    std::cout << "  --prior-p3      <val>    Prior probability of RF->Cancer direct only (default 0.03)\n";
    std::cout << "  --prior-p4      <val>    Prior probability of PP->Cancer only (default 0.02)\n";
    std::cout << "  --prior-p5      <val>    Prior probability of correlated pleiotropy (default 0.02)\n";
    std::cout << "  --sigma2-beta1  <val>    Prior variance for beta1 (RF->PP, default 0.1)\n";
    std::cout << "  --sigma2-beta2  <val>    Prior variance for beta2 (PP->Cancer, default 0.1)\n";
    std::cout << "  --sigma2-beta3  <val>    Prior variance for beta3 (RF->Cancer direct, default 0.1)\n";
    std::cout << "\nInference:\n";
    std::cout << "  --structural-method <mode> legacy-six-state or factorized (default legacy-six-state)\n";
    std::cout << "  --sampling-corr-rf-pqtl <val> Correlation of RF/pQTL estimation errors (default 0)\n";
    std::cout << "  --sampling-corr-rf-outcome <val> Correlation of RF/outcome estimation errors (default 0)\n";
    std::cout << "  --sampling-corr-pqtl-outcome <val> Correlation of pQTL/outcome estimation errors (default 0)\n";
    std::cout << "  --factor-min-set-a <int> Minimum observed Set A instruments (default 2)\n";
    std::cout << "  --factor-min-set-b <int> Minimum independent Set B instruments (default 2)\n";
    std::cout << "  --factor-alpha <val> Per-analysis and BY-FDR threshold (default 0.05)\n";
    std::cout << "  --factor-bf-threshold <val> Minimum BF required for each mediation leg (default 10)\n";
    std::cout << "  --max-cavi-iter <int>    Max CAVI iterations per scenario (default 200)\n";
    std::cout << "  --elbo-tol      <val>    ELBO convergence tolerance (default 1e-6)\n";
    std::cout << "  --max-eb-iter   <int>    Max empirical Bayes iterations when enabled (default 20)\n";
    std::cout << "  --eb-tol        <val>    EB convergence tolerance (default 1e-4)\n";
    std::cout << "  --fixed-priors           Keep all inference hyperparameters fixed (default)\n";
    std::cout << "  --empirical-bayes        Estimate priors/scales from the analyzed proteins\n";
    std::cout << "  --legacy-adaptive-priors Reuse IVW evidence in local priors (not recommended)\n";
    std::cout << "  --eb-prior-strength <val> Scenario-prior pseudo-count strength for EB (default 0)\n";
    std::cout << "  --m1-min-cis-only <int> Minimum Set B instruments required for M1 (default 0)\n";
    std::cout << "  --m1-min-second-stage-z <val> Minimum absolute PP->outcome IVW z for M1 (default 0)\n";
    std::cout << "  --m1-resid-corr-threshold <val> Residual corr threshold for M1 penalty (default disabled)\n";
    std::cout << "  --m1-resid-corr-penalty <val> M1 residual corr penalty strength (default 0)\n";
    std::cout << "  --regional-prior-pp <val> Per-variant protein association prior (default 1e-4)\n";
    std::cout << "  --regional-prior-outcome <val> Per-variant outcome association prior (default 1e-4)\n";
    std::cout << "  --regional-prior-shared <val> Per-variant shared prior (default 1e-8 = product)\n";
    std::cout << "  --regional-min-both <val> Minimum P(shared or distinct regional signal) (default 0.80)\n";
    std::cout << "  --regional-min-shared <val> Minimum P(shared|both) for LD-resolved mediation (default 0.80)\n";
    std::cout << "  --regional-method <mode> Regional model: ld-multisignal or single (default ld-multisignal)\n";
    std::cout << "  --regional-max-signals <int> Maximum conditional signals per trait (default 10)\n";
    std::cout << "  --regional-signal-p <val> Conditional signal inclusion threshold (default 5e-6)\n";
    std::cout << "  --regional-coverage <val> Credible-set coverage (default 0.95)\n";
    std::cout << "  --regional-high-ld-r2 <val> High-LD label threshold (default 0.80)\n";
    std::cout << "  --allow-unresolved-selection Allow confirmatory selection without regional LD resolution\n";
    std::cout << "\nDirection consistency:\n";
    std::cout << "  --direction-mode <mode>  Direction policy: report, prioritize, soft, hard (default report)\n";
    std::cout << "  --direction-weight <val> Soft-mode M1 log-prior penalty weight (default 1.0)\n";
    std::cout << "  --direction-min-prob <val> Hard-mode minimum Pr(direction consistent) (default 0.80)\n";
    std::cout << "\nOther:\n";
    std::cout << "  --threads       <int>    Number of threads (default 1)\n";
    std::cout << "  --verbose                Print detailed log\n";
    std::cout << "  --help                   Print this help message\n";
    std::cout << std::endl;
}

// ============================================================================
// Parse command line arguments
// ============================================================================
void parse_args(int argc, char* argv[], Options& opts) {
    print_banner();

    if (argc < 2) {
        print_usage();
        std::exit(0);
    }

    for (int i = 1; i < argc; i++) {
        std::string flag(argv[i]);

        if (flag == "--help" || flag == "-h") {
            print_usage();
            std::exit(0);
        }
        // Required
        else if (flag == "--rf-sumstat" && i + 1 < argc)
            opts.rf_sumstat_file = argv[++i];
        else if (flag == "--pqtl-sumstat" && i + 1 < argc)
            opts.pqtl_sumstat_file = argv[++i];
        else if (flag == "--protein-gwas-list" && i + 1 < argc)
            opts.protein_gwas_list_file = argv[++i];
        else if (flag == "--cancer-sumstat" && i + 1 < argc)
            opts.cancer_sumstat_file = argv[++i];
        else if (flag == "--protein-info" && i + 1 < argc)
            opts.protein_info_file = argv[++i];
        else if (flag == "--out" && i + 1 < argc)
            opts.out_prefix = argv[++i];
        else if (flag == "--bfile" && i + 1 < argc)
            opts.bfile_prefix = argv[++i];
        // Instrument selection
        else if (flag == "--p-thresh-rf" && i + 1 < argc)
            opts.p_thresh_rf = std::atof(argv[++i]);
        else if (flag == "--p-thresh-cis" && i + 1 < argc)
            opts.p_thresh_cis = std::atof(argv[++i]);
        else if (flag == "--cis-window" && i + 1 < argc)
            opts.cis_window_kb = std::atoi(argv[++i]);
        else if (flag == "--clump-kb" && i + 1 < argc)
            opts.clump_window_kb = std::atoi(argv[++i]);
        else if (flag == "--clump-r2" && i + 1 < argc)
            opts.r2_thresh = std::atof(argv[++i]);
        // Priors
        else if (flag == "--prior-p0" && i + 1 < argc)
            opts.prior_p0 = std::atof(argv[++i]);
        else if (flag == "--prior-p1" && i + 1 < argc)
            opts.prior_p1 = std::atof(argv[++i]);
        else if (flag == "--prior-p2" && i + 1 < argc)
            opts.prior_p2 = std::atof(argv[++i]);
        else if (flag == "--prior-p3" && i + 1 < argc)
            opts.prior_p3 = std::atof(argv[++i]);
        else if (flag == "--prior-p4" && i + 1 < argc)
            opts.prior_p4 = std::atof(argv[++i]);
        else if (flag == "--prior-p5" && i + 1 < argc)
            opts.prior_p5 = std::atof(argv[++i]);
        else if (flag == "--sigma2-beta1" && i + 1 < argc)
            opts.prior_sigma2_beta1 = std::atof(argv[++i]);
        else if (flag == "--sigma2-beta2" && i + 1 < argc)
            opts.prior_sigma2_beta2 = std::atof(argv[++i]);
        else if (flag == "--sigma2-beta3" && i + 1 < argc)
            opts.prior_sigma2_beta3 = std::atof(argv[++i]);
        // Inference
        else if (flag == "--structural-method" && i + 1 < argc)
            opts.structural_method = argv[++i];
        else if (flag == "--sampling-corr-rf-pqtl" && i + 1 < argc)
            opts.sampling_corr_rf_pqtl = std::atof(argv[++i]);
        else if (flag == "--sampling-corr-rf-outcome" && i + 1 < argc)
            opts.sampling_corr_rf_outcome = std::atof(argv[++i]);
        else if (flag == "--sampling-corr-pqtl-outcome" && i + 1 < argc)
            opts.sampling_corr_pqtl_outcome = std::atof(argv[++i]);
        else if (flag == "--factor-min-set-a" && i + 1 < argc)
            opts.factor_min_set_a = std::atoi(argv[++i]);
        else if (flag == "--factor-min-set-b" && i + 1 < argc)
            opts.factor_min_set_b = std::atoi(argv[++i]);
        else if (flag == "--factor-alpha" && i + 1 < argc)
            opts.factor_alpha = std::atof(argv[++i]);
        else if (flag == "--factor-bf-threshold" && i + 1 < argc)
            opts.factor_bf_threshold = std::atof(argv[++i]);
        else if (flag == "--factor-quadrature-points" && i + 1 < argc)
            opts.factor_quadrature_points = std::atoi(argv[++i]);
        else if (flag == "--max-cavi-iter" && i + 1 < argc)
            opts.max_cavi_iter = std::atoi(argv[++i]);
        else if (flag == "--elbo-tol" && i + 1 < argc)
            opts.elbo_tol = std::atof(argv[++i]);
        else if (flag == "--max-eb-iter" && i + 1 < argc)
            opts.max_eb_iter = std::atoi(argv[++i]);
        else if (flag == "--eb-tol" && i + 1 < argc)
            opts.eb_tol = std::atof(argv[++i]);
        else if (flag == "--fixed-priors")
            opts.fixed_priors = true;
        else if (flag == "--empirical-bayes")
            opts.fixed_priors = false;
        else if (flag == "--legacy-adaptive-priors")
            opts.legacy_adaptive_priors = true;
        else if (flag == "--eb-prior-strength" && i + 1 < argc)
            opts.eb_prior_strength = std::atof(argv[++i]);
        else if (flag == "--m1-min-cis-only" && i + 1 < argc)
            opts.m1_min_cis_only = std::atoi(argv[++i]);
        else if (flag == "--m1-min-first-stage-z" && i + 1 < argc)
            opts.m1_min_first_stage_z = std::atof(argv[++i]);
        else if (flag == "--m1-min-second-stage-z" && i + 1 < argc)
            opts.m1_min_second_stage_z = std::atof(argv[++i]);
        else if (flag == "--m1-resid-corr-threshold" && i + 1 < argc)
            opts.m1_resid_corr_threshold = std::atof(argv[++i]);
        else if (flag == "--m1-resid-corr-penalty" && i + 1 < argc)
            opts.m1_resid_corr_penalty = std::atof(argv[++i]);
        else if (flag == "--regional-prior-pp" && i + 1 < argc)
            opts.regional_prior_pp = std::atof(argv[++i]);
        else if (flag == "--regional-prior-outcome" && i + 1 < argc)
            opts.regional_prior_outcome = std::atof(argv[++i]);
        else if (flag == "--regional-prior-shared" && i + 1 < argc)
            opts.regional_prior_shared = std::atof(argv[++i]);
        else if (flag == "--regional-prior-var-pp" && i + 1 < argc)
            opts.regional_prior_var_pp = std::atof(argv[++i]);
        else if (flag == "--regional-prior-var-outcome" && i + 1 < argc)
            opts.regional_prior_var_outcome = std::atof(argv[++i]);
        else if (flag == "--regional-min-both" && i + 1 < argc)
            opts.regional_min_both = std::atof(argv[++i]);
        else if (flag == "--regional-min-shared" && i + 1 < argc)
            opts.regional_min_shared = std::atof(argv[++i]);
        else if (flag == "--regional-method" && i + 1 < argc)
            opts.regional_method = argv[++i];
        else if (flag == "--regional-max-signals" && i + 1 < argc)
            opts.regional_max_signals = std::atoi(argv[++i]);
        else if (flag == "--regional-signal-p" && i + 1 < argc)
            opts.regional_signal_p = std::atof(argv[++i]);
        else if (flag == "--regional-coverage" && i + 1 < argc)
            opts.regional_coverage = std::atof(argv[++i]);
        else if (flag == "--regional-high-ld-r2" && i + 1 < argc)
            opts.regional_high_ld_r2 = std::atof(argv[++i]);
        else if (flag == "--allow-unresolved-selection")
            opts.allow_unresolved_selection = true;
        // Direction consistency
        else if (flag == "--direction-mode" && i + 1 < argc)
            opts.direction_mode = argv[++i];
        else if (flag == "--direction-weight" && i + 1 < argc)
            opts.direction_weight = std::atof(argv[++i]);
        else if (flag == "--direction-min-prob" && i + 1 < argc)
            opts.direction_min_prob = std::atof(argv[++i]);
        // Other
        else if (flag == "--threads" && i + 1 < argc)
            opts.threads = std::atoi(argv[++i]);
        else if (flag == "--verbose")
            opts.verbose = true;
        // Pipeline-specific flags (parsed but stored in opts for reporting)
        else if (flag == "--maf" && i + 1 < argc)
            opts.maf_thresh = std::atof(argv[++i]);
        else if (flag == "--diff-freq" && i + 1 < argc)
            opts.freq_diff_thresh = std::atof(argv[++i]);
        else if (flag == "--n-rf" && i + 1 < argc)
            i++;
        else if (flag == "--n-pqtl" && i + 1 < argc)
            opts.n_pqtl = std::atof(argv[++i]);
        else if (flag == "--n-cancer" && i + 1 < argc)
            opts.n_cancer = std::atof(argv[++i]);
        else if (flag == "--heidi-thresh" && i + 1 < argc)
            opts.heidi_thresh = std::atof(argv[++i]);
        else if (flag == "--heidi-off")
            opts.heidi_flag = false;
        else if (flag == "--heidi-single")
            opts.heidi_global = false;
        else if (flag == "--no-steiger")
            opts.steiger_flag = false;
        else if (flag == "--no-palindromic-remove")
            opts.remove_palindromic = false;
        else if (flag == "--ld-fdr" && i + 1 < argc)
            opts.ld_fdr_thresh = std::atof(argv[++i]);
        else if (flag == "--min-instruments" && i + 1 < argc)
            opts.min_instruments = std::atoi(argv[++i]);
        else {
            std::cerr << "Error: unrecognized option '" << flag << "'\n";
            std::exit(1);
        }
    }

    // Validate required options
    bool ok = true;
    if (opts.rf_sumstat_file.empty())      { std::cerr << "Error: --rf-sumstat is required\n"; ok = false; }
    if (opts.cancer_sumstat_file.empty())   { std::cerr << "Error: --cancer-sumstat is required\n"; ok = false; }
    if (opts.protein_info_file.empty())    { std::cerr << "Error: --protein-info is required\n"; ok = false; }
    if (opts.pqtl_sumstat_file.empty() && opts.protein_gwas_list_file.empty()) {
        std::cerr << "Error: either --pqtl-sumstat or --protein-gwas-list is required\n";
        ok = false;
    }
    if (!opts.pqtl_sumstat_file.empty() && !opts.protein_gwas_list_file.empty()) {
        std::cerr << "Error: use either --pqtl-sumstat or --protein-gwas-list, not both\n";
        ok = false;
    }
    if (!(opts.p_thresh_rf > 0.0 && opts.p_thresh_rf <= 1.0) ||
        !(opts.p_thresh_cis > 0.0 && opts.p_thresh_cis <= 1.0)) {
        std::cerr << "Error: instrument p-value thresholds must be in (0, 1]\n";
        ok = false;
    }
    if (opts.cis_window_kb < 0 || opts.clump_window_kb < 0) {
        std::cerr << "Error: cis and clumping windows must be non-negative\n";
        ok = false;
    }
    if (!(opts.r2_thresh > 0.0 && opts.r2_thresh <= 1.0) ||
        !std::isfinite(opts.r2_thresh)) {
        std::cerr << "Error: --clump-r2 must be finite and in (0, 1]\n";
        ok = false;
    }
    const double scenario_priors[] = {
        opts.prior_p0, opts.prior_p1, opts.prior_p2,
        opts.prior_p3, opts.prior_p4, opts.prior_p5
    };
    double scenario_prior_sum = 0.0;
    for (double prior : scenario_priors) {
        if (!std::isfinite(prior) || prior < 0.0) {
            std::cerr << "Error: scenario priors must be finite and non-negative\n";
            ok = false;
            break;
        }
        scenario_prior_sum += prior;
    }
    if (!(scenario_prior_sum > 0.0) || !std::isfinite(scenario_prior_sum)) {
        std::cerr << "Error: scenario priors must have a positive finite sum\n";
        ok = false;
    }
    if (!(opts.prior_sigma2_beta1 > 0.0) || !std::isfinite(opts.prior_sigma2_beta1) ||
        !(opts.prior_sigma2_beta2 > 0.0) || !std::isfinite(opts.prior_sigma2_beta2) ||
        !(opts.prior_sigma2_beta3 > 0.0) || !std::isfinite(opts.prior_sigma2_beta3)) {
        std::cerr << "Error: effect prior variances must be finite and positive\n";
        ok = false;
    }
    if (opts.max_cavi_iter < 1 || opts.max_eb_iter < 1 || opts.threads < 1) {
        std::cerr << "Error: CAVI iterations, EB iterations, and threads must be positive\n";
        ok = false;
    }
    if (opts.structural_method != "legacy-six-state" &&
        opts.structural_method != "factorized") {
        std::cerr << "Error: --structural-method must be legacy-six-state or factorized\n";
        ok = false;
    }
    if (opts.structural_method == "factorized" && !opts.fixed_priors) {
        std::cerr << "Error: factorized inference requires fixed priors; "
                  << "--empirical-bayes is only available for legacy-six-state\n";
        ok = false;
    }
    const double sampling_corrs[] = {
        opts.sampling_corr_rf_pqtl, opts.sampling_corr_rf_outcome,
        opts.sampling_corr_pqtl_outcome
    };
    for (double corr : sampling_corrs) {
        if (!std::isfinite(corr) || corr <= -1.0 || corr >= 1.0) {
            std::cerr << "Error: sampling error correlations must be finite and in (-1, 1)\n";
            ok = false;
            break;
        }
    }
    if (opts.factor_min_set_a < 1 || opts.factor_min_set_b < 1) {
        std::cerr << "Error: factorized minimum instrument counts must be positive\n";
        ok = false;
    }
    if (!(opts.factor_alpha > 0.0 && opts.factor_alpha < 1.0) ||
        !std::isfinite(opts.factor_alpha)) {
        std::cerr << "Error: --factor-alpha must be finite and in (0, 1)\n";
        ok = false;
    }
    if (!(opts.factor_bf_threshold > 0.0) || !std::isfinite(opts.factor_bf_threshold)) {
        std::cerr << "Error: --factor-bf-threshold must be finite and positive\n";
        ok = false;
    }
    if (opts.factor_quadrature_points < 41) {
        std::cerr << "Error: --factor-quadrature-points must be at least 41\n";
        ok = false;
    }
    if (opts.direction_mode != "report" &&
        opts.direction_mode != "prioritize" &&
        opts.direction_mode != "soft" &&
        opts.direction_mode != "hard") {
        std::cerr << "Error: --direction-mode must be one of report, prioritize, soft, hard\n";
        ok = false;
    }
    if (!(opts.direction_weight >= 0.0) || !std::isfinite(opts.direction_weight)) {
        std::cerr << "Error: --direction-weight must be a finite non-negative value\n";
        ok = false;
    }
    if (!(opts.direction_min_prob >= 0.0 && opts.direction_min_prob <= 1.0) ||
        !std::isfinite(opts.direction_min_prob)) {
        std::cerr << "Error: --direction-min-prob must be between 0 and 1\n";
        ok = false;
    }
    if (!(opts.eb_prior_strength >= 0.0) || !std::isfinite(opts.eb_prior_strength)) {
        std::cerr << "Error: --eb-prior-strength must be a finite non-negative value\n";
        ok = false;
    }
    if (opts.m1_min_cis_only < 0) {
        std::cerr << "Error: --m1-min-cis-only must be non-negative\n";
        ok = false;
    }
    if (!(opts.m1_min_first_stage_z >= 0.0) || !std::isfinite(opts.m1_min_first_stage_z)) {
        std::cerr << "Error: --m1-min-first-stage-z must be a finite non-negative value\n";
        ok = false;
    }
    if (!(opts.m1_min_second_stage_z >= 0.0) || !std::isfinite(opts.m1_min_second_stage_z)) {
        std::cerr << "Error: --m1-min-second-stage-z must be a finite non-negative value\n";
        ok = false;
    }
    if (!(opts.m1_resid_corr_threshold >= 0.0) || !std::isfinite(opts.m1_resid_corr_threshold)) {
        std::cerr << "Error: --m1-resid-corr-threshold must be a finite non-negative value\n";
        ok = false;
    }
    if (!(opts.m1_resid_corr_penalty >= 0.0) || !std::isfinite(opts.m1_resid_corr_penalty)) {
        std::cerr << "Error: --m1-resid-corr-penalty must be a finite non-negative value\n";
        ok = false;
    }
    if (!(opts.regional_prior_pp > 0.0 && opts.regional_prior_pp < 1.0) ||
        !(opts.regional_prior_outcome > 0.0 && opts.regional_prior_outcome < 1.0) ||
        !(opts.regional_prior_shared > 0.0 && opts.regional_prior_shared < 1.0)) {
        std::cerr << "Error: regional configuration priors must be between 0 and 1\n";
        ok = false;
    }
    if (opts.regional_prior_shared > opts.regional_prior_pp ||
        opts.regional_prior_shared > opts.regional_prior_outcome) {
        std::cerr << "Error: --regional-prior-shared cannot exceed either "
                  << "trait-specific regional prior\n";
        ok = false;
    }
    if (!(opts.regional_prior_var_pp > 0.0) || !std::isfinite(opts.regional_prior_var_pp) ||
        !(opts.regional_prior_var_outcome > 0.0) ||
        !std::isfinite(opts.regional_prior_var_outcome)) {
        std::cerr << "Error: regional effect prior variances must be finite and positive\n";
        ok = false;
    }
    if (!(opts.regional_min_shared >= 0.0 && opts.regional_min_shared <= 1.0)) {
        std::cerr << "Error: --regional-min-shared must be between 0 and 1\n";
        ok = false;
    }
    if (!(opts.regional_min_both >= 0.0 && opts.regional_min_both <= 1.0)) {
        std::cerr << "Error: --regional-min-both must be between 0 and 1\n";
        ok = false;
    }
    if (opts.regional_method != "ld-multisignal" && opts.regional_method != "single") {
        std::cerr << "Error: --regional-method must be ld-multisignal or single\n";
        ok = false;
    }
    if (opts.regional_max_signals < 1) {
        std::cerr << "Error: --regional-max-signals must be positive\n";
        ok = false;
    }
    if (!(opts.regional_signal_p > 0.0 && opts.regional_signal_p <= 1.0) ||
        !std::isfinite(opts.regional_signal_p)) {
        std::cerr << "Error: --regional-signal-p must be finite and in (0, 1]\n";
        ok = false;
    }
    if (!(opts.regional_coverage > 0.0 && opts.regional_coverage <= 1.0) ||
        !std::isfinite(opts.regional_coverage)) {
        std::cerr << "Error: --regional-coverage must be finite and in (0, 1]\n";
        ok = false;
    }
    if (!(opts.regional_high_ld_r2 >= 0.0 && opts.regional_high_ld_r2 <= 1.0) ||
        !std::isfinite(opts.regional_high_ld_r2)) {
        std::cerr << "Error: --regional-high-ld-r2 must be finite and between 0 and 1\n";
        ok = false;
    }
    if (!ok) std::exit(1);

    // Normalize priors
    double psum = opts.prior_p0 + opts.prior_p1 + opts.prior_p2
                + opts.prior_p3 + opts.prior_p4 + opts.prior_p5;
    opts.prior_p0 /= psum;
    opts.prior_p1 /= psum;
    opts.prior_p2 /= psum;
    opts.prior_p3 /= psum;
    opts.prior_p4 /= psum;
    opts.prior_p5 /= psum;

    // Log options
    std::cout << "\nOptions:\n";
    std::cout << "  RF summary statistics:      " << opts.rf_sumstat_file << "\n";
    if (!opts.pqtl_sumstat_file.empty())
        std::cout << "  pQTL summary statistics:     " << opts.pqtl_sumstat_file << "\n";
    if (!opts.protein_gwas_list_file.empty())
        std::cout << "  Protein GWAS manifest:       " << opts.protein_gwas_list_file << "\n";
    std::cout << "  Cancer summary statistics:   " << opts.cancer_sumstat_file << "\n";
    std::cout << "  Protein annotation:          " << opts.protein_info_file << "\n";
    if (!opts.bfile_prefix.empty())
        std::cout << "  LD reference prefix:         " << opts.bfile_prefix << "\n";
    std::cout << "  Output prefix:               " << opts.out_prefix << "\n";
    std::cout << "  p-value threshold (RF):      " << opts.p_thresh_rf << "\n";
    std::cout << "  p-value threshold (cis):     " << opts.p_thresh_cis << "\n";
    std::cout << "  cis window (kb):             " << opts.cis_window_kb << "\n";
    std::cout << "  clump window (kb):           " << opts.clump_window_kb << "\n";
    std::cout << "  LD clumping r2:              " << opts.r2_thresh << "\n";
    std::cout << "  Structural method:           " << opts.structural_method << "\n";
    if (opts.structural_method == "factorized") {
        std::cout << "  Factorized min Set A/B:      " << opts.factor_min_set_a
                  << "/" << opts.factor_min_set_b << "\n";
        std::cout << "  Factorized alpha:            " << opts.factor_alpha << "\n";
        std::cout << "  Factorized BF threshold:     " << opts.factor_bf_threshold << "\n";
        std::cout << "  Sampling error corr (XM/XY/MY): ("
                  << opts.sampling_corr_rf_pqtl << ", "
                  << opts.sampling_corr_rf_outcome << ", "
                  << opts.sampling_corr_pqtl_outcome << ")\n";
    }
    std::cout << "  Prior (p0,p1,p2,p3,p4,p5):   ("
              << opts.prior_p0 << ", " << opts.prior_p1 << ", "
              << opts.prior_p2 << ", " << opts.prior_p3 << ", "
              << opts.prior_p4 << ", " << opts.prior_p5 << ")\n";
    std::cout << "  Max CAVI iterations:         " << opts.max_cavi_iter << "\n";
    std::cout << "  ELBO tolerance:              " << opts.elbo_tol << "\n";
    std::cout << "  Max EB iterations:           " << opts.max_eb_iter << "\n";
    std::cout << "  Fixed scenario priors:       " << (opts.fixed_priors ? "YES" : "NO") << "\n";
    std::cout << "  Legacy adaptive priors:      "
              << (opts.legacy_adaptive_priors ? "ON" : "OFF") << "\n";
    std::cout << "  EB prior strength:           " << opts.eb_prior_strength << "\n";
    std::cout << "  M1 min Set B instruments:    " << opts.m1_min_cis_only << "\n";
    std::cout << "  M1 min first-stage z:        " << opts.m1_min_first_stage_z << "\n";
    std::cout << "  M1 min second-stage z:       " << opts.m1_min_second_stage_z << "\n";
    std::cout << "  M1 residual corr penalty:    threshold="
              << opts.m1_resid_corr_threshold
              << ", weight=" << opts.m1_resid_corr_penalty << "\n";
    std::cout << "  Regional min both:           " << opts.regional_min_both << "\n";
    std::cout << "  Regional min shared:         " << opts.regional_min_shared << "\n";
    std::cout << "  Regional method:             " << opts.regional_method << "\n";
    std::cout << "  Regional max signals:        " << opts.regional_max_signals << "\n";
    std::cout << "  Regional signal p:           " << opts.regional_signal_p << "\n";
    std::cout << "  Regional credible coverage:  " << opts.regional_coverage << "\n";
    std::cout << "  Regional high-LD r2:         " << opts.regional_high_ld_r2 << "\n";
    std::cout << "  Allow unresolved selection:  "
              << (opts.allow_unresolved_selection ? "YES" : "NO") << "\n";
    std::cout << "  Direction mode:              " << opts.direction_mode << "\n";
    std::cout << "  Direction weight:            " << opts.direction_weight << "\n";
    std::cout << "  Direction min probability:   " << opts.direction_min_prob << "\n";
    std::cout << "  Threads:                     " << opts.threads << "\n";
    std::cout << "  Verbose logging:             " << (opts.verbose ? "ON" : "OFF") << "\n";
    std::cout << std::endl;
}

// ============================================================================
// Read summary statistics
// Format: SNP A1 A2 FREQ BETA SE P CHR BP
// (tab or space delimited, with header)
// ============================================================================
void read_sumstats(const std::string& file, std::map<std::string, SumStat>& ss) {
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open file " << file << "\n";
        std::exit(1);
    }
    std::string line;
    // Read and skip header
    std::getline(fin, line);

    int n = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        SumStat s;
        double pval;
        if (!(iss >> s.rsid >> s.a1 >> s.a2 >> s.freq
                  >> s.beta >> s.se >> pval >> s.chr >> s.bp)) {
            continue;  // skip malformed lines
        }
        ss[s.rsid] = s;
        n++;
    }
    fin.close();
    std::cout << "  Read " << n << " SNPs from " << file << "\n";
}

void read_sumstats_with_pval(const std::string& file,
                             std::map<std::string, SumStat>& ss,
                             std::map<std::string, double>& pvals) {
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open file " << file << "\n";
        std::exit(1);
    }

    std::string line;
    std::getline(fin, line);

    int n = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        SumStat s;
        double pval;
        if (!(iss >> s.rsid >> s.a1 >> s.a2 >> s.freq
                  >> s.beta >> s.se >> pval >> s.chr >> s.bp)) {
            continue;
        }
        for (auto& c : s.a1) c = toupper(c);
        for (auto& c : s.a2) c = toupper(c);
        ss[s.rsid] = s;
        pvals[s.rsid] = pval;
        n++;
    }
    fin.close();
    std::cout << "  Read " << n << " SNPs from " << file << "\n";
}

static bool is_simple_snp_allele(const std::string& allele) {
    if (allele.size() != 1) return false;
    char base = static_cast<char>(std::toupper(static_cast<unsigned char>(allele[0])));
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

void read_sumstats_with_pval_subset(const std::string& file,
                                    std::map<std::string, SumStat>& ss,
                                    std::map<std::string, double>& pvals,
                                    const std::unordered_set<std::string>& keep_rsids,
                                    bool snps_only,
                                    double min_freq,
                                    bool verbose) {
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open file " << file << "\n";
        std::exit(1);
    }

    std::string line;
    std::getline(fin, line);

    int n_total = 0;
    int n_kept = 0;
    int n_not_needed = 0;
    int n_non_snp = 0;
    int n_low_freq = 0;
    int n_malformed = 0;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        n_total++;

        std::istringstream iss(line);
        SumStat s;
        double pval;
        if (!(iss >> s.rsid >> s.a1 >> s.a2 >> s.freq
                  >> s.beta >> s.se >> pval >> s.chr >> s.bp)) {
            n_malformed++;
            continue;
        }

        if (!keep_rsids.empty() && !keep_rsids.count(s.rsid)) {
            n_not_needed++;
            continue;
        }

        for (auto& c : s.a1) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        for (auto& c : s.a2) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

        if (snps_only && (!is_simple_snp_allele(s.a1) || !is_simple_snp_allele(s.a2))) {
            n_non_snp++;
            continue;
        }
        if (!std::isfinite(s.freq) || s.freq <= min_freq || s.freq >= (1.0 - min_freq)) {
            n_low_freq++;
            continue;
        }

        ss[s.rsid] = s;
        pvals[s.rsid] = pval;
        n_kept++;
    }
    fin.close();

    std::cout << "  Read " << n_kept << " outcome SNPs from " << file
              << " after SNP/AF filtering and needed-SNP subsetting\n";
    if (verbose) {
        std::cout << "    Outcome rows scanned:     " << n_total << "\n";
        std::cout << "    Dropped not needed:       " << n_not_needed << "\n";
        std::cout << "    Dropped non-SNP:          " << n_non_snp << "\n";
        std::cout << "    Dropped AF <= " << min_freq
                  << " or >= " << (1.0 - min_freq) << ": " << n_low_freq << "\n";
        std::cout << "    Dropped malformed:        " << n_malformed << "\n";
    }
}

// ============================================================================
// Read pQTL summary statistics (multi-protein format)
// Format: PROTEIN SNP A1 A2 FREQ BETA SE P CHR BP
// ============================================================================
struct PqtlEntry {
    std::string protein_id;
    SumStat ss;
    double pval;
};

void read_pqtl_sumstats(const std::string& file,
                        std::map<std::string, std::map<std::string, SumStat>>& pqtl_by_protein,
                        std::map<std::string, std::map<std::string, double>>& pqtl_pval) {
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open file " << file << "\n";
        std::exit(1);
    }

    std::string line;
    std::getline(fin, line); // header

    int n = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<std::string> toks;
        std::string tok;
        while (iss >> tok) toks.push_back(tok);
        if (toks.size() < 10) {
            continue;
        }
        std::string prot_id = toks[0];
        SumStat s;
        double pval = 1.0;
        try {
            s.rsid = toks[1];
            s.a1 = toks[2];
            s.a2 = toks[3];
            s.freq = std::stod(toks[4]);
            s.beta = std::stod(toks[5]);
            s.se = std::stod(toks[6]);
            pval = std::stod(toks[7]);
            // Support both legacy 10-column format and the newer 11-column format with N.
            if (toks.size() >= 11) {
                s.chr = std::stoi(toks[9]);
                s.bp = std::stoi(toks[10]);
            } else {
                s.chr = std::stoi(toks[8]);
                s.bp = std::stoi(toks[9]);
            }
        } catch (const std::exception&) {
            continue;
        }
        for (auto& c : s.a1) c = toupper(c);
        for (auto& c : s.a2) c = toupper(c);
        pqtl_by_protein[prot_id][s.rsid] = s;
        pqtl_pval[prot_id][s.rsid] = pval;
        n++;
    }
    fin.close();
    std::cout << "  Read " << n << " pQTL entries for "
              << pqtl_by_protein.size() << " proteins from " << file << "\n";
}

void write_log(const std::string& msg, const Options& opts) {
    if (opts.verbose) {
        std::cout << msg << std::endl;
    }
}

// ============================================================================
// Read protein info
// Format: PROTEIN GENE CHR START END
// ============================================================================
void read_protein_info(const std::string& file,
                       std::vector<ProteinData>& proteins) {
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open file " << file << "\n";
        std::exit(1);
    }

    std::string line;
    std::getline(fin, line); // header

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        ProteinData p;
        if (!(iss >> p.protein_id >> p.gene_name >> p.gene_chr
                  >> p.gene_start >> p.gene_end)) {
            continue;
        }
        proteins.push_back(p);
    }
    fin.close();
    std::cout << "  Read " << proteins.size() << " proteins from " << file << "\n";
}

// ============================================================================
// Build instrument sets A, B, C for each protein
//
// For a given protein j:
//   - Identify RF instruments: SNPs with p < p_thresh_rf in RF GWAS
//   - Identify cis instruments: SNPs with p < p_thresh_cis in pQTL AND
//     within cis_window of gene
//   - Partition into:
//       Set A: RF instrument AND NOT cis
//       Set B: cis instrument AND NOT RF instrument
//       Set C: both RF instrument AND cis
// ============================================================================
void build_instrument_sets(const std::map<std::string, SumStat>& rf_ss,
                           const std::map<std::string, SumStat>& cancer_ss,
                           const std::map<std::string, SumStat>& pqtl_ss_for_protein,
                           const std::map<std::string, double>& pqtl_pval_for_protein,
                           ProteinData& protein,
                           const Options& opts) {
    // Identify RF instruments (p < threshold, approximated via beta/se -> z -> p)
    // User is expected to provide pre-clumped instruments, so we just filter on p
    std::map<std::string, bool> is_rf_instrument;
    for (auto& kv : rf_ss) {
        double z = kv.second.beta / kv.second.se;
        double p = 2.0 * std::erfc(std::fabs(z) / std::sqrt(2.0)) / 2.0;
        // More precise: p = 2 * Phi(-|z|)
        // Using erfc: p = erfc(|z|/sqrt(2))
        p = std::erfc(std::fabs(z) / std::sqrt(2.0));
        if (p < opts.p_thresh_rf) {
            is_rf_instrument[kv.first] = true;
        }
    }

    // Identify cis instruments
    int cis_window_bp = opts.cis_window_kb * 1000;
    std::map<std::string, bool> is_cis_instrument;
    for (auto& kv : pqtl_pval_for_protein) {
        if (kv.second >= opts.p_thresh_cis) continue;
        // Check if SNP is in cis region
        auto it = pqtl_ss_for_protein.find(kv.first);
        if (it == pqtl_ss_for_protein.end()) continue;
        const SumStat& s = it->second;
        if (s.chr == protein.gene_chr &&
            s.bp >= protein.gene_start - cis_window_bp &&
            s.bp <= protein.gene_end + cis_window_bp) {
            is_cis_instrument[kv.first] = true;
        }
    }

    // Partition
    // Set C: both RF and cis
    for (auto& kv : is_rf_instrument) {
        const std::string& snp = kv.first;
        if (is_cis_instrument.count(snp)) {
            // Need: gamma (RF), alpha (PP), Gamma (Cancer)
            auto rf_it = rf_ss.find(snp);
            auto pq_it = pqtl_ss_for_protein.find(snp);
            auto ca_it = cancer_ss.find(snp);
            if (rf_it == rf_ss.end() || pq_it == pqtl_ss_for_protein.end() ||
                ca_it == cancer_ss.end()) continue;

            protein.setC_rsid.push_back(snp);
            protein.setC_gamma.push_back(rf_it->second.beta);
            protein.setC_se_gamma.push_back(rf_it->second.se);
            protein.setC_alpha.push_back(pq_it->second.beta);
            protein.setC_se_alpha.push_back(pq_it->second.se);
            protein.setC_Gamma.push_back(ca_it->second.beta);
            protein.setC_se_Gamma.push_back(ca_it->second.se);
        }
    }

    // Set A: RF only (not cis)
    for (auto& kv : is_rf_instrument) {
        const std::string& snp = kv.first;
        if (is_cis_instrument.count(snp)) continue; // already in Set C

        auto rf_it = rf_ss.find(snp);
        auto pq_it = pqtl_ss_for_protein.find(snp);
        auto ca_it = cancer_ss.find(snp);
        if (rf_it == rf_ss.end() || ca_it == cancer_ss.end()) continue;

        // For Set A, we need the pQTL effect of this SNP on the protein
        // even if it's not significant — it's used in the likelihood
        double alpha_val = std::numeric_limits<double>::quiet_NaN();
        double se_alpha_val = std::numeric_limits<double>::quiet_NaN();
        bool alpha_observed = false;
        if (pq_it != pqtl_ss_for_protein.end()) {
            alpha_val = pq_it->second.beta;
            se_alpha_val = pq_it->second.se;
            alpha_observed = true;
        }

        protein.setA_rsid.push_back(snp);
        protein.setA_gamma.push_back(rf_it->second.beta);
        protein.setA_se_gamma.push_back(rf_it->second.se);
        protein.setA_alpha.push_back(alpha_val);
        protein.setA_se_alpha.push_back(se_alpha_val);
        protein.setA_alpha_observed.push_back(alpha_observed);
        protein.setA_Gamma.push_back(ca_it->second.beta);
        protein.setA_se_Gamma.push_back(ca_it->second.se);
    }

    // Set B: cis only (not RF)
    for (auto& kv : is_cis_instrument) {
        const std::string& snp = kv.first;
        if (is_rf_instrument.count(snp)) continue; // already in Set C

        auto pq_it = pqtl_ss_for_protein.find(snp);
        auto ca_it = cancer_ss.find(snp);
        if (pq_it == pqtl_ss_for_protein.end() || ca_it == cancer_ss.end())
            continue;

        protein.setB_rsid.push_back(snp);
        protein.setB_alpha_cis.push_back(pq_it->second.beta);
        protein.setB_se_alpha_cis.push_back(pq_it->second.se);
        protein.setB_Gamma_cis.push_back(ca_it->second.beta);
        protein.setB_se_Gamma_cis.push_back(ca_it->second.se);
    }
}

// ============================================================================
// Write results (GCTA-style tab-delimited output)
// ============================================================================
void write_results(const std::vector<ProteinResult>& results,
                   const Hyperparams& hyp,
                   const Options& opts) {
    // Main results file: .mediation
    std::string fname = opts.out_prefix + ".mediation";
    std::ofstream fout(fname);
    if (!fout.is_open()) {
        std::cerr << "Error: cannot open output file " << fname << "\n";
        return;
    }

    fout << std::fixed << std::setprecision(6);
    fout << "Protein\tGene\tnA\tnB\tnC\tnC_exact\tnC_proxy\tn_rf_to_pp_obs\trf_to_pp_identifiable\t"
         << "P_M0\tP_M1\tP_M2\tP_M3\tP_M4\tP_M5\t"
         << "P_mediator\tP_mediator_ld_resolved\tP_mediator_identified\t"
         << "P_protein_disease\tP_rf_responsive\tP_rf_direct\t"
         << "regional_n_variants\tregional_PP_shared\tregional_PP_distinct\t"
         << "regional_shared_given_both\tregional_method\t"
         << "regional_protein_signals\tregional_outcome_signals\t"
         << "regional_signal_pairs\tregional_max_credible_set_pair_r2\t"
         << "mediation_identifiability\t"
         << "beta1\tse_beta1\tbeta2\tse_beta2\tbeta3\tse_beta3\t"
         << "ivw_rf_to_pp_beta\tivw_rf_to_pp_se\tivw_rf_to_pp_p\t"
         << "ivw_pp_to_outcome_beta\tivw_pp_to_outcome_se\tivw_pp_to_outcome_p\t"
         << "ivw_rf_to_outcome_beta\tivw_rf_to_outcome_se\tivw_rf_to_outcome_p\t"
         << "indirect_direction\trf_to_outcome_direction\tdirection_consistent\t"
         << "direction_consistency_prob\tproportion_mediated\tdirectional_mediator_prob\t"
         << "selection_probability\tselection_local_fdr\tselection_cum_fdr\tselection_rank\t"
         << "posterior_local_fdr\ttarget_local_fdr\tposterior_cum_fdr\tposterior_cum_fdr5\tmediation_rank\tselected_fdr10\tselected_fdr5\tevidence_tier\t"
         << "factor_beta1\tfactor_beta1_se\tfactor_p_XM\tfactor_log_BF_XM\t"
         << "factor_beta2\tfactor_beta2_se\tfactor_p_MY\tfactor_log_BF_MY\t"
         << "factor_beta3\tfactor_beta3_se\tfactor_p_XY\tfactor_log_BF_XY\t"
         << "factor_indirect\tfactor_indirect_se\tfactor_conjunction_p\t"
         << "factor_conjunction_q_BY\tfactor_min_log_BF\t"
         << "factor_pleiotropy_rho\tfactor_pleiotropy_p\t"
         << "factor_nA\tfactor_nB\tfactor_ld_source\tfactor_pattern\tfactor_mediation_status\t"
         << "factor_frequentist_status\t"
         << "mediated_effect\tse_mediated\t"
         << "ELBO_M0\tELBO_M1\tELBO_M2\tELBO_M3\tELBO_M4\tELBO_M5\t"
         << "converged\n";

    for (auto& r : results) {
        fout << r.protein_id << "\t"
             << r.gene_name << "\t"
             << r.nA << "\t" << r.nB << "\t" << r.nC << "\t"
             << r.nC_exact << "\t" << r.nC_proxy << "\t"
             << r.n_rf_to_pp_obs << "\t"
             << (r.rf_to_pp_identifiable ? "YES" : "NO") << "\t"
             << r.prob_M0 << "\t" << r.prob_M1 << "\t"
             << r.prob_M2 << "\t" << r.prob_M3 << "\t"
             << r.prob_M4 << "\t" << r.prob_M5 << "\t"
             << r.prob_mediator << "\t" << r.prob_mediator_ld_resolved << "\t"
             << r.prob_mediator_identified << "\t"
             << r.prob_protein_disease << "\t"
             << r.prob_rf_responsive << "\t" << r.prob_rf_direct << "\t"
             << r.regional_n_variants << "\t"
             << r.regional_pp_shared << "\t" << r.regional_pp_distinct << "\t"
             << r.regional_shared_given_both << "\t"
             << r.regional_method << "\t"
             << r.regional_protein_signals << "\t"
             << r.regional_outcome_signals << "\t"
             << r.regional_signal_pair_count << "\t"
             << r.regional_max_credible_set_pair_r2 << "\t"
             << r.mediation_identifiability << "\t"
             << r.beta1_est << "\t" << r.beta1_se << "\t"
             << r.beta2_est << "\t" << r.beta2_se << "\t"
             << r.beta3_est << "\t" << r.beta3_se << "\t"
             << r.ivw_rf_to_pp_beta << "\t" << r.ivw_rf_to_pp_se << "\t"
             << r.ivw_rf_to_pp_p << "\t"
             << r.ivw_pp_to_outcome_beta << "\t" << r.ivw_pp_to_outcome_se << "\t"
             << r.ivw_pp_to_outcome_p << "\t"
             << r.ivw_rf_to_outcome_beta << "\t" << r.ivw_rf_to_outcome_se << "\t"
             << r.ivw_rf_to_outcome_p << "\t"
             << r.indirect_direction << "\t" << r.rf_to_outcome_direction << "\t"
             << r.direction_consistent << "\t"
             << r.direction_consistency_prob << "\t"
             << r.proportion_mediated << "\t"
             << r.directional_mediator_prob << "\t"
             << r.selection_probability << "\t"
             << r.selection_local_fdr << "\t"
             << r.selection_cum_fdr << "\t"
             << r.selection_rank << "\t"
             << r.posterior_local_fdr << "\t" << r.target_local_fdr << "\t"
             << r.posterior_cum_fdr << "\t"
             << r.posterior_cum_fdr5 << "\t"
             << r.mediation_rank << "\t"
             << (r.selected_fdr_10 ? "YES" : "NO") << "\t"
             << (r.selected_fdr_5 ? "YES" : "NO") << "\t"
             << r.evidence_tier << "\t"
             << r.factor_beta1 << "\t" << r.factor_beta1_se << "\t"
             << r.factor_p_xm << "\t" << r.factor_log_bf_xm << "\t"
             << r.factor_beta2 << "\t" << r.factor_beta2_se << "\t"
             << r.factor_p_my << "\t" << r.factor_log_bf_my << "\t"
             << r.factor_beta3 << "\t" << r.factor_beta3_se << "\t"
             << r.factor_p_xy << "\t" << r.factor_log_bf_xy << "\t"
             << r.factor_indirect << "\t" << r.factor_indirect_se << "\t"
             << r.factor_conjunction_p << "\t" << r.factor_conjunction_q_by << "\t"
             << r.factor_min_log_bf << "\t"
             << r.factor_pleiotropy_rho << "\t" << r.factor_pleiotropy_p << "\t"
             << r.factor_nA << "\t" << r.factor_nB << "\t"
             << r.factor_ld_source << "\t" << r.factor_pattern << "\t"
             << r.factor_mediation_status << "\t" << r.factor_frequentist_status << "\t"
             << r.mediated_effect << "\t" << r.mediated_effect_se << "\t"
             << r.elbo_M0 << "\t" << r.elbo_M1 << "\t"
             << r.elbo_M2 << "\t" << r.elbo_M3 << "\t"
             << r.elbo_M4 << "\t" << r.elbo_M5 << "\t"
             << (r.converged ? "YES" : "NO") << "\n";
    }
    fout.close();
    std::cout << "Results saved to " << fname << "\n";

    std::string regional_fname = opts.out_prefix + ".regional";
    std::ofstream regional_out(regional_fname);
    regional_out << std::fixed << std::setprecision(6);
    regional_out
        << "Protein\tGene\tprotein_signal\toutcome_signal\tprotein_lead\toutcome_lead\t"
        << "PP_H0\tPP_H1\tPP_H2\tPP_H3\tPP_H4\tshared_given_both\t"
        << "lead_pair_r2\tmax_credible_set_pair_r2\tinterpretation\n";
    for (const auto& result : results) {
        for (const auto& pair : result.regional_signal_pairs) {
            regional_out << result.protein_id << "\t" << result.gene_name << "\t"
                         << pair.protein_signal << "\t" << pair.outcome_signal << "\t"
                         << pair.protein_lead << "\t" << pair.outcome_lead << "\t"
                         << pair.pp_h0 << "\t" << pair.pp_h1 << "\t" << pair.pp_h2 << "\t"
                         << pair.pp_h3 << "\t" << pair.pp_h4 << "\t"
                         << pair.shared_given_both << "\t" << pair.lead_pair_r2 << "\t"
                         << pair.max_credible_set_pair_r2 << "\t"
                         << pair.interpretation << "\n";
        }
    }
    regional_out.close();
    std::cout << "Regional signal pairs saved to " << regional_fname << "\n";

    // Hyperparameter file: .hyp
    std::string hfname = opts.out_prefix + ".hyp";
    std::ofstream hout(hfname);
    hout << std::fixed << std::setprecision(12);
    hout << "# BMEDIATOR priors and final inference hyperparameters\n";
    hout << "p0\t" << hyp.p0 << "\n";
    hout << "p1\t" << hyp.p1 << "\n";
    hout << "p2\t" << hyp.p2 << "\n";
    hout << "p3\t" << hyp.p3 << "\n";
    hout << "p4\t" << hyp.p4 << "\n";
    hout << "p5\t" << hyp.p5 << "\n";
    hout << "sigma2_beta1\t" << hyp.sigma2_beta1 << "\n";
    hout << "sigma2_beta2\t" << hyp.sigma2_beta2 << "\n";
    hout << "sigma2_beta3\t" << hyp.sigma2_beta3 << "\n";
    hout << "pi1\t" << hyp.pi1 << "\n";
    hout << "tau2_1\t" << hyp.tau2_1 << "\n";
    hout << "pi2_cis\t" << hyp.pi2_cis << "\n";
    hout << "tau2_2_cis\t" << hyp.tau2_2_cis << "\n";
    hout << "pi3\t" << hyp.pi3 << "\n";
    hout << "tau2_3\t" << hyp.tau2_3 << "\n";
    hout << "fixed_priors\t" << (opts.fixed_priors ? "YES" : "NO") << "\n";
    hout << "legacy_adaptive_priors\t"
         << (opts.legacy_adaptive_priors ? "YES" : "NO") << "\n";
    hout << "eb_prior_strength\t" << opts.eb_prior_strength << "\n";
    hout << "m1_min_cis_only\t" << opts.m1_min_cis_only << "\n";
    hout << "m1_min_first_stage_z\t" << opts.m1_min_first_stage_z << "\n";
    hout << "m1_min_second_stage_z\t" << opts.m1_min_second_stage_z << "\n";
    hout << "m1_resid_corr_threshold\t" << opts.m1_resid_corr_threshold << "\n";
    hout << "m1_resid_corr_penalty\t" << opts.m1_resid_corr_penalty << "\n";
    hout << "clump_r2\t" << opts.r2_thresh << "\n";
    hout << "regional_prior_pp\t" << opts.regional_prior_pp << "\n";
    hout << "regional_prior_outcome\t" << opts.regional_prior_outcome << "\n";
    hout << "regional_prior_shared\t" << opts.regional_prior_shared << "\n";
    hout << "regional_prior_var_pp\t" << opts.regional_prior_var_pp << "\n";
    hout << "regional_prior_var_outcome\t" << opts.regional_prior_var_outcome << "\n";
    hout << "regional_min_both\t" << opts.regional_min_both << "\n";
    hout << "regional_min_shared\t" << opts.regional_min_shared << "\n";
    hout << "regional_method\t" << opts.regional_method << "\n";
    hout << "regional_max_signals\t" << opts.regional_max_signals << "\n";
    hout << "regional_signal_p\t" << opts.regional_signal_p << "\n";
    hout << "regional_coverage\t" << opts.regional_coverage << "\n";
    hout << "regional_high_ld_r2\t" << opts.regional_high_ld_r2 << "\n";
    hout << "allow_unresolved_selection\t"
         << (opts.allow_unresolved_selection ? "YES" : "NO") << "\n";
    hout << "structural_method\t" << opts.structural_method << "\n";
    hout << "sampling_corr_rf_pqtl\t" << opts.sampling_corr_rf_pqtl << "\n";
    hout << "sampling_corr_rf_outcome\t" << opts.sampling_corr_rf_outcome << "\n";
    hout << "sampling_corr_pqtl_outcome\t" << opts.sampling_corr_pqtl_outcome << "\n";
    hout << "factor_min_set_a\t" << opts.factor_min_set_a << "\n";
    hout << "factor_min_set_b\t" << opts.factor_min_set_b << "\n";
    hout << "factor_alpha\t" << opts.factor_alpha << "\n";
    hout << "factor_bf_threshold\t" << opts.factor_bf_threshold << "\n";
    hout << "factor_quadrature_points\t" << opts.factor_quadrature_points << "\n";
    hout << "direction_mode\t" << opts.direction_mode << "\n";
    hout << "direction_weight\t" << opts.direction_weight << "\n";
    hout << "direction_min_prob\t" << opts.direction_min_prob << "\n";
    hout.close();
    std::cout << "Hyperparameters saved to " << hfname << "\n";

    // Summary to log
    int n_mediators = 0;
    int n_conditionally_identified = 0;
    int n_targets = 0;
    int n_fdr10 = 0, n_fdr5 = 0;
    int n_factor_supported = 0;
    for (auto& r : results) {
        if (r.prob_M1 > 0.5) n_mediators++;
        if (r.prob_mediator_ld_resolved > 0.5) n_conditionally_identified++;
        if (r.prob_protein_disease > 0.5) n_targets++;
        if (r.selected_fdr_10) n_fdr10++;
        if (r.selected_fdr_5) n_fdr5++;
        if (r.factor_mediation_status == "SUPPORTED_CONDITIONAL") n_factor_supported++;
    }
    std::cout << "\nSummary:\n";
    std::cout << "  Total proteins analyzed:     " << results.size() << "\n";
    std::cout << "  Proteins with P(M1) > 0.5:   " << n_mediators << "\n";
    std::cout << "  LD-resolved conditional P(M1) > 0.5: "
              << n_conditionally_identified << "\n";
    std::cout << "  Proteins with P(M1+M4) > 0.5:" << n_targets << "\n";
    std::cout << "  Proteins selected at 5% mode-specific FDR:  " << n_fdr5 << "\n";
    std::cout << "  Proteins selected at 10% mode-specific FDR: " << n_fdr10 << "\n";
    std::cout << "  Direction mode:              " << opts.direction_mode << "\n";
    if (opts.structural_method == "factorized") {
        std::cout << "  Factorized conditional mediators with both BFs >= "
                  << opts.factor_bf_threshold << ": " << n_factor_supported << "\n";
    }
    std::cout << "  Final prior (p0,p1,p2,p3,p4,p5): ("
              << hyp.p0 << ", " << hyp.p1 << ", "
              << hyp.p2 << ", " << hyp.p3 << ", "
              << hyp.p4 << ", " << hyp.p5 << ")\n";
}

} // namespace bmediator
