#include "bmediator.h"

namespace bmediator {

namespace {

struct IvwSummary {
    double beta;
    double se;
    double p;
};

struct SoftPriors {
    double p0, p1, p2, p3, p4, p5;
};

struct LocalScales {
    double sigma2_beta1;
    double sigma2_beta2;
    double sigma2_beta3;
};

struct RegionalEvidence {
    int n_variants = 0;
    double pp_shared = 0.0;
    double pp_distinct = 0.0;
    double shared_given_both = 0.0;
};

inline double log_sum_exp_values(const std::vector<double>& values) {
    if (values.empty()) return -std::numeric_limits<double>::infinity();
    double max_value = *std::max_element(values.begin(), values.end());
    if (!std::isfinite(max_value)) return max_value;
    double total = 0.0;
    for (double value : values) total += std::exp(value - max_value);
    return max_value + std::log(total);
}

inline double wakefield_log_abf(double beta, double se, double prior_variance) {
    if (!std::isfinite(beta) || !std::isfinite(se) || se <= 0.0 || prior_variance <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    double variance = se * se;
    double shrinkage = prior_variance / (variance + prior_variance);
    double z = beta / se;
    return 0.5 * (std::log(1.0 - shrinkage) + shrinkage * z * z);
}

inline RegionalEvidence compute_regional_evidence(const ProteinData& prot,
                                                   const Options& opts) {
    RegionalEvidence result;
    if (prot.regional_multisignal_evaluated) {
        result.n_variants = static_cast<int>(prot.regional_cis_rsid.size());
        result.pp_shared = prot.regional_best_pp_shared;
        result.pp_distinct = prot.regional_best_pp_distinct;
        result.shared_given_both = prot.regional_best_shared_given_both;
        return result;
    }
    if (!prot.regional_data_complete) return result;

    std::vector<double> log_bf_pp;
    std::vector<double> log_bf_outcome;
    std::vector<double> log_bf_shared;
    size_t count = std::min(
        std::min(prot.regional_pp_beta.size(), prot.regional_pp_se.size()),
        std::min(prot.regional_outcome_beta.size(), prot.regional_outcome_se.size())
    );
    log_bf_pp.reserve(count);
    log_bf_outcome.reserve(count);
    log_bf_shared.reserve(count);
    for (size_t idx = 0; idx < count; ++idx) {
        double pp = wakefield_log_abf(prot.regional_pp_beta[idx], prot.regional_pp_se[idx],
                                      opts.regional_prior_var_pp);
        double outcome = wakefield_log_abf(prot.regional_outcome_beta[idx],
                                           prot.regional_outcome_se[idx],
                                           opts.regional_prior_var_outcome);
        if (!std::isfinite(pp) || !std::isfinite(outcome)) continue;
        log_bf_pp.push_back(pp);
        log_bf_outcome.push_back(outcome);
        log_bf_shared.push_back(pp + outcome);
    }
    result.n_variants = static_cast<int>(log_bf_pp.size());
    if (result.n_variants < 2) return result;

    double log_sum_pp = log_sum_exp_values(log_bf_pp);
    double log_sum_outcome = log_sum_exp_values(log_bf_outcome);
    double log_sum_shared = log_sum_exp_values(log_bf_shared);
    double log_all_pairs = log_sum_pp + log_sum_outcome;
    double ratio = std::exp(std::min(0.0, log_sum_shared - log_all_pairs));
    double log_sum_distinct = ratio >= 1.0
                                ? -std::numeric_limits<double>::infinity()
                                : log_all_pairs + std::log1p(-ratio);

    std::vector<double> log_hypotheses = {
        0.0,
        std::log(opts.regional_prior_pp) + log_sum_pp,
        std::log(opts.regional_prior_outcome) + log_sum_outcome,
        std::log(opts.regional_prior_pp) + std::log(opts.regional_prior_outcome) + log_sum_distinct,
        std::log(opts.regional_prior_shared) + log_sum_shared
    };
    double normalizer = log_sum_exp_values(log_hypotheses);
    result.pp_distinct = std::exp(log_hypotheses[3] - normalizer);
    result.pp_shared = std::exp(log_hypotheses[4] - normalizer);
    double both = result.pp_distinct + result.pp_shared;
    result.shared_given_both = both > 0.0 ? result.pp_shared / both : 0.0;
    return result;
}

inline IvwSummary compute_ivw_summary(const std::vector<double>& bx,
                                      const std::vector<double>& by,
                                      const std::vector<double>& se_by) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (bx.size() != by.size() || bx.size() != se_by.size() || bx.empty()) {
        return {nan, nan, nan};
    }

    double num = 0.0;
    double den = 0.0;
    for (size_t i = 0; i < bx.size(); ++i) {
        if (!std::isfinite(bx[i]) || !std::isfinite(by[i]) ||
            !std::isfinite(se_by[i]) || se_by[i] <= 0.0) {
            continue;
        }
        double w = 1.0 / (se_by[i] * se_by[i]);
        num += bx[i] * by[i] * w;
        den += bx[i] * bx[i] * w;
    }
    if (!(den > 0.0) || !std::isfinite(num)) {
        return {nan, nan, nan};
    }

    double beta = num / den;
    double se = std::sqrt(1.0 / den);
    if (!(se > 0.0) || !std::isfinite(se)) {
        return {beta, nan, nan};
    }
    double z = beta / se;
    double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
    return {beta, se, p};
}

inline double normal_cdf(double x) {
    if (x > 8.0) return 1.0;
    if (x < -8.0) return 0.0;
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

inline std::string effect_direction(double beta) {
    if (!std::isfinite(beta)) return "NA";
    if (beta > 0.0) return "positive";
    if (beta < 0.0) return "negative";
    return "zero";
}

inline bool has_nonzero_finite_sign(double x) {
    return std::isfinite(x) && x != 0.0;
}

inline double product_sign_consistency_prob(double beta1,
                                            double se1,
                                            double beta2,
                                            double se2,
                                            double total_beta) {
    if (!has_nonzero_finite_sign(total_beta) ||
        !std::isfinite(beta1) || !std::isfinite(beta2) ||
        !(std::isfinite(se1) && se1 > 0.0) ||
        !(std::isfinite(se2) && se2 > 0.0)) {
        return 0.5;
    }

    double p1_pos = normal_cdf(beta1 / se1);
    double p2_pos = normal_cdf(beta2 / se2);
    double p1_neg = 1.0 - p1_pos;
    double p2_neg = 1.0 - p2_pos;
    double p_product_pos = p1_pos * p2_pos + p1_neg * p2_neg;
    double p_product_neg = p1_pos * p2_neg + p1_neg * p2_pos;
    double p = (total_beta > 0.0) ? p_product_pos : p_product_neg;
    return std::min(std::max(p, 0.0), 1.0);
}

inline void populate_direction_descriptives(ProteinResult& res) {
    res.indirect_direction = effect_direction(res.mediated_effect);
    res.rf_to_outcome_direction = effect_direction(res.ivw_rf_to_outcome_beta);
    res.direction_consistency_prob = product_sign_consistency_prob(
        res.beta1_est, res.beta1_se,
        res.beta2_est, res.beta2_se,
        res.ivw_rf_to_outcome_beta
    );
    if (has_nonzero_finite_sign(res.mediated_effect) &&
        has_nonzero_finite_sign(res.ivw_rf_to_outcome_beta)) {
        bool consistent = (res.mediated_effect > 0.0) == (res.ivw_rf_to_outcome_beta > 0.0);
        res.direction_consistent = consistent ? "YES" : "NO";
    } else {
        res.direction_consistent = "NA";
    }

    if (has_nonzero_finite_sign(res.ivw_rf_to_outcome_beta) &&
        std::isfinite(res.mediated_effect)) {
        res.proportion_mediated = res.mediated_effect / res.ivw_rf_to_outcome_beta;
    } else {
        res.proportion_mediated = std::numeric_limits<double>::quiet_NaN();
    }
}

inline void finalize_direction_metrics(ProteinResult& res) {
    double directional_prob = res.prob_M1 * res.direction_consistency_prob;
    res.directional_mediator_prob = std::min(std::max(directional_prob, 0.0), 1.0);
}

inline double selection_probability_for_mode(const ProteinResult& r,
                                             const Options& opts) {
    double mediation_probability = opts.allow_unresolved_selection
                                      ? r.prob_M1 : r.prob_mediator_ld_resolved;
    if (opts.direction_mode == "prioritize") {
        return mediation_probability * r.direction_consistency_prob;
    }
    return mediation_probability;
}

inline double alpha_weight(const ProteinData& prot, int idx) {
    double reliability = 1.0;
    if (idx >= 0 && idx < prot.nA()) {
        if (idx < (int)prot.setA_alpha_reliability.size()) {
            reliability = prot.setA_alpha_reliability[idx];
        }
    } else if (idx >= prot.nA()) {
        int c = idx - prot.nA();
        if (c >= 0 && c < (int)prot.setC_alpha_reliability.size()) {
            reliability = prot.setC_alpha_reliability[c];
        }
    }
    reliability = std::max(0.05, reliability);
    if (idx >= 0 && idx < (int)prot.ld_weight_alpha_ac.size()) {
        return prot.ld_weight_alpha_ac[idx] * reliability * reliability;
    }
    return reliability * reliability;
}

inline double cancer_weight(const ProteinData& prot, int idx) {
    if (idx >= 0 && idx < (int)prot.ld_weight_cancer_union.size()) {
        return prot.ld_weight_cancer_union[idx];
    }
    return 1.0;
}

inline double shrink_prior_variance(double sigma2, double floor) {
    return std::max(sigma2, floor);
}

inline double bounded_prob(double x, double lo = 1e-4, double hi = 1.0 - 1e-4) {
    return std::min(std::max(x, lo), hi);
}

inline double clamp_value(double x, double lo, double hi) {
    return std::min(std::max(x, lo), hi);
}

inline bool scenario_beta1_free(int scenario) {
    return scenario == 1 || scenario == 2 || scenario == 5;
}

inline bool scenario_beta2_free(int scenario) {
    return scenario == 1 || scenario == 4;
}

inline bool scenario_beta3_free(int scenario) {
    return scenario == 1 || scenario == 3 || scenario == 5;
}

inline bool scenario_correlated_pleiotropy(int scenario) {
    return scenario == 5;
}

inline bool scenario_delta_free(int scenario) {
    (void)scenario;
    return true;
}

inline bool scenario_phi_free(int scenario) {
    (void)scenario;
    return true;
}

inline bool scenario_psi_free(int scenario) {
    (void)scenario;
    return true;
}

inline double stabilized_scale_update(double current,
                                      double estimate,
                                      double anchor,
                                      double floor,
                                      double max_multiple,
                                      double blend = 0.25) {
    double lower = std::max(floor, anchor * 0.25);
    double upper = std::max(lower, anchor * max_multiple);
    double clipped_estimate = clamp_value(estimate, lower, upper);
    double updated = (1.0 - blend) * current + blend * clipped_estimate;
    return clamp_value(updated, lower, upper);
}

inline SoftPriors build_soft_priors(const Hyperparams& hyp,
                                    const IvwSummary& rf_ivw,
                                    const IvwSummary& out_ivw,
                                    const IvwSummary& direct_ivw,
                                    bool adaptive) {
    if (!adaptive) {
        double total = hyp.p0 + hyp.p1 + hyp.p2 + hyp.p3 + hyp.p4 + hyp.p5;
        return {
            hyp.p0 / total, hyp.p1 / total, hyp.p2 / total,
            hyp.p3 / total, hyp.p4 / total, hyp.p5 / total
        };
    }
    double z_rf = (std::isfinite(rf_ivw.beta) && std::isfinite(rf_ivw.se) && rf_ivw.se > 0.0)
                    ? std::fabs(rf_ivw.beta / rf_ivw.se) : 0.0;
    double z_out = (std::isfinite(out_ivw.beta) && std::isfinite(out_ivw.se) && out_ivw.se > 0.0)
                    ? std::fabs(out_ivw.beta / out_ivw.se) : 0.0;
    double z_direct = (std::isfinite(direct_ivw.beta) && std::isfinite(direct_ivw.se) && direct_ivw.se > 0.0)
                    ? std::fabs(direct_ivw.beta / direct_ivw.se) : 0.0;

    // Soft evidence gates avoid hard scenario emphasis when mediator effects are weak.
    double g_rf = sigmoid((z_rf - 1.0) / 1.5);
    double g_out = sigmoid((z_out - 1.0) / 1.5);
    double g_direct = sigmoid((z_direct - 1.0) / 1.5);
    double g_both = std::sqrt(g_rf * g_out);

    SoftPriors sp;
    sp.p1 = hyp.p1 * (0.50 + 0.50 * g_both);
    sp.p2 = hyp.p2 * (0.50 + 0.50 * g_rf);
    sp.p3 = hyp.p3 * (0.50 + 0.50 * g_direct);
    sp.p4 = hyp.p4 * (0.50 + 0.50 * g_out);
    sp.p5 = hyp.p5 * (0.50 + 0.50 * std::max(g_both, g_direct));
    sp.p0 = hyp.p0;
    double s = sp.p0 + sp.p1 + sp.p2 + sp.p3 + sp.p4 + sp.p5;
    sp.p0 /= s; sp.p1 /= s; sp.p2 /= s; sp.p3 /= s; sp.p4 /= s; sp.p5 /= s;
    return sp;
}

inline LocalScales build_local_scales(const Hyperparams& hyp,
                                      const IvwSummary& rf_ivw,
                                      const IvwSummary& out_ivw,
                                      bool adaptive) {
    if (!adaptive) {
        return {hyp.sigma2_beta1, hyp.sigma2_beta2, hyp.sigma2_beta3};
    }
    double z_rf = (std::isfinite(rf_ivw.beta) && std::isfinite(rf_ivw.se) && rf_ivw.se > 0.0)
                    ? std::fabs(rf_ivw.beta / rf_ivw.se) : 0.0;
    double z_out = (std::isfinite(out_ivw.beta) && std::isfinite(out_ivw.se) && out_ivw.se > 0.0)
                    ? std::fabs(out_ivw.beta / out_ivw.se) : 0.0;
    double g_rf = sigmoid((z_rf - 1.0) / 1.5);
    double g_out = sigmoid((z_out - 1.0) / 1.5);
    LocalScales ls;
    double beta1_small = std::max(1e-4, hyp.sigma2_beta1 * 0.5);
    double beta1_large = std::max(beta1_small, hyp.sigma2_beta1 * 4.0);
    double beta2_small = std::max(1e-4, hyp.sigma2_beta2 * 0.5);
    double beta2_large = std::max(beta2_small, hyp.sigma2_beta2 * 4.0);
    ls.sigma2_beta1 = beta1_small + (beta1_large - beta1_small) * g_rf;
    ls.sigma2_beta2 = beta2_small + (beta2_large - beta2_small) * std::sqrt(g_rf * g_out);
    ls.sigma2_beta3 = hyp.sigma2_beta3;
    return ls;
}

inline void summarize_spike_slab(const std::vector<double>& omega,
                                 const std::vector<double>& mu,
                                 const std::vector<double>& s2,
                                 double& mean_pi,
                                 double& mean_second_moment) {
    if (omega.empty()) {
        mean_pi = 0.0;
        mean_second_moment = 0.0;
        return;
    }
    double s_pi = 0.0;
    double s_m2 = 0.0;
    for (size_t i = 0; i < omega.size(); ++i) {
        double oi = std::isfinite(omega[i]) ? omega[i] : 0.0;
        double mui = std::isfinite(mu[i]) ? mu[i] : 0.0;
        double s2i = std::isfinite(s2[i]) ? s2[i] : 0.0;
        s_pi += oi;
        s_m2 += oi * (s2i + mui * mui);
    }
    mean_pi = s_pi / omega.size();
    mean_second_moment = s_m2 / omega.size();
}

inline bool has_setA_alpha(const ProteinData& prot, int idx);
inline int observed_rf_to_pp_count(const ProteinData& prot);

inline double clipped_corr(double x) {
    return clamp_value(x, -0.995, 0.995);
}

inline double compute_m3_residual_corr(const ProteinData& prot,
                                       const VarParams& vp) {
    std::vector<double> ra;
    std::vector<double> rg;
    std::vector<double> ww;
    ra.reserve(prot.nA() + prot.nC());
    rg.reserve(prot.nA() + prot.nC());
    ww.reserve(prot.nA() + prot.nC());

    for (int k = 0; k < prot.nA(); ++k) {
        if (!has_setA_alpha(prot, k)) continue;
        double alpha_resid = (prot.setA_alpha[k] - vp.mu_beta1 * prot.setA_gamma[k]) /
                             std::max(prot.setA_se_alpha[k], 1e-8);
        double gamma_resid = (prot.setA_Gamma[k] - vp.mu_beta3 * prot.setA_gamma[k]) /
                             std::max(prot.setA_se_Gamma[k], 1e-8);
        ra.push_back(alpha_resid);
        rg.push_back(gamma_resid);
        ww.push_back(std::sqrt(std::max(alpha_weight(prot, k) * cancer_weight(prot, k), 1e-8)));
    }

    int offset = prot.nA();
    for (int c = 0; c < prot.nC(); ++c) {
        double alpha_resid = (prot.setC_alpha[c] - vp.mu_beta1 * prot.setC_gamma[c]) /
                             std::max(prot.setC_se_alpha[c], 1e-8);
        double gamma_resid = (prot.setC_Gamma[c] - vp.mu_beta3 * prot.setC_gamma[c]) /
                             std::max(prot.setC_se_Gamma[c], 1e-8);
        ra.push_back(alpha_resid);
        rg.push_back(gamma_resid);
        ww.push_back(std::sqrt(std::max(alpha_weight(prot, offset + c) *
                                        cancer_weight(prot, offset + c), 1e-8)));
    }

    if (ra.size() < 2) return 0.0;
    double wsum = std::accumulate(ww.begin(), ww.end(), 0.0);
    if (!(wsum > 0.0)) return 0.0;

    double ma = 0.0, mg = 0.0;
    for (size_t i = 0; i < ra.size(); ++i) {
        ma += ww[i] * ra[i];
        mg += ww[i] * rg[i];
    }
    ma /= wsum;
    mg /= wsum;

    double cov = 0.0, va = 0.0, vg = 0.0;
    for (size_t i = 0; i < ra.size(); ++i) {
        double da = ra[i] - ma;
        double dg = rg[i] - mg;
        cov += ww[i] * da * dg;
        va += ww[i] * da * da;
        vg += ww[i] * dg * dg;
    }
    if (!(va > 0.0) || !(vg > 0.0)) return 0.0;
    return clipped_corr(cov / std::sqrt(va * vg));
}

inline double compute_m3_correlation_bonus(const ProteinData& prot,
                                           const VarParams& vp) {
    int n_obs = observed_rf_to_pp_count(prot);
    if (n_obs < 2) return 0.0;
    double corr = std::fabs(compute_m3_residual_corr(prot, vp));
    double rho = clamp_value(std::fabs(vp.rho_delta_psi), 0.05, 0.995);
    double n_eff = std::min<double>(n_obs, 12.0);
    double score = 0.5 * n_eff * (2.0 * rho * corr - rho * rho);
    return clamp_value(score, -4.0, 4.0);
}

inline bool has_setA_alpha(const ProteinData& prot, int idx) {
    return idx >= 0 && idx < (int)prot.setA_alpha_observed.size() && prot.setA_alpha_observed[idx];
}

inline int observed_rf_to_pp_count(const ProteinData& prot) {
    int n = prot.nC();
    for (int k = 0; k < prot.nA(); ++k) {
        if (has_setA_alpha(prot, k)) n++;
    }
    return n;
}

} // namespace

// ============================================================================
// Initialize variational parameters for a given scenario
// ============================================================================
VarParams init_var_params(const ProteinData& prot, int scenario,
                          const Hyperparams& hyp) {
    VarParams vp;
    int nAC = prot.nA() + prot.nC();  // combined Set A + Set C count

    // Causal effects: init at zero with prior variance
    vp.mu_beta1 = 0.0;  vp.s2_beta1 = hyp.sigma2_beta1;
    vp.mu_beta2 = 0.0;  vp.s2_beta2 = hyp.sigma2_beta2;
    vp.mu_beta3 = 0.0;  vp.s2_beta3 = hyp.sigma2_beta3;
    vp.prior_sigma2_beta1 = hyp.sigma2_beta1;
    vp.prior_sigma2_beta2 = hyp.sigma2_beta2;
    vp.prior_sigma2_beta3 = hyp.sigma2_beta3;

    // Apply scenario constraints.
    if (!scenario_beta1_free(scenario)) vp.s2_beta1 = 0.0;
    if (!scenario_beta2_free(scenario)) vp.s2_beta2 = 0.0;
    if (!scenario_beta3_free(scenario)) vp.s2_beta3 = 0.0;

    // Spike-and-slab for delta (Set A + C)
    vp.mu_delta.assign(nAC, 0.0);
    vp.s2_delta.assign(nAC, hyp.tau2_1);
    vp.omega_delta.assign(nAC, hyp.pi1);

    // Spike-and-slab for phi (Set B)
    vp.mu_phi.assign(prot.nB(), 0.0);
    vp.s2_phi.assign(prot.nB(), hyp.tau2_2_cis);
    vp.omega_phi.assign(prot.nB(), hyp.pi2_cis);

    // Spike-and-slab for psi (Set A + C)
    vp.mu_psi.assign(nAC, 0.0);
    vp.s2_psi.assign(nAC, hyp.tau2_3);
    vp.omega_psi.assign(nAC, hyp.pi3);

    // Correlated-pleiotropy scenario correlation parameter
    vp.rho_delta_psi = scenario_correlated_pleiotropy(scenario) ? hyp.rho_prior : 0.0;
    vp.omega_joint.assign(nAC, hyp.pi1);  // joint inclusion for bivariate
    vp.cov_delta_psi.assign(nAC, 0.0);

    vp.elbo = -1e30;
    return vp;
}

// ============================================================================
// CAVI update for beta1 (RF -> PP)
//
// Identified from Set A and Set C via the alpha equation:
//   alpha_k = beta1 * gamma_k + delta_k + e
//
// Posterior: N(mu_beta1, s2_beta1)
// ============================================================================
void update_beta1(const ProteinData& prot, VarParams& vp,
                  const Hyperparams& hyp) {
    (void)hyp;
    double prec = 1.0 / shrink_prior_variance(vp.prior_sigma2_beta1, 1e-4);  // prior precision
    double weighted_sum = 0.0;

    // Set A contributions
    for (int k = 0; k < prot.nA(); k++) {
        if (!has_setA_alpha(prot, k)) continue;
        double se2 = prot.setA_se_alpha[k] * prot.setA_se_alpha[k];
        double gamma_k = prot.setA_gamma[k];
        double E_delta_k = vp.omega_delta[k] * vp.mu_delta[k];
        double w = alpha_weight(prot, k);

        prec += w * gamma_k * gamma_k / se2;
        weighted_sum += w * gamma_k * (prot.setA_alpha[k] - E_delta_k) / se2;
    }

    // Set C contributions (same equation structure, offset index)
    int offset = prot.nA();
    for (int c = 0; c < prot.nC(); c++) {
        double se2 = prot.setC_se_alpha[c] * prot.setC_se_alpha[c];
        double gamma_c = prot.setC_gamma[c];
        double E_delta_c = vp.omega_delta[offset + c] * vp.mu_delta[offset + c];
        double w = alpha_weight(prot, offset + c);

        prec += w * gamma_c * gamma_c / se2;
        weighted_sum += w * gamma_c * (prot.setC_alpha[c] - E_delta_c) / se2;
    }

    vp.s2_beta1 = 1.0 / prec;
    vp.mu_beta1 = vp.s2_beta1 * weighted_sum;
}

// ============================================================================
// CAVI update for beta2 (PP -> Cancer)
//
// Identified from Set B (cis instruments) and Set C (overlap instruments):
//
// Set B: Gamma_l^cis = beta2 * alpha_l^cis + phi_l + e
// Set C: Gamma_c     = beta3 * gamma_c + beta2 * alpha_c + psi_c + e
// ============================================================================
void update_beta2(const ProteinData& prot, VarParams& vp,
                  const Hyperparams& hyp) {
    (void)hyp;
    double prec = 1.0 / shrink_prior_variance(vp.prior_sigma2_beta2, 1e-4);
    double weighted_sum = 0.0;

    // Set B contributions
    for (int l = 0; l < prot.nB(); l++) {
        double se2 = prot.setB_se_Gamma_cis[l] * prot.setB_se_Gamma_cis[l];
        double alpha_l = prot.setB_alpha_cis[l];
        double E_phi_l = vp.omega_phi[l] * vp.mu_phi[l];
        double w = cancer_weight(prot, prot.nA() + prot.nC() + l);

        prec += w * alpha_l * alpha_l / se2;
        weighted_sum += w * alpha_l * (prot.setB_Gamma_cis[l] - E_phi_l) / se2;
    }

    // Set C contributions
    int offset = prot.nA();
    for (int c = 0; c < prot.nC(); c++) {
        double se2 = prot.setC_se_Gamma[c] * prot.setC_se_Gamma[c];
        double alpha_c = prot.setC_alpha[c];
        double E_beta3 = vp.mu_beta3;
        double gamma_c = prot.setC_gamma[c];
        double E_psi_c = vp.omega_psi[offset + c] * vp.mu_psi[offset + c];
        double w = cancer_weight(prot, offset + c);

        prec += w * alpha_c * alpha_c / se2;
        weighted_sum += w * alpha_c * (prot.setC_Gamma[c] - E_beta3 * gamma_c - E_psi_c) / se2;
    }

    vp.s2_beta2 = 1.0 / prec;
    vp.mu_beta2 = vp.s2_beta2 * weighted_sum;
}

// ============================================================================
// CAVI update for beta3 (RF -> Cancer direct)
//
// Set A: Gamma_k = (beta3 + beta1*beta2)*gamma_k + beta2*delta_k + psi_k + e
// Set C: Gamma_c = beta3*gamma_c + beta2*alpha_c + psi_c + e
// ============================================================================
void update_beta3(const ProteinData& prot, VarParams& vp,
                  const Hyperparams& hyp) {
    (void)hyp;
    double prec = 1.0 / shrink_prior_variance(vp.prior_sigma2_beta3, 1e-4);
    double weighted_sum = 0.0;

    double E_b1b2 = vp.mu_beta1 * vp.mu_beta2;  // mean-field factorization
    double E_beta2 = vp.mu_beta2;

    // Set A contributions
    for (int k = 0; k < prot.nA(); k++) {
        double se2 = prot.setA_se_Gamma[k] * prot.setA_se_Gamma[k];
        double gamma_k = prot.setA_gamma[k];
        double E_delta_k = vp.omega_delta[k] * vp.mu_delta[k];
        double E_psi_k = vp.omega_psi[k] * vp.mu_psi[k];
        double w = cancer_weight(prot, k);

        prec += w * gamma_k * gamma_k / se2;
        double residual = prot.setA_Gamma[k]
                        - E_b1b2 * gamma_k
                        - E_beta2 * E_delta_k
                        - E_psi_k;
        weighted_sum += w * gamma_k * residual / se2;
    }

    // Set C contributions
    int offset = prot.nA();
    for (int c = 0; c < prot.nC(); c++) {
        double se2 = prot.setC_se_Gamma[c] * prot.setC_se_Gamma[c];
        double gamma_c = prot.setC_gamma[c];
        double alpha_c = prot.setC_alpha[c];
        double E_psi_c = vp.omega_psi[offset + c] * vp.mu_psi[offset + c];
        double w = cancer_weight(prot, offset + c);

        prec += w * gamma_c * gamma_c / se2;
        double residual = prot.setC_Gamma[c]
                        - E_beta2 * alpha_c
                        - E_psi_c;
        weighted_sum += w * gamma_c * residual / se2;
    }

    vp.s2_beta3 = 1.0 / prec;
    vp.mu_beta3 = vp.s2_beta3 * weighted_sum;
}

// ============================================================================
// Spike-and-slab update for delta_k (stage 1 pleiotropy, Set A + C)
//
// Residual from alpha equation: r_k = alpha_k - beta1*gamma_k
// delta_k ~ (1-omega)*delta_0 + omega*N(mu, s2)
// ============================================================================
void update_spike_slab_delta(const ProteinData& prot, int idx,
                             VarParams& vp, const Hyperparams& hyp) {
    // idx runs 0..nA-1 for Set A, nA..nA+nC-1 for Set C
    double alpha_obs, se_alpha, gamma_k;

    if (idx < prot.nA()) {
        if (!has_setA_alpha(prot, idx)) {
            vp.omega_delta[idx] = hyp.pi1;
            vp.mu_delta[idx] = 0.0;
            vp.s2_delta[idx] = hyp.tau2_1;
            return;
        }
        alpha_obs = prot.setA_alpha[idx];
        se_alpha  = prot.setA_se_alpha[idx];
        gamma_k   = prot.setA_gamma[idx];
    } else {
        int c = idx - prot.nA();
        alpha_obs = prot.setC_alpha[c];
        se_alpha  = prot.setC_se_alpha[c];
        gamma_k   = prot.setC_gamma[c];
    }

    double se2 = se_alpha * se_alpha;
    double residual = alpha_obs - vp.mu_beta1 * gamma_k;
    double w = alpha_weight(prot, idx);

    // Slab posterior
    double prec_slab = 1.0 / hyp.tau2_1 + w / se2;
    double s2 = 1.0 / prec_slab;
    double mu = s2 * (w * residual / se2);

    // Log Bayes factor for slab vs spike
    double log_bf = 0.5 * std::log(s2 / hyp.tau2_1)
                  + 0.5 * mu * mu / s2;

    // Posterior inclusion probability
    double log_odds_prior = std::log(hyp.pi1 / (1.0 - hyp.pi1 + TINY));
    vp.omega_delta[idx] = sigmoid(log_odds_prior + log_bf);
    vp.mu_delta[idx] = mu;
    vp.s2_delta[idx] = s2;
}

// ============================================================================
// Spike-and-slab update for phi_l (cis pleiotropy, Set B)
//
// Residual from cis cancer equation: r_l = Gamma_l - beta2*alpha_l
// ============================================================================
void update_spike_slab_phi(const ProteinData& prot, int l,
                           VarParams& vp, const Hyperparams& hyp) {
    double se2 = prot.setB_se_Gamma_cis[l] * prot.setB_se_Gamma_cis[l];
    double residual = prot.setB_Gamma_cis[l]
                    - vp.mu_beta2 * prot.setB_alpha_cis[l];
    double w = cancer_weight(prot, prot.nA() + prot.nC() + l);

    double prec_slab = 1.0 / hyp.tau2_2_cis + w / se2;
    double s2 = 1.0 / prec_slab;
    double mu = s2 * (w * residual / se2);

    double log_bf = 0.5 * std::log(s2 / hyp.tau2_2_cis)
                  + 0.5 * mu * mu / s2;

    double log_odds_prior = std::log(hyp.pi2_cis / (1.0 - hyp.pi2_cis + TINY));
    vp.omega_phi[l] = sigmoid(log_odds_prior + log_bf);
    vp.mu_phi[l] = mu;
    vp.s2_phi[l] = s2;
}

// ============================================================================
// Spike-and-slab update for psi_k (cancer pleiotropy, Set A + C)
//
// Set A residual: r_k = Gamma_k - (beta3+beta1*beta2)*gamma_k - beta2*delta_k
// Set C residual: r_c = Gamma_c - beta3*gamma_c - beta2*alpha_c
// ============================================================================
void update_spike_slab_psi(const ProteinData& prot, int idx,
                           VarParams& vp, const Hyperparams& hyp) {
    double se_Gamma, residual;

    double E_beta2 = vp.mu_beta2;
    double E_beta3 = vp.mu_beta3;
    double E_b1b2  = vp.mu_beta1 * vp.mu_beta2;

    if (idx < prot.nA()) {
        se_Gamma = prot.setA_se_Gamma[idx];
        double gamma_k = prot.setA_gamma[idx];
        double E_delta_k = vp.omega_delta[idx] * vp.mu_delta[idx];
        residual = prot.setA_Gamma[idx]
                 - (E_beta3 + E_b1b2) * gamma_k
                 - E_beta2 * E_delta_k;
    } else {
        int c = idx - prot.nA();
        se_Gamma = prot.setC_se_Gamma[c];
        double gamma_c = prot.setC_gamma[c];
        double alpha_c = prot.setC_alpha[c];
        residual = prot.setC_Gamma[c]
                 - E_beta3 * gamma_c
                 - E_beta2 * alpha_c;
    }

    double se2 = se_Gamma * se_Gamma;
    double w = cancer_weight(prot, idx);
    double prec_slab = 1.0 / hyp.tau2_3 + w / se2;
    double s2 = 1.0 / prec_slab;
    double mu = s2 * (w * residual / se2);

    double log_bf = 0.5 * std::log(s2 / hyp.tau2_3)
                  + 0.5 * mu * mu / s2;

    double log_odds_prior = std::log(hyp.pi3 / (1.0 - hyp.pi3 + TINY));
    vp.omega_psi[idx] = sigmoid(log_odds_prior + log_bf);
    vp.mu_psi[idx] = mu;
    vp.s2_psi[idx] = s2;
}

// ============================================================================
// Scenario M5: Bivariate spike-and-slab update for (delta_k, psi_k)
//
// Under M=3, beta2=0. The two observation equations for Set A instrument k:
//
//   alpha_k = beta1*gamma_k + delta_k + e_alpha    (e_alpha ~ N(0, se_alpha^2))
//   Gamma_k = beta3*gamma_k         + psi_k + e_Gamma  (e_Gamma ~ N(0, se_Gamma^2))
//
// Prior on (delta_k, psi_k) | z_k=1:
//   N_2(0, Sigma_prior)  where Sigma_prior = [[tau1^2,    rho*tau1*tau3],
//                                              [rho*tau1*tau3, tau3^2     ]]
//
// Residuals:
//   r_alpha = alpha_k - beta1*gamma_k
//   r_Gamma = Gamma_k - beta3*gamma_k
//
// Likelihood for (delta_k, psi_k):
//   r_alpha ~ N(delta_k, se_alpha^2)
//   r_Gamma ~ N(psi_k,   se_Gamma^2)
//   These are independent given (delta, psi) because alpha and Gamma have
//   independent measurement errors.
//
// Posterior (bivariate normal):
//   Precision_post = Sigma_prior^{-1} + diag(1/se_alpha^2, 1/se_Gamma^2)
//   mu_post = Precision_post^{-1} * [r_alpha/se_alpha^2, r_Gamma/se_Gamma^2]'
//
// Log Bayes factor for slab vs spike:
//   log BF = 0.5 * log(|Sigma_post| / |Sigma_prior|)
//          + 0.5 * mu_post' * Precision_post * mu_post
//
// The key difference from independent updates: the posterior mu_delta and
// mu_psi are coupled through the off-diagonal of the prior precision,
// which propagates information from the cancer residual to the delta estimate
// and vice versa. This is exactly how correlated pleiotropy "leaks" signal.
// ============================================================================
void update_bivariate_delta_psi(const ProteinData& prot, int idx,
                                VarParams& vp, const Hyperparams& hyp) {
    // Get residuals
    double r_alpha, se_alpha, r_Gamma, se_Gamma;

    if (idx < prot.nA()) {
        se_Gamma = prot.setA_se_Gamma[idx];
        r_Gamma = prot.setA_Gamma[idx] - vp.mu_beta3 * prot.setA_gamma[idx];
        if (!has_setA_alpha(prot, idx)) {
            double sg2 = se_Gamma * se_Gamma;
            double t1 = hyp.tau2_1;
            double t3 = hyp.tau2_3;
            double rho = vp.rho_delta_psi;
            double cov_prior = rho * std::sqrt(t1 * t3);
            double det_prior = t1 * t3 - cov_prior * cov_prior;
            if (det_prior < 1e-20) det_prior = 1e-20;
            double P22_prior = t1 / det_prior;
            double P22_post = P22_prior + 1.0 / sg2;
            double S22_post = 1.0 / P22_post;
            double mu_psi_post = S22_post * (r_Gamma / sg2);
            double log_bf = 0.5 * std::log(S22_post / t3)
                          + 0.5 * mu_psi_post * mu_psi_post / S22_post;
            double pi_joint = std::max(hyp.pi1, hyp.pi3);
            double log_odds_prior = std::log(pi_joint / (1.0 - pi_joint + TINY));
            double omega = sigmoid(log_odds_prior + log_bf);
            vp.omega_joint[idx] = omega;
            vp.mu_delta[idx] = 0.0;
            vp.s2_delta[idx] = hyp.tau2_1;
            vp.mu_psi[idx] = mu_psi_post;
            vp.s2_psi[idx] = S22_post;
            vp.omega_psi[idx] = omega;
            vp.cov_delta_psi[idx] = 0.0;
            return;
        }
        se_alpha = prot.setA_se_alpha[idx];
        r_alpha = prot.setA_alpha[idx] - vp.mu_beta1 * prot.setA_gamma[idx];
    } else {
        int c = idx - prot.nA();
        se_alpha = prot.setC_se_alpha[c];
        r_alpha = prot.setC_alpha[c] - vp.mu_beta1 * prot.setC_gamma[c];
        se_Gamma = prot.setC_se_Gamma[c];
        // For Set C under M=3 (beta2=0): Gamma = beta3*gamma + psi
        r_Gamma = prot.setC_Gamma[c] - vp.mu_beta3 * prot.setC_gamma[c];
    }

    double sa2 = se_alpha * se_alpha;
    double sg2 = se_Gamma * se_Gamma;

    // Prior covariance matrix and its inverse
    double t1 = hyp.tau2_1;
    double t3 = hyp.tau2_3;
    double rho = vp.rho_delta_psi;
    double cov_prior = rho * std::sqrt(t1 * t3);

    // Sigma_prior = [[t1, cov_prior], [cov_prior, t3]]
    // det(Sigma_prior) = t1*t3 - cov_prior^2 = t1*t3*(1-rho^2)
    double det_prior = t1 * t3 - cov_prior * cov_prior;
    if (det_prior < 1e-20) det_prior = 1e-20;  // numerical safety

    // Sigma_prior^{-1}
    double P11_prior = t3 / det_prior;
    double P22_prior = t1 / det_prior;
    double P12_prior = -cov_prior / det_prior;

    // Posterior precision = prior precision + data precision
    double P11_post = P11_prior + 1.0 / sa2;
    double P22_post = P22_prior + 1.0 / sg2;
    double P12_post = P12_prior;  // data precisions are diagonal (no cross-term)

    // Posterior covariance (invert 2x2)
    double det_post = P11_post * P22_post - P12_post * P12_post;
    if (det_post < 1e-20) det_post = 1e-20;

    double S11_post = P22_post / det_post;
    double S22_post = P11_post / det_post;
    double S12_post = -P12_post / det_post;

    // Posterior mean
    double b1 = r_alpha / sa2;  // data contribution to delta
    double b2 = r_Gamma / sg2;  // data contribution to psi
    double mu_delta_post = S11_post * b1 + S12_post * b2;
    double mu_psi_post   = S12_post * b1 + S22_post * b2;

    // Log Bayes factor for bivariate slab vs spike
    // BF = p(r_alpha, r_Gamma | slab) / p(r_alpha, r_Gamma | spike)
    //
    // Under spike: delta=psi=0, so
    //   log p(data|spike) = -0.5*[r_alpha^2/sa2 + r_Gamma^2/sg2 + log(sa2) + log(sg2) + 2*log(2pi)]
    //
    // Under slab: marginalize over (delta, psi) with bivariate normal prior
    //   log p(data|slab) = log p(data|spike) + log BF
    //
    // log BF = 0.5 * [log(det_post^{-1}) - log(det_prior)] + 0.5 * mu_post' * P_post * mu_post
    //        = 0.5 * [-log(det_post) - log(det_prior)] + 0.5 * quadratic form
    //
    // Actually more precisely:
    // log BF = 0.5 * log(det(Sigma_post) / det(Sigma_prior))
    //        + 0.5 * mu_post' * Sigma_post^{-1} * mu_post
    //
    // det(Sigma_post) = 1/det_post, det(Sigma_prior) = det_prior
    double log_det_ratio = std::log(1.0 / (det_post * det_prior) + TINY);
    double quad_form = P11_post * mu_delta_post * mu_delta_post
                     + P22_post * mu_psi_post * mu_psi_post
                     + 2.0 * P12_post * mu_delta_post * mu_psi_post;
    double log_bf = 0.5 * log_det_ratio + 0.5 * quad_form;

    // Joint inclusion probability
    // Use a joint pi: probability that both delta_k and psi_k are non-null
    // Under M=3, we use a shared pi for the joint pair
    double pi_joint = std::max(hyp.pi1, hyp.pi3);  // or could be a separate hyp
    double log_odds_prior = std::log(pi_joint / (1.0 - pi_joint + TINY));
    double omega = sigmoid(log_odds_prior + log_bf);

    // Store
    vp.omega_joint[idx] = omega;
    vp.omega_delta[idx] = omega;  // marginal inclusion for delta = joint inclusion
    vp.omega_psi[idx]   = omega;  // marginal inclusion for psi = joint inclusion
    vp.mu_delta[idx]    = mu_delta_post;
    vp.s2_delta[idx]    = S11_post;
    vp.mu_psi[idx]      = mu_psi_post;
    vp.s2_psi[idx]      = S22_post;
    vp.cov_delta_psi[idx] = S12_post;
}

// ============================================================================
// Compute ELBO
// ============================================================================
double compute_elbo(const ProteinData& prot, int scenario,
                    const VarParams& vp, const Hyperparams& hyp) {
    double elbo = 0.0;
    int nAC = prot.nA() + prot.nC();

    bool beta1_free = scenario_beta1_free(scenario);
    bool beta2_free = scenario_beta2_free(scenario);
    bool beta3_free = scenario_beta3_free(scenario);
    bool delta_free = scenario_delta_free(scenario);
    bool phi_free = scenario_phi_free(scenario);
    bool psi_free = scenario_psi_free(scenario);

    // --- Expected log-likelihood contributions ---

    // Set A: alpha equation
    // alpha_k = beta1 * gamma_k + delta_k + e
    // Under M=0: beta1=0, no delta -> alpha_k = e (pure noise)
    for (int k = 0; k < prot.nA(); k++) {
        if (!has_setA_alpha(prot, k)) continue;
        double se2 = prot.setA_se_alpha[k] * prot.setA_se_alpha[k];
        double E_delta = delta_free ? (vp.omega_delta[k] * vp.mu_delta[k]) : 0.0;
        double r = prot.setA_alpha[k] - vp.mu_beta1 * prot.setA_gamma[k] - E_delta;
        double var_term = 0.0;
        double w = alpha_weight(prot, k);
        if (delta_free) {
            var_term += vp.s2_beta1 * prot.setA_gamma[k] * prot.setA_gamma[k];
            var_term += vp.omega_delta[k] * (vp.s2_delta[k] + vp.mu_delta[k] * vp.mu_delta[k])
                      - E_delta * E_delta;
        }
        elbo += -0.5 * w * (LOG2PI + std::log(se2) + (r * r + var_term) / se2);
    }

    // Set A: cancer equation
    // Gamma_k = (beta3 + beta1*beta2)*gamma_k + beta2*delta_k + psi_k + e
    // Under M=0: Gamma_k = beta3*gamma_k + psi_k + e
    // Under M=2,3: beta2=0 -> Gamma_k = (beta3)*gamma_k + psi_k + e (but beta1 is free)
    double E_b1b2 = vp.mu_beta1 * vp.mu_beta2;
    for (int k = 0; k < prot.nA(); k++) {
        double se2 = prot.setA_se_Gamma[k] * prot.setA_se_Gamma[k];
        double gk = prot.setA_gamma[k];
        double E_delta = delta_free ? (vp.omega_delta[k] * vp.mu_delta[k]) : 0.0;
        double E_psi = psi_free ? (vp.omega_psi[k] * vp.mu_psi[k]) : 0.0;
        double w = cancer_weight(prot, k);

        double mean_term;
        if (beta2_free) {
            mean_term = ((beta3_free ? vp.mu_beta3 : 0.0) + E_b1b2) * gk
                      + vp.mu_beta2 * E_delta + E_psi;
        } else {
            // beta2 = 0: no mediation term, no delta*beta2 term
            mean_term = (beta3_free ? vp.mu_beta3 : 0.0) * gk + E_psi;
        }
        double r = prot.setA_Gamma[k] - mean_term;
        double var_term = 0.0;
        if (beta3_free) {
            var_term += vp.s2_beta3 * gk * gk;
        }
        if (psi_free) {
            var_term += vp.omega_psi[k] * (vp.s2_psi[k] + vp.mu_psi[k] * vp.mu_psi[k])
                      - E_psi * E_psi;
        }
        elbo += -0.5 * w * (LOG2PI + std::log(se2) + (r * r + var_term) / se2);
    }

    // Set C: alpha equation (same structure as Set A alpha)
    int offset = prot.nA();
    for (int c = 0; c < prot.nC(); c++) {
        double se2 = prot.setC_se_alpha[c] * prot.setC_se_alpha[c];
        double E_delta = delta_free ? (vp.omega_delta[offset + c] * vp.mu_delta[offset + c]) : 0.0;
        double r = prot.setC_alpha[c] - vp.mu_beta1 * prot.setC_gamma[c] - E_delta;
        double var_term = 0.0;
        double w = alpha_weight(prot, offset + c);
        if (delta_free) {
            var_term += vp.s2_beta1 * prot.setC_gamma[c] * prot.setC_gamma[c];
            var_term += vp.omega_delta[offset + c] * (vp.s2_delta[offset + c]
                        + vp.mu_delta[offset + c] * vp.mu_delta[offset + c])
                      - E_delta * E_delta;
        }
        elbo += -0.5 * w * (LOG2PI + std::log(se2) + (r * r + var_term) / se2);
    }

    // Set C: cancer equation
    for (int c = 0; c < prot.nC(); c++) {
        double se2 = prot.setC_se_Gamma[c] * prot.setC_se_Gamma[c];
        double gc = prot.setC_gamma[c];
        double ac = prot.setC_alpha[c];
        double E_psi = psi_free ? (vp.omega_psi[offset + c] * vp.mu_psi[offset + c]) : 0.0;
        double w = cancer_weight(prot, offset + c);
        double mean_term;
        if (beta2_free) {
            mean_term = (beta3_free ? vp.mu_beta3 : 0.0) * gc + vp.mu_beta2 * ac + E_psi;
        } else {
            // beta2=0: for Set C, cancer equation becomes beta3*gamma + psi
            mean_term = (beta3_free ? vp.mu_beta3 : 0.0) * gc + E_psi;
        }
        double r = prot.setC_Gamma[c] - mean_term;
        double var_term = 0.0;
        if (beta3_free) {
            var_term += vp.s2_beta3 * gc * gc;
        }
        if (psi_free) {
            var_term += vp.omega_psi[offset + c] * (vp.s2_psi[offset + c]
                          + vp.mu_psi[offset + c] * vp.mu_psi[offset + c])
                        - E_psi * E_psi;
        }
        if (beta2_free) {
            var_term += vp.s2_beta2 * ac * ac;
        }
        elbo += -0.5 * w * (LOG2PI + std::log(se2) + (r * r + var_term) / se2);
    }

    // Set B: cancer equation
    // Gamma_l = beta2 * alpha_l + phi_l + e
    // Under M=0,2,3: beta2=0 -> Gamma_l = phi_l + e (or pure noise if phi not modeled)
    for (int l = 0; l < prot.nB(); l++) {
        double se2 = prot.setB_se_Gamma_cis[l] * prot.setB_se_Gamma_cis[l];
        double al = prot.setB_alpha_cis[l];
        double E_phi = phi_free ? (vp.omega_phi[l] * vp.mu_phi[l]) : 0.0;
        double w = cancer_weight(prot, prot.nA() + prot.nC() + l);
        double mean_term;
        if (beta2_free) {
            mean_term = vp.mu_beta2 * al + E_phi;
        } else {
            mean_term = 0.0;  // pure noise: cancer effect is just noise
        }
        double r = prot.setB_Gamma_cis[l] - mean_term;
        double var_term = 0.0;
        if (phi_free) {
            var_term += vp.s2_beta2 * al * al;
            var_term += vp.omega_phi[l] * (vp.s2_phi[l] + vp.mu_phi[l] * vp.mu_phi[l])
                      - E_phi * E_phi;
        }
        elbo += -0.5 * w * (LOG2PI + std::log(se2) + (r * r + var_term) / se2);
    }

    // --- KL divergences: -KL(q || prior) = E_q[log p(theta)] - E_q[log q(theta)] ---

    // KL for beta1
    if (beta1_free) {
        double prior_var = shrink_prior_variance(vp.prior_sigma2_beta1, 1e-4);
        // -KL(N(mu,s2) || N(0, sigma2))
        elbo += -0.5 * (vp.mu_beta1 * vp.mu_beta1 + vp.s2_beta1) / prior_var
                + 0.5 * std::log(vp.s2_beta1 / prior_var + TINY) + 0.5;
    }

    // KL for beta2
    if (beta2_free) {
        double prior_var = shrink_prior_variance(vp.prior_sigma2_beta2, 1e-4);
        elbo += -0.5 * (vp.mu_beta2 * vp.mu_beta2 + vp.s2_beta2) / prior_var
                + 0.5 * std::log(vp.s2_beta2 / prior_var + TINY) + 0.5;
    }

    // KL for beta3
    if (beta3_free) {
        double prior_var = shrink_prior_variance(vp.prior_sigma2_beta3, 1e-4);
        elbo += -0.5 * (vp.mu_beta3 * vp.mu_beta3 + vp.s2_beta3) / prior_var
                + 0.5 * std::log(vp.s2_beta3 / prior_var + TINY) + 0.5;
    }

    // KL for spike-and-slab terms (Bernoulli-Gaussian)
    auto kl_spike_slab = [](double omega, double mu, double s2,
                            double pi_prior, double tau2) -> double {
        double kl = 0.0;
        if (omega > EPS && omega < 1.0 - EPS) {
            kl += omega * std::log(omega / (pi_prior + TINY))
                + (1.0 - omega) * std::log((1.0 - omega) / (1.0 - pi_prior + TINY));
            kl += omega * 0.5 * ((mu * mu + s2) / tau2
                                 - 1.0 - std::log(s2 / tau2 + TINY));
        }
        return -kl;
    };

    // delta and psi KL divergences
    if (scenario_correlated_pleiotropy(scenario)) {
        // Correlated pleiotropy: bivariate KL for joint (delta_k, psi_k)
        // -KL(q(delta,psi,z) || p(delta,psi,z))
        //
        // q: Bernoulli(omega_joint) * N_2(mu_post, Sigma_post)
        // p: Bernoulli(pi_joint) * N_2(0, Sigma_prior)
        //
        // -KL = -KL_Bernoulli - omega * KL_bivariate_normal
        double t1 = hyp.tau2_1;
        double t3 = hyp.tau2_3;
        double rho = vp.rho_delta_psi;
        double cov_prior = rho * std::sqrt(t1 * t3);
        double det_prior = t1 * t3 - cov_prior * cov_prior;
        if (det_prior < 1e-20) det_prior = 1e-20;

        double pi_joint = std::max(hyp.pi1, hyp.pi3);

        for (int k = 0; k < nAC; k++) {
            double omega = vp.omega_joint[k];
            if (omega < EPS || omega > 1.0 - EPS) continue;

            // Bernoulli KL
            double kl_bern = omega * std::log(omega / (pi_joint + TINY))
                           + (1.0 - omega) * std::log((1.0 - omega) / (1.0 - pi_joint + TINY));

            // Bivariate normal KL: KL(N(mu_post, Sigma_post) || N(0, Sigma_prior))
            // = 0.5 * [tr(Sigma_prior^{-1} Sigma_post) + mu' Sigma_prior^{-1} mu
            //          - 2 + log(det_prior / det_post)]
            double S11 = vp.s2_delta[k];
            double S22 = vp.s2_psi[k];
            double S12 = vp.cov_delta_psi[k];
            double det_post = S11 * S22 - S12 * S12;
            if (det_post < 1e-20) det_post = 1e-20;

            // Sigma_prior^{-1} entries
            double P11_pr = t3 / det_prior;
            double P22_pr = t1 / det_prior;
            double P12_pr = -cov_prior / det_prior;

            // tr(Sigma_prior^{-1} * Sigma_post)
            double trace = P11_pr * S11 + P22_pr * S22 + 2.0 * P12_pr * S12;

            // mu' Sigma_prior^{-1} mu
            double md = vp.mu_delta[k];
            double mp = vp.mu_psi[k];
            double quad = P11_pr * md * md + P22_pr * mp * mp + 2.0 * P12_pr * md * mp;

            double kl_gauss = 0.5 * (trace + quad - 2.0 + std::log(det_prior / det_post + TINY));

            elbo += -(kl_bern + omega * kl_gauss);
        }
    } else {
        // Independent KL for delta and psi when those terms are active.
        for (int k = 0; k < nAC; k++) {
            if (delta_free) {
                elbo += kl_spike_slab(vp.omega_delta[k], vp.mu_delta[k], vp.s2_delta[k],
                                      hyp.pi1, hyp.tau2_1);
            }
            if (psi_free) {
                elbo += kl_spike_slab(vp.omega_psi[k], vp.mu_psi[k], vp.s2_psi[k],
                                      hyp.pi3, hyp.tau2_3);
            }
        }
    }
    // phi: only when protein->disease is active
    for (int l = 0; l < prot.nB(); l++) {
        if (phi_free) {
            elbo += kl_spike_slab(vp.omega_phi[l], vp.mu_phi[l], vp.s2_phi[l],
                                  hyp.pi2_cis, hyp.tau2_2_cis);
        }
    }

    if (scenario_correlated_pleiotropy(scenario)) {
        elbo += compute_m3_correlation_bonus(prot, vp);
    }

    return elbo;
}

// ============================================================================
// Run CAVI for a single protein under a single scenario
// Returns the final ELBO
// ============================================================================
double run_cavi(const ProteinData& prot, int scenario,
                VarParams& vp, const Hyperparams& hyp,
                const Options& opts) {
    int nAC = prot.nA() + prot.nC();
    double prev_elbo = -1e30;

    for (int iter = 0; iter < opts.max_cavi_iter; iter++) {
        // Update causal effects, respecting scenario constraints:
        // M0 null; M1 partial mediation; M2 RF->PP only; M3 RF->disease only;
        // M4 PP->disease only; M5 correlated/shared pleiotropy.
        if (scenario_beta1_free(scenario)) {
            update_beta1(prot, vp, hyp);
        }
        if (scenario_beta2_free(scenario)) {
            update_beta2(prot, vp, hyp);
        }
        if (scenario_beta3_free(scenario)) {
            update_beta3(prot, vp, hyp);
        }

        // Update spike-and-slab terms
        for (int k = 0; k < nAC; k++) {
            if (scenario_correlated_pleiotropy(scenario)) {
                // Correlated/shared pleiotropy: joint bivariate update for (delta_k, psi_k)
                // with correlated pleiotropy prior
                update_bivariate_delta_psi(prot, k, vp, hyp);
            } else {
                if (scenario_delta_free(scenario)) {
                    update_spike_slab_delta(prot, k, vp, hyp);
                }
                if (scenario_psi_free(scenario)) {
                    update_spike_slab_psi(prot, k, vp, hyp);
                }
            }
        }
        // phi is a nuisance cis-to-outcome path and is available in every state.
        for (int l = 0; l < prot.nB(); l++) {
            if (scenario_phi_free(scenario)) {
                update_spike_slab_phi(prot, l, vp, hyp);
            }
        }

        // Compute ELBO
        double elbo = compute_elbo(prot, scenario, vp, hyp);
        if (!std::isfinite(elbo)) {
            vp.elbo = -1e30;
            return vp.elbo;
        }
        vp.elbo = elbo;

        // Check convergence
        if (std::fabs(elbo - prev_elbo) < opts.elbo_tol) {
            break;
        }
        prev_elbo = elbo;
    }

    return vp.elbo;
}

// ============================================================================
// Analyze a single protein: run CAVI under all scenarios, compute posteriors
// ============================================================================
ProteinResult analyze_protein(const ProteinData& prot,
                              const Hyperparams& hyp,
                              const Options& opts) {
    ProteinResult res;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    res.protein_id = prot.protein_id;
    res.gene_name = prot.gene_name;
    res.nA = prot.nA();
    res.nB = prot.nB();
    res.nC = prot.nC();
    res.nC_exact = prot.nC_exact;
    res.nC_proxy = prot.nC_proxy;
    res.converged = true;
    res.ivw_rf_to_pp_beta = nan;
    res.ivw_rf_to_pp_se = nan;
    res.ivw_rf_to_pp_p = nan;
    res.ivw_pp_to_outcome_beta = nan;
    res.ivw_pp_to_outcome_se = nan;
    res.ivw_pp_to_outcome_p = nan;
    res.ivw_rf_to_outcome_beta = nan;
    res.ivw_rf_to_outcome_se = nan;
    res.ivw_rf_to_outcome_p = nan;
    res.indirect_direction = "NA";
    res.rf_to_outcome_direction = "NA";
    res.direction_consistent = "NA";
    res.direction_consistency_prob = 0.5;
    res.proportion_mediated = nan;
    res.directional_mediator_prob = nan;
    res.selection_probability = nan;
    res.selection_local_fdr = nan;
    res.selection_cum_fdr = nan;
    res.selection_rank = -1;
    res.posterior_local_fdr = nan;
    res.target_local_fdr = nan;
    res.posterior_cum_fdr = nan;
    res.posterior_cum_fdr5 = nan;
    res.mediation_rank = -1;
    res.selected_fdr_10 = false;
    res.selected_fdr_5 = false;
    res.evidence_tier = "none";
    res.eb_beta1_second_moment = 0.0;
    res.eb_beta2_second_moment = 0.0;
    res.eb_beta3_second_moment = 0.0;
    res.eb_delta_pi = 0.0;
    res.eb_delta_second_moment = 0.0;
    res.eb_phi_pi = 0.0;
    res.eb_phi_second_moment = 0.0;
    res.eb_psi_pi = 0.0;
    res.eb_psi_second_moment = 0.0;
    res.eb_m3_resid_corr = 0.0;
    int n_rf_to_pp = observed_rf_to_pp_count(prot);
    res.n_rf_to_pp_obs = n_rf_to_pp;
    res.rf_to_pp_identifiable = n_rf_to_pp > 0;
    RegionalEvidence regional = compute_regional_evidence(prot, opts);
    res.regional_n_variants = regional.n_variants;
    res.regional_pp_shared = regional.pp_shared;
    res.regional_pp_distinct = regional.pp_distinct;
    res.regional_shared_given_both = regional.shared_given_both;
    res.regional_method = prot.regional_multisignal_evaluated
        ? "ld-multisignal" : "single";
    res.regional_protein_signals = prot.regional_protein_signals;
    res.regional_outcome_signals = prot.regional_outcome_signals;
    res.regional_signal_pair_count = static_cast<int>(prot.regional_signal_pairs.size());
    res.regional_max_credible_set_pair_r2 = prot.regional_best_cs_pair_r2;
    res.regional_signal_pairs = prot.regional_signal_pairs;
    if (!prot.ld_reference_used) {
        res.mediation_identifiability = "UNRESOLVED_NO_LD_REFERENCE";
    } else if (prot.regional_multisignal_evaluated) {
        const std::string& interpretation = prot.regional_multisignal_interpretation;
        if (interpretation == "SHARED_SIGNAL_SUPPORTED") {
            if (n_rf_to_pp >= 2 && prot.nB() >= 1) {
                res.mediation_identifiability =
                    "LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL";
            } else if (n_rf_to_pp < 2 && prot.nB() < 1) {
                res.mediation_identifiability =
                    "UNRESOLVED_INSUFFICIENT_RF_AND_CIS_INSTRUMENTS";
            } else if (n_rf_to_pp < 2) {
                res.mediation_identifiability =
                    "UNRESOLVED_INSUFFICIENT_RF_INSTRUMENTS";
            } else {
                res.mediation_identifiability =
                    "UNRESOLVED_INSUFFICIENT_CIS_INSTRUMENTS";
            }
        } else if (interpretation == "DISTINCT_SIGNALS_HIGH_LD" ||
                   interpretation == "DISTINCT_SIGNALS_LOW_MODERATE_LD") {
            res.mediation_identifiability = "LD_DISTINCT_SUPPORTED";
        } else if (interpretation == "NO_OUTCOME_CREDIBLE_SET") {
            res.mediation_identifiability = "UNRESOLVED_NO_OUTCOME_SIGNAL";
        } else if (interpretation == "NO_PROTEIN_CREDIBLE_SET") {
            res.mediation_identifiability = "UNRESOLVED_NO_PROTEIN_SIGNAL";
        } else if (interpretation == "NO_REGIONAL_DATA") {
            res.mediation_identifiability = "UNRESOLVED_NO_REGIONAL_DATA";
        } else {
            res.mediation_identifiability = "LD_CONFIGURATION_AMBIGUOUS";
        }
    } else if (!prot.regional_data_complete || regional.n_variants < 2) {
        res.mediation_identifiability = "UNRESOLVED_NO_REGIONAL_DATA";
    } else if (regional.pp_shared + regional.pp_distinct < opts.regional_min_both) {
        res.mediation_identifiability = "UNRESOLVED_WEAK_REGIONAL_EVIDENCE";
    } else if (regional.shared_given_both >= opts.regional_min_shared) {
        if (n_rf_to_pp >= 2 && prot.nB() >= 1) {
            res.mediation_identifiability =
                "LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL";
        } else if (n_rf_to_pp < 2 && prot.nB() < 1) {
            res.mediation_identifiability =
                "UNRESOLVED_INSUFFICIENT_RF_AND_CIS_INSTRUMENTS";
        } else if (n_rf_to_pp < 2) {
            res.mediation_identifiability =
                "UNRESOLVED_INSUFFICIENT_RF_INSTRUMENTS";
        } else {
            res.mediation_identifiability =
                "UNRESOLVED_INSUFFICIENT_CIS_INSTRUMENTS";
        }
    } else if (regional.shared_given_both <= 1.0 - opts.regional_min_shared) {
        res.mediation_identifiability = "LD_DISTINCT_SUPPORTED";
    } else {
        res.mediation_identifiability = "LD_CONFIGURATION_AMBIGUOUS";
    }
    bool has_first_stage = n_rf_to_pp > 0;
    bool has_second_stage = (prot.nB() + prot.nC()) > 0;
    bool has_rf_outcome = (prot.nA() + prot.nC()) > 0;

    {
        std::vector<double> rf_to_pp_bx;
        std::vector<double> rf_to_pp_by;
        std::vector<double> rf_to_pp_se;
        rf_to_pp_bx.reserve(prot.nA() + prot.nC());
        rf_to_pp_by.reserve(prot.nA() + prot.nC());
        rf_to_pp_se.reserve(prot.nA() + prot.nC());
        for (int k = 0; k < prot.nA(); ++k) {
            if (!has_setA_alpha(prot, k)) continue;
            rf_to_pp_bx.push_back(prot.setA_gamma[k]);
            rf_to_pp_by.push_back(prot.setA_alpha[k]);
            rf_to_pp_se.push_back(prot.setA_se_alpha[k]);
        }
        for (int c = 0; c < prot.nC(); ++c) {
            rf_to_pp_bx.push_back(prot.setC_gamma[c]);
            rf_to_pp_by.push_back(prot.setC_alpha[c]);
            rf_to_pp_se.push_back(prot.setC_se_alpha[c]);
        }
        IvwSummary ivw = compute_ivw_summary(rf_to_pp_bx, rf_to_pp_by, rf_to_pp_se);
        res.ivw_rf_to_pp_beta = ivw.beta;
        res.ivw_rf_to_pp_se = ivw.se;
        res.ivw_rf_to_pp_p = ivw.p;
    }

    {
        std::vector<double> pp_to_outcome_bx;
        std::vector<double> pp_to_outcome_by;
        std::vector<double> pp_to_outcome_se;
        pp_to_outcome_bx.reserve(prot.nB() + prot.nC());
        pp_to_outcome_by.reserve(prot.nB() + prot.nC());
        pp_to_outcome_se.reserve(prot.nB() + prot.nC());
        for (int l = 0; l < prot.nB(); ++l) {
            pp_to_outcome_bx.push_back(prot.setB_alpha_cis[l]);
            pp_to_outcome_by.push_back(prot.setB_Gamma_cis[l]);
            pp_to_outcome_se.push_back(prot.setB_se_Gamma_cis[l]);
        }
        for (int c = 0; c < prot.nC(); ++c) {
            pp_to_outcome_bx.push_back(prot.setC_alpha[c]);
            pp_to_outcome_by.push_back(prot.setC_Gamma[c]);
            pp_to_outcome_se.push_back(prot.setC_se_Gamma[c]);
        }
        IvwSummary ivw = compute_ivw_summary(pp_to_outcome_bx, pp_to_outcome_by, pp_to_outcome_se);
        res.ivw_pp_to_outcome_beta = ivw.beta;
        res.ivw_pp_to_outcome_se = ivw.se;
        res.ivw_pp_to_outcome_p = ivw.p;
    }

    {
        std::vector<double> rf_to_outcome_bx;
        std::vector<double> rf_to_outcome_by;
        std::vector<double> rf_to_outcome_se;
        rf_to_outcome_bx.reserve(prot.nA() + prot.nC());
        rf_to_outcome_by.reserve(prot.nA() + prot.nC());
        rf_to_outcome_se.reserve(prot.nA() + prot.nC());
        for (int k = 0; k < prot.nA(); ++k) {
            rf_to_outcome_bx.push_back(prot.setA_gamma[k]);
            rf_to_outcome_by.push_back(prot.setA_Gamma[k]);
            rf_to_outcome_se.push_back(prot.setA_se_Gamma[k]);
        }
        for (int c = 0; c < prot.nC(); ++c) {
            rf_to_outcome_bx.push_back(prot.setC_gamma[c]);
            rf_to_outcome_by.push_back(prot.setC_Gamma[c]);
            rf_to_outcome_se.push_back(prot.setC_se_Gamma[c]);
        }
        IvwSummary ivw = compute_ivw_summary(rf_to_outcome_bx, rf_to_outcome_by, rf_to_outcome_se);
        res.ivw_rf_to_outcome_beta = ivw.beta;
        res.ivw_rf_to_outcome_se = ivw.se;
        res.ivw_rf_to_outcome_p = ivw.p;
    }

    // Skip proteins with no instruments at all
    if (prot.nTotal() == 0) {
        res.prob_M0 = 1.0; res.prob_M1 = 0.0;
        res.prob_M2 = 0.0; res.prob_M3 = 0.0;
        res.prob_M4 = 0.0; res.prob_M5 = 0.0;
        res.prob_mediator = 0.0;
        res.prob_mediator_ld_resolved = 0.0;
        res.prob_mediator_identified = 0.0;
        res.prob_protein_disease = 0.0;
        res.prob_rf_responsive = 0.0;
        res.prob_rf_direct = 0.0;
        res.beta1_est = 0; res.beta1_se = 0;
        res.beta2_est = 0; res.beta2_se = 0;
        res.beta3_est = 0; res.beta3_se = 0;
        res.mediated_effect = 0; res.mediated_effect_se = 0;
        populate_direction_descriptives(res);
        res.directional_mediator_prob = 0.0;
        res.selection_probability = 0.0;
        res.selection_local_fdr = 1.0;
        res.selection_cum_fdr = nan;
        res.selection_rank = -1;
        res.elbo_M0 = 0; res.elbo_M1 = 0;
        res.elbo_M2 = 0; res.elbo_M3 = 0;
        res.elbo_M4 = 0; res.elbo_M5 = 0;
        res.converged = false;
        return res;
    }

    // Run CAVI under each scenario
    constexpr int N_SCENARIOS = 6;
    double elbos[N_SCENARIOS];
    double m1_resid_corr = 0.0;
    LocalScales local_scales = build_local_scales(hyp,
                                                  {res.ivw_rf_to_pp_beta, res.ivw_rf_to_pp_se, res.ivw_rf_to_pp_p},
                                                  {res.ivw_pp_to_outcome_beta, res.ivw_pp_to_outcome_se, res.ivw_pp_to_outcome_p},
                                                  opts.legacy_adaptive_priors);

    for (int m = 0; m < N_SCENARIOS; m++) {
        VarParams vp = init_var_params(prot, m, hyp);
        vp.prior_sigma2_beta1 = local_scales.sigma2_beta1;
        vp.prior_sigma2_beta2 = local_scales.sigma2_beta2;
        vp.prior_sigma2_beta3 = local_scales.sigma2_beta3;
        vp.s2_beta1 = scenario_beta1_free(m) ? vp.prior_sigma2_beta1 : 0.0;
        vp.s2_beta2 = scenario_beta2_free(m) ? vp.prior_sigma2_beta2 : 0.0;
        vp.s2_beta3 = scenario_beta3_free(m) ? vp.prior_sigma2_beta3 : 0.0;
        elbos[m] = run_cavi(prot, m, vp, hyp, opts);

        // Store M=1 estimates
        if (m == 1) {
            res.beta1_est = vp.mu_beta1;
            res.beta1_se  = std::sqrt(vp.s2_beta1);
            res.beta2_est = vp.mu_beta2;
            res.beta2_se  = std::sqrt(vp.s2_beta2);
            res.beta3_est = vp.mu_beta3;
            res.beta3_se  = std::sqrt(vp.s2_beta3);

            // Mediated effect: beta1 * beta2
            res.mediated_effect = vp.mu_beta1 * vp.mu_beta2;
            // Delta method SE: sqrt(beta2^2*var_beta1 + beta1^2*var_beta2)
            res.mediated_effect_se = std::sqrt(
                vp.mu_beta2 * vp.mu_beta2 * vp.s2_beta1 +
                vp.mu_beta1 * vp.mu_beta1 * vp.s2_beta2
            );
            res.eb_beta1_second_moment = vp.mu_beta1 * vp.mu_beta1 + vp.s2_beta1;
            res.eb_beta2_second_moment = vp.mu_beta2 * vp.mu_beta2 + vp.s2_beta2;
            res.eb_beta3_second_moment = vp.mu_beta3 * vp.mu_beta3 + vp.s2_beta3;
            summarize_spike_slab(vp.omega_delta, vp.mu_delta, vp.s2_delta,
                                 res.eb_delta_pi, res.eb_delta_second_moment);
            summarize_spike_slab(vp.omega_phi, vp.mu_phi, vp.s2_phi,
                                 res.eb_phi_pi, res.eb_phi_second_moment);
            summarize_spike_slab(vp.omega_psi, vp.mu_psi, vp.s2_psi,
                                 res.eb_psi_pi, res.eb_psi_second_moment);
            m1_resid_corr = std::fabs(compute_m3_residual_corr(prot, vp));
        }
        if (m == 3) {
            res.eb_m3_resid_corr = std::fabs(compute_m3_residual_corr(prot, vp));
        }
    }

    res.elbo_M0 = elbos[0];
    res.elbo_M1 = elbos[1];
    res.elbo_M2 = elbos[2];
    res.elbo_M3 = elbos[3];
    res.elbo_M4 = elbos[4];
    res.elbo_M5 = elbos[5];
    populate_direction_descriptives(res);

    // Compute posterior scenario probabilities
    SoftPriors soft_priors = build_soft_priors(hyp,
                                               {res.ivw_rf_to_pp_beta, res.ivw_rf_to_pp_se, res.ivw_rf_to_pp_p},
                                               {res.ivw_pp_to_outcome_beta, res.ivw_pp_to_outcome_se, res.ivw_pp_to_outcome_p},
                                               {res.ivw_rf_to_outcome_beta, res.ivw_rf_to_outcome_se, res.ivw_rf_to_outcome_p},
                                               opts.legacy_adaptive_priors);
    double log_priors[N_SCENARIOS] = {
        std::log(soft_priors.p0 + TINY),
        std::log(soft_priors.p1 + TINY),
        std::log(soft_priors.p2 + TINY),
        std::log(soft_priors.p3 + TINY),
        std::log(soft_priors.p4 + TINY),
        std::log(soft_priors.p5 + TINY)
    };

    double log_posts[N_SCENARIOS];
    for (int m = 0; m < N_SCENARIOS; m++) {
        log_posts[m] = log_priors[m] + elbos[m];
        if (!std::isfinite(log_posts[m])) {
            log_posts[m] = -1e30;
            res.converged = false;
        }
    }
    if (!has_second_stage) {
        log_posts[1] = -1e30;  // cannot support mediation without PP->outcome instruments
        log_posts[4] = -1e30;  // cannot support PP->disease only
    }
    if (!has_first_stage) {
        log_posts[1] = -1e30;  // cannot support mediation without RF->PP information
        log_posts[2] = -1e30;  // cannot support RF->PP only
        log_posts[5] = -1e30;  // cannot support shared RF/protein pleiotropy
    }
    if (!has_rf_outcome) {
        log_posts[3] = -1e30;  // cannot support RF direct without RF outcome instruments
        log_posts[5] = -1e30;  // cannot support correlated pleiotropy without RF outcome instruments
    }
    if (opts.m1_min_cis_only > 0 && prot.nB() < opts.m1_min_cis_only) {
        log_posts[1] = -1e30;
    }
    if (opts.m1_min_first_stage_z > 0.0) {
        double z_first = (std::isfinite(res.ivw_rf_to_pp_beta) &&
                          std::isfinite(res.ivw_rf_to_pp_se) &&
                          res.ivw_rf_to_pp_se > 0.0)
                            ? std::fabs(res.ivw_rf_to_pp_beta / res.ivw_rf_to_pp_se)
                            : 0.0;
        if (z_first < opts.m1_min_first_stage_z) {
            log_posts[1] = -1e30;
        }
    }
    if (opts.m1_min_second_stage_z > 0.0) {
        double z_second = (std::isfinite(res.ivw_pp_to_outcome_beta) &&
                           std::isfinite(res.ivw_pp_to_outcome_se) &&
                           res.ivw_pp_to_outcome_se > 0.0)
                            ? std::fabs(res.ivw_pp_to_outcome_beta / res.ivw_pp_to_outcome_se)
                            : 0.0;
        if (z_second < opts.m1_min_second_stage_z) {
            log_posts[1] = -1e30;
        }
    }
    if (opts.m1_resid_corr_penalty > 0.0 &&
        m1_resid_corr > opts.m1_resid_corr_threshold &&
        log_posts[1] > -1e20) {
        double excess = m1_resid_corr - opts.m1_resid_corr_threshold;
        double n_eff = std::max(1.0, static_cast<double>(prot.nA() + prot.nC()));
        log_posts[1] -= opts.m1_resid_corr_penalty * n_eff * excess * excess;
    }
    if (opts.direction_mode == "soft") {
        double pdir = std::min(std::max(res.direction_consistency_prob, 1e-8), 1.0);
        log_posts[1] += opts.direction_weight * std::log(pdir);
    } else if (opts.direction_mode == "hard") {
        if (res.direction_consistency_prob < opts.direction_min_prob) {
            log_posts[1] = -1e30;
        }
    }

    // Softmax
    std::vector<double> lp_vec(log_posts, log_posts + N_SCENARIOS);
    double lse = log_sum_exp(lp_vec);
    res.prob_M0 = std::exp(log_posts[0] - lse);
    res.prob_M1 = std::exp(log_posts[1] - lse);
    res.prob_M2 = std::exp(log_posts[2] - lse);
    res.prob_M3 = std::exp(log_posts[3] - lse);
    res.prob_M4 = std::exp(log_posts[4] - lse);
    res.prob_M5 = std::exp(log_posts[5] - lse);
    res.prob_mediator = res.prob_M1;
    res.prob_protein_disease = res.prob_M1 + res.prob_M4;
    res.prob_rf_responsive = res.prob_M1 + res.prob_M2 + res.prob_M5;
    res.prob_rf_direct = res.prob_M1 + res.prob_M3 + res.prob_M5;
    res.posterior_local_fdr = 1.0 - res.prob_M1;
    res.target_local_fdr = 1.0 - res.prob_protein_disease;
    finalize_direction_metrics(res);
    bool ld_resolved = res.mediation_identifiability ==
                       "LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL";
    res.prob_mediator_ld_resolved = ld_resolved ? res.prob_M1 : 0.0;
    res.prob_mediator_identified = res.prob_mediator_ld_resolved;
    res.selection_probability = selection_probability_for_mode(res, opts);
    res.selection_local_fdr = 1.0 - res.selection_probability;

    bool strong_rf = std::isfinite(res.ivw_rf_to_pp_p) && res.ivw_rf_to_pp_p < 0.05;
    bool strong_out = std::isfinite(res.ivw_pp_to_outcome_p) && res.ivw_pp_to_outcome_p < 0.05;
    if (res.prob_mediator_ld_resolved >= 0.5 && strong_rf && strong_out) {
        res.evidence_tier = "high";
    } else if (res.prob_M1 >= 0.5 && !ld_resolved) {
        res.evidence_tier = "unresolved";
    } else if (res.prob_M1 >= 0.1 || (strong_rf && strong_out)) {
        res.evidence_tier = "moderate";
    } else if (res.prob_M1 >= 0.01 || strong_rf || strong_out) {
        res.evidence_tier = "exploratory";
    } else {
        res.evidence_tier = "none";
    }

    if (opts.verbose) {
        std::cout << "    " << prot.protein_id
                  << " -> P(M1)=" << std::fixed << std::setprecision(4) << res.prob_M1
                  << ", P(M4)=" << res.prob_M4
                  << " [nA=" << prot.nA() << ", nB=" << prot.nB() << ", nC=" << prot.nC() << "]\n";
    }

    return res;
}

// ============================================================================
// Empirical Bayes: update hyperparameters across all proteins
// ============================================================================
void run_empirical_bayes(std::vector<ProteinData>& proteins,
                         Hyperparams& hyp,
                         std::vector<ProteinResult>& results,
                         const Options& opts) {
    int J = (int)proteins.size();
    hyp.p0 = opts.prior_p0;
    hyp.p1 = opts.prior_p1;
    hyp.p2 = opts.prior_p2;
    hyp.p3 = opts.prior_p3;
    hyp.p4 = opts.prior_p4;
    hyp.p5 = opts.prior_p5;
    const double init_p0 = hyp.p0, init_p1 = hyp.p1, init_p2 = hyp.p2;
    const double init_p3 = hyp.p3, init_p4 = hyp.p4, init_p5 = hyp.p5;
    hyp.sigma2_beta1 = opts.prior_sigma2_beta1;
    hyp.sigma2_beta2 = opts.prior_sigma2_beta2;
    hyp.sigma2_beta3 = opts.prior_sigma2_beta3;
    const double sigma_floor = opts.prior_sigma2_floor;
    const double sigma1_anchor = std::max(opts.prior_sigma2_beta1, sigma_floor);
    const double sigma2_anchor = std::max(opts.prior_sigma2_beta2, sigma_floor);
    const double sigma3_anchor = std::max(opts.prior_sigma2_beta3, sigma_floor);

    for (int eb_iter = 0; eb_iter < opts.max_eb_iter; eb_iter++) {
        std::cout << "EB iteration " << (eb_iter + 1) << " / " << opts.max_eb_iter << " ...";
        std::cout.flush();

        double old_p0 = hyp.p0, old_p1 = hyp.p1;
        double old_p2 = hyp.p2, old_p3 = hyp.p3;
        double old_p4 = hyp.p4, old_p5 = hyp.p5;

        results.resize(J);

        // E-step: run CAVI for each protein
        // (could be parallelized with OpenMP here)
#ifdef _OPENMP
        #pragma omp parallel for num_threads(opts.threads) schedule(dynamic)
#endif
        for (int j = 0; j < J; j++) {
            results[j] = analyze_protein(proteins[j], hyp, opts);
            if (opts.verbose) {
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    if ((j + 1) % 100 == 0 || j + 1 == J) {
                        std::cout << "      analyzed " << (j + 1) << " / " << J
                                  << " proteins in EB iteration " << (eb_iter + 1) << "\n";
                    }
                }
            }
        }

        if (opts.fixed_priors) {
            std::cout << " fixed hyperparameters\n";
            break;
        }

        // M-step: update prior mixture weights
        double sum_p0 = 0, sum_p1 = 0, sum_p2 = 0, sum_p3 = 0, sum_p4 = 0, sum_p5 = 0;
        for (int j = 0; j < J; j++) {
            sum_p0 += results[j].prob_M0;
            sum_p1 += results[j].prob_M1;
            sum_p2 += results[j].prob_M2;
            sum_p3 += results[j].prob_M3;
            sum_p4 += results[j].prob_M4;
            sum_p5 += results[j].prob_M5;
        }
        double den = static_cast<double>(J) + opts.eb_prior_strength;
        hyp.p0 = (sum_p0 + opts.eb_prior_strength * init_p0) / den;
        hyp.p1 = (sum_p1 + opts.eb_prior_strength * init_p1) / den;
        hyp.p2 = (sum_p2 + opts.eb_prior_strength * init_p2) / den;
        hyp.p3 = (sum_p3 + opts.eb_prior_strength * init_p3) / den;
        hyp.p4 = (sum_p4 + opts.eb_prior_strength * init_p4) / den;
        hyp.p5 = (sum_p5 + opts.eb_prior_strength * init_p5) / den;

        // Ensure minimum probability
        double pmin = 1e-4;
        hyp.p0 = std::max(hyp.p0, pmin);
        hyp.p1 = std::max(hyp.p1, pmin);
        hyp.p2 = std::max(hyp.p2, pmin);
        hyp.p3 = std::max(hyp.p3, pmin);
        hyp.p4 = std::max(hyp.p4, pmin);
        hyp.p5 = std::max(hyp.p5, pmin);
        double psum = hyp.p0 + hyp.p1 + hyp.p2 + hyp.p3 + hyp.p4 + hyp.p5;
        hyp.p0 /= psum; hyp.p1 /= psum; hyp.p2 /= psum; hyp.p3 /= psum;
        hyp.p4 /= psum; hyp.p5 /= psum;

        // Update continuous effect and pleiotropy scales via EB moment matching.
        auto safe_weighted_mean = [](double num, double den, double fallback) {
            return (den > 0.0 && std::isfinite(num)) ? (num / den) : fallback;
        };
        double beta1_num = 0.0, beta1_den = 0.0;
        double beta2_num = 0.0, beta2_den = 0.0;
        double beta3_num = 0.0, beta3_den = 0.0;
        double delta_pi_num = 0.0, delta_den = 0.0, delta_m2_num = 0.0;
        double phi_pi_num = 0.0, phi_den = 0.0, phi_m2_num = 0.0;
        double psi_pi_num = 0.0, psi_den = 0.0, psi_m2_num = 0.0;
        double rho_num = 0.0, rho_den = 0.0;
        double beta1_cap = std::max(sigma1_anchor * 25.0, 1.0);
        double beta2_cap = std::max(sigma2_anchor * 10.0, 0.25);
        double beta3_cap = std::max(sigma3_anchor * 10.0, 0.25);
        for (const auto& r : results) {
            double info1 = std::sqrt(std::max(r.n_rf_to_pp_obs, 0));
            double w_beta1 = (r.prob_M1 + r.prob_M2 + r.prob_M5) * std::max(1.0, info1);
            double w_beta2 = (r.prob_M1 + r.prob_M4) * std::max(1.0, info1);
            double w_beta3 = (r.prob_M1 + r.prob_M3 + r.prob_M5) * std::max(1.0, info1);
            beta1_num += w_beta1 * std::min(r.eb_beta1_second_moment, beta1_cap);
            beta1_den += w_beta1;
            beta2_num += w_beta2 * std::min(r.eb_beta2_second_moment, beta2_cap);
            beta2_den += w_beta2;
            beta3_num += w_beta3 * std::min(r.eb_beta3_second_moment, beta3_cap);
            beta3_den += w_beta3;

            double w_pleio = std::max((r.prob_M1 + r.prob_M2 + r.prob_M3 + r.prob_M4 + r.prob_M5) *
                                      std::max(1.0, info1), 1e-8);
            delta_pi_num += w_pleio * r.eb_delta_pi;
            delta_m2_num += w_pleio * r.eb_delta_second_moment;
            delta_den += w_pleio;
            phi_pi_num += w_pleio * r.eb_phi_pi;
            phi_m2_num += w_pleio * r.eb_phi_second_moment;
            phi_den += w_pleio;
            psi_pi_num += w_pleio * r.eb_psi_pi;
            psi_m2_num += w_pleio * r.eb_psi_second_moment;
            psi_den += w_pleio;

            double w_m3 = std::max(r.prob_M5 * std::max(1.0, info1), 0.0);
            rho_num += w_m3 * r.eb_m3_resid_corr;
            rho_den += w_m3;
        }
        double beta1_est = safe_weighted_mean(beta1_num, beta1_den, hyp.sigma2_beta1);
        double beta2_est = safe_weighted_mean(beta2_num, beta2_den, hyp.sigma2_beta2);
        double beta3_est = safe_weighted_mean(beta3_num, beta3_den, hyp.sigma2_beta3);
        hyp.sigma2_beta1 = stabilized_scale_update(hyp.sigma2_beta1, beta1_est,
                                                   sigma1_anchor, sigma_floor, 25.0);
        hyp.sigma2_beta2 = stabilized_scale_update(hyp.sigma2_beta2, beta2_est,
                                                   sigma2_anchor, sigma_floor, 10.0);
        hyp.sigma2_beta3 = stabilized_scale_update(hyp.sigma2_beta3, beta3_est,
                                                   sigma3_anchor, sigma_floor, 10.0);
        hyp.pi1 = bounded_prob(safe_weighted_mean(delta_pi_num, delta_den, hyp.pi1), 1e-3, 0.5);
        hyp.tau2_1 = stabilized_scale_update(hyp.tau2_1,
                                             safe_weighted_mean(delta_m2_num, delta_den, hyp.tau2_1),
                                             0.01, 1e-5, 5.0, 0.20);
        hyp.pi2_cis = bounded_prob(safe_weighted_mean(phi_pi_num, phi_den, hyp.pi2_cis), 1e-3, 0.5);
        hyp.tau2_2_cis = stabilized_scale_update(hyp.tau2_2_cis,
                                                 safe_weighted_mean(phi_m2_num, phi_den, hyp.tau2_2_cis),
                                                 0.01, 1e-5, 5.0, 0.20);
        hyp.pi3 = bounded_prob(safe_weighted_mean(psi_pi_num, psi_den, hyp.pi3), 1e-3, 0.5);
        hyp.tau2_3 = stabilized_scale_update(hyp.tau2_3,
                                             safe_weighted_mean(psi_m2_num, psi_den, hyp.tau2_3),
                                             0.01, 1e-5, 5.0, 0.20);
        {
            double rho_est = safe_weighted_mean(rho_num, rho_den, hyp.rho_prior);
            rho_est = clamp_value(rho_est, 0.05, 0.995);
            hyp.rho_prior = clamp_value(0.80 * hyp.rho_prior + 0.20 * rho_est, 0.05, 0.995);
        }

        std::cout << " priors = (" << std::fixed << std::setprecision(4)
                  << hyp.p0 << ", " << hyp.p1 << ", "
                  << hyp.p2 << ", " << hyp.p3 << ", "
                  << hyp.p4 << ", " << hyp.p5 << ")\n";
        if (opts.verbose) {
            int strong_m1 = 0;
            for (const auto& r : results) {
                if (r.prob_M1 > 0.5) strong_m1++;
            }
            std::cout << "      proteins with P(M1) > 0.5 this iteration: "
                      << strong_m1 << "\n"
                      << "      sigma2 = (" << hyp.sigma2_beta1 << ", "
                      << hyp.sigma2_beta2 << ", " << hyp.sigma2_beta3 << ")\n"
                      << "      rho_prior = " << hyp.rho_prior << "\n";
        }

        // Check convergence
        double delta = std::fabs(hyp.p0 - old_p0) + std::fabs(hyp.p1 - old_p1)
                     + std::fabs(hyp.p2 - old_p2) + std::fabs(hyp.p3 - old_p3)
                     + std::fabs(hyp.p4 - old_p4) + std::fabs(hyp.p5 - old_p5);
        if (!opts.fixed_priors && delta < opts.eb_tol) {
            std::cout << "EB converged after " << (eb_iter + 1) << " iterations.\n";
            break;
        }
    }

    std::vector<int> order(results.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return results[a].prob_M1 > results[b].prob_M1;
    });
    double cum_err = 0.0;
    for (size_t rank = 0; rank < order.size(); ++rank) {
        auto& r = results[order[rank]];
        r.mediation_rank = (int)rank + 1;
        r.posterior_local_fdr = 1.0 - r.prob_M1;
        cum_err += r.posterior_local_fdr;
        r.posterior_cum_fdr = cum_err / (rank + 1.0);
        r.posterior_cum_fdr5 = r.posterior_cum_fdr;
    }

    std::vector<int> selection_order(results.size());
    std::iota(selection_order.begin(), selection_order.end(), 0);
    std::sort(selection_order.begin(), selection_order.end(), [&](int a, int b) {
        return results[a].selection_probability > results[b].selection_probability;
    });
    double selection_cum_err = 0.0;
    for (size_t rank = 0; rank < selection_order.size(); ++rank) {
        auto& r = results[selection_order[rank]];
        r.selection_rank = (int)rank + 1;
        r.selection_local_fdr = 1.0 - r.selection_probability;
        selection_cum_err += r.selection_local_fdr;
        r.selection_cum_fdr = selection_cum_err / (rank + 1.0);
        bool has_leg_estimates = std::isfinite(r.ivw_rf_to_pp_beta) &&
                                 std::isfinite(r.ivw_pp_to_outcome_beta);
        r.selected_fdr_10 = r.selection_cum_fdr <= 0.10 && has_leg_estimates;
        r.selected_fdr_5 = r.selection_cum_fdr <= 0.05 && has_leg_estimates;
    }
}

} // namespace bmediator
