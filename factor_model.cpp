#include "bmediator.h"

#include <functional>

namespace bmediator {
namespace {

using Matrix = std::vector<std::vector<double>>;

struct FactorFit {
    double beta = std::numeric_limits<double>::quiet_NaN();
    double se = std::numeric_limits<double>::quiet_NaN();
    double p = std::numeric_limits<double>::quiet_NaN();
    double log_bf = std::numeric_limits<double>::quiet_NaN();
    int n = 0;
    bool valid = false;
};

double clamp_corr(double value) {
    return std::max(-0.999, std::min(0.999, value));
}

double beta_continued_fraction(double a, double b, double x) {
    constexpr int max_iter = 200;
    constexpr double fp_min = 1e-300;
    constexpr double tolerance = 3e-14;
    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < fp_min) d = fp_min;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= max_iter; ++m) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < fp_min) d = fp_min;
        c = 1.0 + aa / c;
        if (std::fabs(c) < fp_min) c = fp_min;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < fp_min) d = fp_min;
        c = 1.0 + aa / c;
        if (std::fabs(c) < fp_min) c = fp_min;
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::fabs(delta - 1.0) < tolerance) break;
    }
    return h;
}

double regularized_incomplete_beta(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double front = std::exp(std::lgamma(a + b) - std::lgamma(a) -
                                  std::lgamma(b) + a * std::log(x) +
                                  b * std::log1p(-x));
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return front * beta_continued_fraction(a, b, x) / a;
    }
    return 1.0 - front * beta_continued_fraction(b, a, 1.0 - x) / b;
}

double two_sided_student_t_p(double statistic, double df) {
    if (!(df > 0.0) || !std::isfinite(statistic)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double x = df / (df + statistic * statistic);
    return regularized_incomplete_beta(0.5 * df, 0.5, x);
}

Matrix identity_matrix(int n) {
    Matrix result(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) result[i][i] = 1.0;
    return result;
}

bool valid_ld(const Matrix& ld, int n) {
    if (static_cast<int>(ld.size()) != n) return false;
    for (const auto& row : ld) {
        if (static_cast<int>(row.size()) != n) return false;
    }
    return true;
}

Matrix subset_matrix(const Matrix& matrix, const std::vector<int>& keep) {
    if (!valid_ld(matrix, static_cast<int>(matrix.size()))) {
        return identity_matrix(static_cast<int>(keep.size()));
    }
    Matrix result(keep.size(), std::vector<double>(keep.size(), 0.0));
    for (size_t i = 0; i < keep.size(); ++i) {
        for (size_t j = 0; j < keep.size(); ++j) {
            result[i][j] = matrix[keep[i]][keep[j]];
        }
    }
    return result;
}

// Cholesky factorization with a small diagonal ridge for reference-panel noise.
bool cholesky(const Matrix& input, Matrix& lower, double& log_det) {
    const int n = static_cast<int>(input.size());
    for (int attempt = 0; attempt < 8; ++attempt) {
        const double ridge = attempt == 0 ? 0.0 : std::pow(10.0, attempt - 11);
        lower.assign(n, std::vector<double>(n, 0.0));
        bool ok = true;
        for (int i = 0; i < n && ok; ++i) {
            for (int j = 0; j <= i; ++j) {
                double value = input[i][j];
                if (i == j) value += ridge;
                for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
                if (i == j) {
                    if (!(value > 0.0) || !std::isfinite(value)) {
                        ok = false;
                        break;
                    }
                    lower[i][j] = std::sqrt(value);
                } else {
                    lower[i][j] = value / lower[j][j];
                }
            }
        }
        if (ok) {
            log_det = 0.0;
            for (int i = 0; i < n; ++i) log_det += 2.0 * std::log(lower[i][i]);
            return true;
        }
    }
    return false;
}

double quadratic_form_from_cholesky(const Matrix& lower,
                                    const std::vector<double>& residual) {
    const int n = static_cast<int>(residual.size());
    std::vector<double> z(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double value = residual[i];
        for (int j = 0; j < i; ++j) value -= lower[i][j] * z[j];
        z[i] = value / lower[i][i];
    }
    return std::inner_product(z.begin(), z.end(), z.begin(), 0.0);
}

std::vector<double> whiten(const Matrix& lower, const std::vector<double>& value) {
    const int n = static_cast<int>(value.size());
    std::vector<double> result(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double current = value[i];
        for (int j = 0; j < i; ++j) current -= lower[i][j] * result[j];
        result[i] = current / lower[i][i];
    }
    return result;
}

// Null score test conditional on the selected exposure associations. With
// independent exposure and outcome GWAS errors this is Gaussian under H0,
// even when instruments are weak. Multiplicative overdispersion is only
// allowed to increase the variance.
double gls_score_p(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& se_y,
                   const Matrix& ld) {
    const int n = static_cast<int>(x.size());
    Matrix covariance(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            covariance[i][j] = ld[i][j] * se_y[i] * se_y[j];
        }
    }
    Matrix lower;
    double log_det = 0.0;
    if (!cholesky(covariance, lower, log_det)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const auto wx = whiten(lower, x);
    const auto wy = whiten(lower, y);
    const double q = std::inner_product(wx.begin(), wx.end(), wx.begin(), 0.0);
    const double t = std::inner_product(wx.begin(), wx.end(), wy.begin(), 0.0);
    if (!(q > 0.0) || !std::isfinite(q) || !std::isfinite(t)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double beta = t / q;
    double residual_q = 0.0;
    for (int i = 0; i < n; ++i) {
        const double residual = wy[i] - beta * wx[i];
        residual_q += residual * residual;
    }
    const double estimated_dispersion = n > 1 ? residual_q / (n - 1.0) : 1.0;
    const double dispersion = std::max(1.0, estimated_dispersion);
    const double z = t / std::sqrt(q * dispersion);
    return n > 1
        ? two_sided_student_t_p(z, n - 1.0)
        : std::erfc(std::fabs(z) / std::sqrt(2.0));
}

double eiv_log_likelihood(double beta,
                          const std::vector<double>& x,
                          const std::vector<double>& se_x,
                          const std::vector<double>& y,
                          const std::vector<double>& se_y,
                          const Matrix& ld,
                          double sampling_corr) {
    const int n = static_cast<int>(x.size());
    Matrix covariance(n, std::vector<double>(n, 0.0));
    std::vector<double> residual(n, 0.0);
    sampling_corr = clamp_corr(sampling_corr);
    for (int i = 0; i < n; ++i) {
        residual[i] = y[i] - beta * x[i];
        for (int j = 0; j < n; ++j) {
            const double r = ld[i][j];
            const double cross = sampling_corr * r *
                (se_x[i] * se_y[j] + se_y[i] * se_x[j]);
            covariance[i][j] = r * (se_y[i] * se_y[j] +
                                     beta * beta * se_x[i] * se_x[j]) -
                               beta * cross;
        }
    }
    Matrix lower;
    double log_det = 0.0;
    if (!cholesky(covariance, lower, log_det)) {
        return -std::numeric_limits<double>::infinity();
    }
    const double quadratic = quadratic_form_from_cholesky(lower, residual);
    return -0.5 * (n * LOG2PI + log_det + quadratic);
}

FactorFit fit_eiv(const std::vector<double>& x,
                  const std::vector<double>& se_x,
                  const std::vector<double>& y,
                  const std::vector<double>& se_y,
                  const Matrix& supplied_ld,
                  double sampling_corr,
                  double prior_variance,
                  int quadrature_points) {
    FactorFit result;
    result.n = static_cast<int>(x.size());
    if (x.empty() || x.size() != y.size() || x.size() != se_x.size() ||
        x.size() != se_y.size()) return result;
    for (size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i]) ||
            !(se_x[i] > 0.0) || !(se_y[i] > 0.0)) return result;
    }
    const Matrix ld = valid_ld(supplied_ld, result.n)
        ? supplied_ld : identity_matrix(result.n);
    auto likelihood = [&](double beta) {
        return eiv_log_likelihood(beta, x, se_x, y, se_y, ld, sampling_corr);
    };

    double numerator = 0.0, denominator = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        const double weight = 1.0 / (se_y[i] * se_y[i]);
        numerator += weight * x[i] * y[i];
        denominator += weight * x[i] * x[i];
    }
    const double initial = denominator > 0.0 ? numerator / denominator : 0.0;
    const double prior_sd = std::sqrt(prior_variance);
    const double bound = std::max({2.0, 8.0 * prior_sd, 4.0 * std::fabs(initial) + 1.0});
    const int scan_points = 121;
    double best_beta = 0.0;
    double best_ll = -std::numeric_limits<double>::infinity();
    int best_index = 0;
    for (int i = 0; i < scan_points; ++i) {
        const double beta = -bound + 2.0 * bound * i / (scan_points - 1.0);
        const double ll = likelihood(beta);
        if (ll > best_ll) {
            best_ll = ll;
            best_beta = beta;
            best_index = i;
        }
    }
    double left = -bound + 2.0 * bound * std::max(0, best_index - 1) /
                              (scan_points - 1.0);
    double right = -bound + 2.0 * bound * std::min(scan_points - 1, best_index + 1) /
                               (scan_points - 1.0);
    const double golden = 0.5 * (std::sqrt(5.0) - 1.0);
    double c = right - golden * (right - left);
    double d = left + golden * (right - left);
    double fc = likelihood(c), fd = likelihood(d);
    for (int iter = 0; iter < 80; ++iter) {
        if (fc > fd) {
            right = d; d = c; fd = fc;
            c = right - golden * (right - left); fc = likelihood(c);
        } else {
            left = c; c = d; fc = fd;
            d = left + golden * (right - left); fd = likelihood(d);
        }
    }
    best_beta = 0.5 * (left + right);
    best_ll = likelihood(best_beta);

    const double h = std::max(1e-5, std::max(std::fabs(best_beta), 1.0) * 1e-4);
    const double curvature = -(likelihood(best_beta + h) - 2.0 * best_ll +
                               likelihood(best_beta - h)) / (h * h);
    if (curvature > 0.0 && std::isfinite(curvature)) {
        result.se = 1.0 / std::sqrt(curvature);
    }

    // Deterministic Simpson quadrature of L(beta) under N(0, prior_variance).
    int points = std::max(41, quadrature_points);
    if (points % 2 == 0) points++;
    const double q_bound = 8.0 * prior_sd;
    const double step = 2.0 * q_bound / (points - 1.0);
    std::vector<double> terms(points, -std::numeric_limits<double>::infinity());
    double max_term = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < points; ++i) {
        const double beta = -q_bound + i * step;
        terms[i] = likelihood(beta) - 0.5 * beta * beta / prior_variance -
                   0.5 * std::log(2.0 * M_PI * prior_variance);
        max_term = std::max(max_term, terms[i]);
    }
    double weighted_sum = 0.0;
    for (int i = 0; i < points; ++i) {
        const double weight = (i == 0 || i == points - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        weighted_sum += weight * std::exp(terms[i] - max_term);
    }
    const double log_marginal = max_term + std::log(step * weighted_sum / 3.0);
    result.beta = best_beta;
    result.p = gls_score_p(x, y, se_y, ld);
    result.log_bf = log_marginal - likelihood(0.0);
    result.valid = std::isfinite(result.beta) && std::isfinite(result.se) &&
                   std::isfinite(result.p) && std::isfinite(result.log_bf);
    return result;
}

std::string factor_pattern(bool xm, bool my, bool xy, bool pleiotropy) {
    std::string result = xm ? "XM+" : "XM-";
    result += my ? "_MY+" : "_MY-";
    result += xy ? "_XY+" : "_XY-";
    result += pleiotropy ? "_PLEIO+" : "_PLEIO-";
    return result;
}

} // namespace

void run_factorized_inference(const ProteinData& prot,
                              ProteinResult& result,
                              const Options& opts) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    result.factor_beta1 = result.factor_beta1_se = result.factor_p_xm = nan;
    result.factor_log_bf_xm = nan;
    result.factor_beta2 = result.factor_beta2_se = result.factor_p_my = nan;
    result.factor_log_bf_my = nan;
    result.factor_beta3 = result.factor_beta3_se = result.factor_p_xy = nan;
    result.factor_log_bf_xy = nan;
    result.factor_indirect = result.factor_indirect_se = nan;
    result.factor_conjunction_p = result.factor_conjunction_q_by = nan;
    result.factor_min_log_bf = nan;
    result.factor_pleiotropy_rho = result.factor_pleiotropy_p = nan;
    result.factor_nA = 0;
    result.factor_nB = prot.nB();
    result.factor_ld_source = prot.ld_reference_used ? "reference" : "identity";
    result.factor_pattern = "NOT_RUN";
    result.factor_mediation_status = "NOT_RUN";
    result.factor_frequentist_status = "NOT_RUN";
    if (opts.structural_method != "factorized") return;

    std::vector<int> keep_a;
    std::vector<double> gamma, se_gamma, alpha, se_alpha, outcome, se_outcome;
    for (int i = 0; i < prot.nA(); ++i) {
        const bool observed = i < static_cast<int>(prot.setA_alpha_observed.size()) &&
                              prot.setA_alpha_observed[i];
        const double reliability = i < static_cast<int>(prot.setA_alpha_reliability.size())
            ? prot.setA_alpha_reliability[i] : 1.0;
        // Proxy-projected protein effects have additional LD-estimation error
        // that is not represented by the reported pQTL SE. Keep them in the
        // legacy model, but exclude them from confirmatory factor inference.
        if (!observed || reliability < 0.999) continue;
        keep_a.push_back(i);
        gamma.push_back(prot.setA_gamma[i]);
        se_gamma.push_back(prot.setA_se_gamma[i]);
        alpha.push_back(prot.setA_alpha[i]);
        se_alpha.push_back(prot.setA_se_alpha[i]);
        outcome.push_back(prot.setA_Gamma[i]);
        se_outcome.push_back(prot.setA_se_Gamma[i]);
    }
    result.factor_nA = static_cast<int>(keep_a.size());
    const Matrix ld_a = valid_ld(prot.setA_ld, prot.nA())
        ? subset_matrix(prot.setA_ld, keep_a) : identity_matrix(result.factor_nA);
    const Matrix ld_b = valid_ld(prot.setB_ld, prot.nB())
        ? prot.setB_ld : identity_matrix(prot.nB());

    const FactorFit stage1 = fit_eiv(
        gamma, se_gamma, alpha, se_alpha, ld_a,
        opts.sampling_corr_rf_pqtl, opts.prior_sigma2_beta1,
        opts.factor_quadrature_points);
    const FactorFit stage2 = fit_eiv(
        prot.setB_alpha_cis, prot.setB_se_alpha_cis,
        prot.setB_Gamma_cis, prot.setB_se_Gamma_cis, ld_b,
        opts.sampling_corr_pqtl_outcome, opts.prior_sigma2_beta2,
        opts.factor_quadrature_points);

    if (stage1.valid) {
        result.factor_beta1 = stage1.beta;
        result.factor_beta1_se = stage1.se;
        result.factor_p_xm = stage1.p;
        result.factor_log_bf_xm = stage1.log_bf;
    }
    if (stage2.valid) {
        result.factor_beta2 = stage2.beta;
        result.factor_beta2_se = stage2.se;
        result.factor_p_my = stage2.p;
        result.factor_log_bf_my = stage2.log_bf;
    }

    FactorFit direct;
    if (stage2.valid && !gamma.empty()) {
        std::vector<double> residual_outcome(outcome.size());
        std::vector<double> residual_se(outcome.size());
        for (size_t i = 0; i < outcome.size(); ++i) {
            residual_outcome[i] = outcome[i] - stage2.beta * alpha[i];
            double variance = se_outcome[i] * se_outcome[i] +
                              stage2.beta * stage2.beta * se_alpha[i] * se_alpha[i] +
                              alpha[i] * alpha[i] * stage2.se * stage2.se -
                              2.0 * stage2.beta * opts.sampling_corr_pqtl_outcome *
                                  se_outcome[i] * se_alpha[i];
            residual_se[i] = std::sqrt(std::max(variance, 1e-12));
        }
        direct = fit_eiv(gamma, se_gamma, residual_outcome, residual_se, ld_a,
                         0.0, opts.prior_sigma2_beta3,
                         opts.factor_quadrature_points);
        if (direct.valid) {
            result.factor_beta3 = direct.beta;
            result.factor_beta3_se = direct.se;
            result.factor_p_xy = direct.p;
            result.factor_log_bf_xy = direct.log_bf;
        }
    }

    if (stage1.valid && stage2.valid) {
        result.factor_indirect = stage1.beta * stage2.beta;
        result.factor_indirect_se = std::sqrt(
            stage2.beta * stage2.beta * stage1.se * stage1.se +
            stage1.beta * stage1.beta * stage2.se * stage2.se);
        result.factor_conjunction_p = std::max(stage1.p, stage2.p);
        result.factor_min_log_bf = std::min(stage1.log_bf, stage2.log_bf);
    }

    // A de-noised residual-correlation diagnostic for coexisting pleiotropy.
    // It is deliberately not part of the confirmatory mediation p-value.
    if (stage1.valid && stage2.valid && direct.valid && gamma.size() >= 4) {
        double cross = 0.0, var_m = 0.0, var_y = 0.0;
        for (size_t i = 0; i < gamma.size(); ++i) {
            const double rm = alpha[i] - stage1.beta * gamma[i];
            const double ry = outcome[i] - stage2.beta * alpha[i] - direct.beta * gamma[i];
            const double vm = se_alpha[i] * se_alpha[i] +
                              stage1.beta * stage1.beta * se_gamma[i] * se_gamma[i] -
                              2.0 * stage1.beta * opts.sampling_corr_rf_pqtl *
                                  se_alpha[i] * se_gamma[i];
            const double vy = se_outcome[i] * se_outcome[i] +
                              stage2.beta * stage2.beta * se_alpha[i] * se_alpha[i] +
                              direct.beta * direct.beta * se_gamma[i] * se_gamma[i];
            const double cmy = opts.sampling_corr_pqtl_outcome * se_alpha[i] * se_outcome[i] -
                stage2.beta * se_alpha[i] * se_alpha[i] -
                direct.beta * opts.sampling_corr_rf_pqtl * se_alpha[i] * se_gamma[i] -
                stage1.beta * opts.sampling_corr_rf_outcome * se_gamma[i] * se_outcome[i] +
                stage1.beta * stage2.beta * opts.sampling_corr_rf_pqtl *
                    se_gamma[i] * se_alpha[i] +
                stage1.beta * direct.beta * se_gamma[i] * se_gamma[i];
            cross += rm * ry - cmy;
            var_m += std::max(0.0, rm * rm - vm);
            var_y += std::max(0.0, ry * ry - vy);
        }
        if (var_m > 0.0 && var_y > 0.0) {
            const double rho = clamp_corr(cross / std::sqrt(var_m * var_y));
            result.factor_pleiotropy_rho = rho;
            const double z = std::atanh(rho) * std::sqrt(gamma.size() - 3.0);
            result.factor_pleiotropy_p = std::erfc(std::fabs(z) / std::sqrt(2.0));
        }
    }

    const bool xm = stage1.valid && stage1.p <= opts.factor_alpha;
    const bool my = stage2.valid && stage2.p <= opts.factor_alpha;
    const bool xy = direct.valid && direct.p <= opts.factor_alpha;
    const bool pleio = std::isfinite(result.factor_pleiotropy_p) &&
                       result.factor_pleiotropy_p <= opts.factor_alpha;
    result.factor_pattern = factor_pattern(xm, my, xy, pleio);

    auto structural_gate = [&](bool evidence, const std::string& no_evidence,
                               const std::string& pending) {
        if (result.factor_nA < opts.factor_min_set_a && result.factor_nB < opts.factor_min_set_b)
            return std::string("INSUFFICIENT_SET_A_AND_SET_B");
        if (result.factor_nA < opts.factor_min_set_a) return std::string("INSUFFICIENT_SET_A");
        if (result.factor_nB < opts.factor_min_set_b) return std::string("INSUFFICIENT_SET_B");
        if (!stage1.valid || !stage2.valid) return std::string("NUMERICAL_FAILURE");
        if (!evidence) return no_evidence;
        if (std::fabs(opts.sampling_corr_rf_pqtl) > EPS ||
            std::fabs(opts.sampling_corr_pqtl_outcome) > EPS)
            return std::string("UNRESOLVED_SAMPLE_OVERLAP");
        if (!prot.ld_reference_used) return std::string("UNRESOLVED_NO_LD_REFERENCE");
        if (result.mediation_identifiability == "LD_DISTINCT_SUPPORTED")
            return std::string("REJECTED_DISTINCT_REGIONAL_SIGNALS");
        if (result.mediation_identifiability ==
            "LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL") return pending;
        return std::string("UNRESOLVED_REGIONAL_CONFIGURATION");
    };
    const bool bayes_evidence = std::isfinite(result.factor_min_log_bf) &&
        result.factor_min_log_bf >= std::log(opts.factor_bf_threshold);
    result.factor_mediation_status = structural_gate(
        bayes_evidence, "NO_TWO_STAGE_BAYES_EVIDENCE", "SUPPORTED_CONDITIONAL");
    const bool frequentist_evidence = std::isfinite(result.factor_conjunction_p) &&
        result.factor_conjunction_p <= opts.factor_alpha;
    result.factor_frequentist_status = structural_gate(
        frequentist_evidence, "NO_TWO_STAGE_FREQUENTIST_EVIDENCE",
        "PENDING_MULTIPLE_TESTING");
}

void finalize_factorized_multiple_testing(std::vector<ProteinResult>& results,
                                          const Options& opts) {
    if (opts.structural_method != "factorized") return;
    std::vector<int> order;
    for (int i = 0; i < static_cast<int>(results.size()); ++i) {
        if (std::isfinite(results[i].factor_conjunction_p)) order.push_back(i);
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return results[a].factor_conjunction_p < results[b].factor_conjunction_p;
    });
    const int m = static_cast<int>(order.size());
    double harmonic = 0.0;
    for (int i = 1; i <= m; ++i) harmonic += 1.0 / i;
    double running = 1.0;
    for (int rank = m; rank >= 1; --rank) {
        ProteinResult& result = results[order[rank - 1]];
        const double raw = result.factor_conjunction_p * m * harmonic / rank;
        running = std::min(running, std::min(1.0, raw));
        result.factor_conjunction_q_by = running;
    }
    for (auto& result : results) {
        if (result.factor_frequentist_status != "PENDING_MULTIPLE_TESTING") continue;
        result.factor_frequentist_status = result.factor_conjunction_q_by <= opts.factor_alpha
            ? "SUPPORTED_CONDITIONAL" : "NOT_SELECTED_BY_FDR";
    }
}

} // namespace bmediator
