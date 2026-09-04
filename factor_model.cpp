#include "bmediator.h"

#include <functional>

namespace bmediator {
namespace {

using Matrix = std::vector<std::vector<double>>;

struct FactorFit {
    double beta = std::numeric_limits<double>::quiet_NaN();
    double se = std::numeric_limits<double>::quiet_NaN();
    double p = std::numeric_limits<double>::quiet_NaN();
    double p_strict = std::numeric_limits<double>::quiet_NaN();
    double log_bf = std::numeric_limits<double>::quiet_NaN();
    double log_bf_heterogeneity = std::numeric_limits<double>::quiet_NaN();
    double log_bf_directional = std::numeric_limits<double>::quiet_NaN();
    double log_bf_slope_only = std::numeric_limits<double>::quiet_NaN();
    double log_bf_directional_only = std::numeric_limits<double>::quiet_NaN();
    double log_bf_slope_directional = std::numeric_limits<double>::quiet_NaN();
    double pp_slope = std::numeric_limits<double>::quiet_NaN();
    double pp_directional = std::numeric_limits<double>::quiet_NaN();
    double directional_intercept = std::numeric_limits<double>::quiet_NaN();
    double directional_intercept_se = std::numeric_limits<double>::quiet_NaN();
    double directional_collinearity = std::numeric_limits<double>::quiet_NaN();
    double log_e = std::numeric_limits<double>::quiet_NaN();
    double tau = std::numeric_limits<double>::quiet_NaN();
    int n = 0;
    bool effect_valid = false;
    bool valid = false;
};

struct ScoreFit {
    double p_robust = std::numeric_limits<double>::quiet_NaN();
    double p_strict = std::numeric_limits<double>::quiet_NaN();
    double dispersion = std::numeric_limits<double>::quiet_NaN();
    double statistic = std::numeric_limits<double>::quiet_NaN();
    double log_e = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

double clamp_corr(double value) {
    return std::max(-0.999, std::min(0.999, value));
}

double log_sum_exp(const std::vector<double>& values) {
    if (values.empty()) return -std::numeric_limits<double>::infinity();
    const double maximum = *std::max_element(values.begin(), values.end());
    if (!std::isfinite(maximum)) return maximum;
    double total = 0.0;
    for (double value : values) total += std::exp(value - maximum);
    return maximum + std::log(total);
}

double log_weighted_pair(double first, double first_weight,
                         double second, double second_weight) {
    return log_sum_exp({first + std::log(first_weight),
                        second + std::log(second_weight)});
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

double student_t_critical(double two_sided_alpha, double df) {
    if (!(two_sided_alpha > 0.0 && two_sided_alpha < 1.0) || !(df > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double lower = 0.0;
    double upper = 2.0;
    while (two_sided_student_t_p(upper, df) > two_sided_alpha && upper < 1e6) {
        upper *= 2.0;
    }
    for (int iter = 0; iter < 80; ++iter) {
        const double middle = 0.5 * (lower + upper);
        if (two_sided_student_t_p(middle, df) > two_sided_alpha) lower = middle;
        else upper = middle;
    }
    return 0.5 * (lower + upper);
}

double student_t_log_density(double statistic, double df) {
    if (!(df > 0.0) || !std::isfinite(statistic)) {
        return -std::numeric_limits<double>::infinity();
    }
    return std::lgamma(0.5 * (df + 1.0)) - std::lgamma(0.5 * df) -
           0.5 * std::log(df * M_PI) -
           0.5 * (df + 1.0) * std::log1p(statistic * statistic / df);
}

// A fixed proper alternative density on the exact nuisance-free t statistic.
// Every shifted or scaled component integrates to one, so the density ratio
// has null expectation one for every common Gaussian scale and oriented
// intercept. The symmetric grid is prespecified and never fit to analyzed data.
double student_t_mixture_log_e(double statistic, double df) {
    const double log_null = student_t_log_density(statistic, df);
    if (!std::isfinite(log_null)) return std::numeric_limits<double>::quiet_NaN();
    std::vector<double> terms;
    for (double shift : {2.0, 4.0, 6.0}) {
        terms.push_back(log_weighted_pair(
            student_t_log_density(statistic - shift, df), 0.5,
            student_t_log_density(statistic + shift, df), 0.5));
    }
    for (double scale : {2.0, 4.0, 8.0}) {
        terms.push_back(student_t_log_density(statistic / scale, df) -
                        std::log(scale));
    }
    return log_sum_exp(terms) - std::log(static_cast<double>(terms.size())) -
           log_null;
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

// Null score test conditional on independently selected exposure associations.
// The oriented intercept is an unrestricted nuisance term, so the test asks
// whether a proportional slope remains after allele-oriented pleiotropy is
// projected out. Under independent exposure/outcome errors, the strict score
// is Gaussian under H0 even when the exposure associations are weak.
ScoreFit gls_score(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& se_y,
                   const Matrix& ld) {
    ScoreFit result;
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
        return result;
    }
    const auto wx = whiten(lower, x);
    const auto wy = whiten(lower, y);
    std::vector<double> orientation(n, 1.0);
    for (int i = 0; i < n; ++i) orientation[i] = x[i] < 0.0 ? -1.0 : 1.0;
    const auto ws = whiten(lower, orientation);
    const double ss = std::inner_product(ws.begin(), ws.end(), ws.begin(), 0.0);
    const double xs = std::inner_product(wx.begin(), wx.end(), ws.begin(), 0.0);
    const double ys = std::inner_product(wy.begin(), wy.end(), ws.begin(), 0.0);
    if (!(ss > 0.0) || !std::isfinite(ss) || !std::isfinite(xs) ||
        !std::isfinite(ys)) {
        return result;
    }
    std::vector<double> rx(n), ry(n);
    for (int i = 0; i < n; ++i) {
        rx[i] = wx[i] - (xs / ss) * ws[i];
        ry[i] = wy[i] - (ys / ss) * ws[i];
    }
    const double q = std::inner_product(rx.begin(), rx.end(), rx.begin(), 0.0);
    const double t = std::inner_product(rx.begin(), rx.end(), ry.begin(), 0.0);
    const double total_x = std::inner_product(wx.begin(), wx.end(), wx.begin(), 0.0);
    if (!std::isfinite(q) || !std::isfinite(t) ||
        q <= 1e-12 * std::max(1.0, total_x)) {
        result.p_strict = 1.0;
        result.p_robust = 1.0;
        result.dispersion = 1.0;
        result.statistic = 0.0;
        result.log_e = n > 2 ? student_t_mixture_log_e(0.0, n - 2.0) :
            std::numeric_limits<double>::quiet_NaN();
        result.valid = std::isfinite(result.log_e);
        return result;
    }
    const double beta = t / q;
    double residual_q = 0.0;
    for (int i = 0; i < n; ++i) {
        const double residual = ry[i] - beta * rx[i];
        residual_q += residual * residual;
    }
    const double df = n - 2.0;
    const double estimated_dispersion = df > 0.0 ? residual_q / df :
        std::numeric_limits<double>::quiet_NaN();
    const double dispersion = std::max(estimated_dispersion, 1e-300);
    const double strict_z = t / std::sqrt(q);
    const double robust_z = strict_z / std::sqrt(dispersion);
    result.p_strict = std::erfc(std::fabs(strict_z) / std::sqrt(2.0));
    result.p_robust = df > 0.0
        ? two_sided_student_t_p(robust_z, df)
        : result.p_strict;
    result.dispersion = dispersion;
    result.statistic = robust_z;
    result.log_e = student_t_mixture_log_e(robust_z, df);
    result.valid = std::isfinite(result.p_strict) &&
                   std::isfinite(result.p_robust) && std::isfinite(result.log_e);
    return result;
}

double p_to_e_log_mixture(double p_value) {
    if (!(p_value > 0.0) || p_value > 1.0 || !std::isfinite(p_value)) {
        if (p_value == 0.0) p_value = std::numeric_limits<double>::min();
        else return std::numeric_limits<double>::quiet_NaN();
    }
    const std::vector<double> kappas = {0.10, 0.25, 0.50, 0.75};
    std::vector<double> terms;
    terms.reserve(kappas.size());
    for (double kappa : kappas) {
        terms.push_back(std::log(kappa) + (kappa - 1.0) * std::log(p_value));
    }
    return log_sum_exp(terms) - std::log(static_cast<double>(kappas.size()));
}

double posterior_from_log_bf(double log_bf, double prior_probability) {
    if (!std::isfinite(log_bf) || !(prior_probability > 0.0) ||
        !(prior_probability < 1.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double log_odds = log_bf + std::log(prior_probability) -
                            std::log1p(-prior_probability);
    if (log_odds >= 0.0) {
        return 1.0 / (1.0 + std::exp(-log_odds));
    }
    const double odds = std::exp(log_odds);
    return odds / (1.0 + odds);
}

double eiv_log_likelihood(double beta,
                          const std::vector<double>& x,
                          const std::vector<double>& se_x,
                          const std::vector<double>& y,
                          const std::vector<double>& se_y,
                          const Matrix& ld,
                          double sampling_corr,
                          double pleiotropy_sd) {
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
            if (i == j) covariance[i][j] += pleiotropy_sd * pleiotropy_sd;
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

double eiv_directional_marginal_log_likelihood(
        double beta,
        const std::vector<double>& x,
        const std::vector<double>& se_x,
        const std::vector<double>& y,
        const std::vector<double>& se_y,
        const Matrix& ld,
        double sampling_corr,
        double pleiotropy_sd,
        double directional_variance) {
    const int n = static_cast<int>(x.size());
    Matrix covariance(n, std::vector<double>(n, 0.0));
    std::vector<double> residual(n, 0.0);
    std::vector<double> orientation(n, 1.0);
    sampling_corr = clamp_corr(sampling_corr);
    for (int i = 0; i < n; ++i) {
        residual[i] = y[i] - beta * x[i];
        orientation[i] = x[i] < 0.0 ? -1.0 : 1.0;
        for (int j = 0; j < n; ++j) {
            const double r = ld[i][j];
            const double cross = sampling_corr * r *
                (se_x[i] * se_y[j] + se_y[i] * se_x[j]);
            covariance[i][j] = r * (se_y[i] * se_y[j] +
                                     beta * beta * se_x[i] * se_x[j]) -
                               beta * cross;
            if (i == j) covariance[i][j] += pleiotropy_sd * pleiotropy_sd;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            covariance[i][j] += directional_variance *
                                orientation[i] * orientation[j];
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

// Profile likelihood for beta after eliminating the variant-specific true
// exposure associations. The residual determinant is omitted when beta is
// optimized; retaining it creates weak-instrument attenuation. For fixed beta,
// the full residual likelihood remains appropriate for estimating tau.
double eiv_adjusted_profile_log_likelihood(
        double beta,
        const std::vector<double>& x,
        const std::vector<double>& se_x,
        const std::vector<double>& y,
        const std::vector<double>& se_y,
        const Matrix& ld,
        double sampling_corr,
        double pleiotropy_sd) {
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
            if (i == j) covariance[i][j] += pleiotropy_sd * pleiotropy_sd;
        }
    }
    Matrix lower;
    double log_det = 0.0;
    if (!cholesky(covariance, lower, log_det)) {
        return -std::numeric_limits<double>::infinity();
    }
    return -0.5 * quadratic_form_from_cholesky(lower, residual);
}

double eiv_fixed_directional_log_likelihood(
        double beta,
        double eta,
        const std::vector<double>& x,
        const std::vector<double>& se_x,
        const std::vector<double>& y,
        const std::vector<double>& se_y,
        const Matrix& ld,
        double sampling_corr,
        double pleiotropy_sd,
        bool include_determinant) {
    const int n = static_cast<int>(x.size());
    Matrix covariance(n, std::vector<double>(n, 0.0));
    std::vector<double> residual(n, 0.0);
    sampling_corr = clamp_corr(sampling_corr);
    for (int i = 0; i < n; ++i) {
        const double orientation = x[i] < 0.0 ? -1.0 : 1.0;
        residual[i] = y[i] - beta * x[i] - eta * orientation;
        for (int j = 0; j < n; ++j) {
            const double r = ld[i][j];
            const double cross = sampling_corr * r *
                (se_x[i] * se_y[j] + se_y[i] * se_x[j]);
            covariance[i][j] = r * (se_y[i] * se_y[j] +
                                     beta * beta * se_x[i] * se_x[j]) -
                               beta * cross;
            if (i == j) covariance[i][j] += pleiotropy_sd * pleiotropy_sd;
        }
    }
    Matrix lower;
    double log_det = 0.0;
    if (!cholesky(covariance, lower, log_det)) {
        return -std::numeric_limits<double>::infinity();
    }
    const double quadratic = quadratic_form_from_cholesky(lower, residual);
    if (!include_determinant) return -0.5 * quadratic;
    return -0.5 * (n * LOG2PI + log_det + quadratic);
}

FactorFit fit_eiv(const std::vector<double>& x,
                  const std::vector<double>& se_x,
                  const std::vector<double>& y,
                  const std::vector<double>& se_y,
                  const Matrix& supplied_ld,
                  double sampling_corr,
                  double prior_variance,
                  int quadrature_points,
                  double pleiotropy_prior_sd,
                  int pleiotropy_quadrature_points,
                  double slope_prior_probability,
                  double directional_prior_probability,
                  double directional_variance) {
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
    Matrix null_covariance(result.n, std::vector<double>(result.n, 0.0));
    for (int i = 0; i < result.n; ++i) {
        for (int j = 0; j < result.n; ++j) {
            null_covariance[i][j] = ld[i][j] * se_y[i] * se_y[j];
        }
    }
    Matrix null_lower;
    double null_log_det = 0.0;
    std::vector<double> orientation(result.n, 1.0);
    for (int i = 0; i < result.n; ++i) {
        orientation[i] = x[i] < 0.0 ? -1.0 : 1.0;
    }
    if (cholesky(null_covariance, null_lower, null_log_det)) {
        const std::vector<double> wx = whiten(null_lower, x);
        const std::vector<double> ws = whiten(null_lower, orientation);
        const double cross = std::inner_product(wx.begin(), wx.end(), ws.begin(), 0.0);
        const double xx = std::inner_product(wx.begin(), wx.end(), wx.begin(), 0.0);
        const double ss = std::inner_product(ws.begin(), ws.end(), ws.begin(), 0.0);
        if (xx > 0.0 && ss > 0.0) {
            result.directional_collinearity =
                std::min(1.0, std::fabs(cross) / std::sqrt(xx * ss));
        }
    }
    auto likelihood = [&](double beta, double tau) {
        return eiv_log_likelihood(beta, x, se_x, y, se_y, ld, sampling_corr, tau);
    };
    auto directional_likelihood = [&](double beta, double tau) {
        return eiv_directional_marginal_log_likelihood(
            beta, x, se_x, y, se_y, ld, sampling_corr, tau,
            directional_variance);
    };
    bool independent_errors = std::fabs(sampling_corr) < 1e-12;
    for (int i = 0; i < result.n; ++i) {
        for (int j = 0; j < result.n; ++j) {
            if (i != j && std::fabs(ld[i][j]) > 1e-12) {
                independent_errors = false;
            }
        }
    }

    double numerator = 0.0, denominator = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        const double weight = 1.0 / (se_y[i] * se_y[i]);
        numerator += weight * x[i] * y[i];
        denominator += weight * x[i] * x[i];
    }
    const double initial = denominator > 0.0 ? numerator / denominator : 0.0;
    const double prior_sd = std::sqrt(prior_variance);
    const double bound = std::max({2.0, 8.0 * prior_sd, 4.0 * std::fabs(initial) + 1.0});
    const double golden = 0.5 * (std::sqrt(5.0) - 1.0);
    auto maximize_interval = [&](const std::function<double(double)>& objective,
                                 double lower, double upper, int scan_points) {
        double best_value = -std::numeric_limits<double>::infinity();
        int best_index = 0;
        double best_grid_point = lower;
        for (int i = 0; i < scan_points; ++i) {
            const double value = lower + (upper - lower) * i / (scan_points - 1.0);
            const double objective_value = objective(value);
            if (objective_value > best_value) {
                best_value = objective_value;
                best_index = i;
                best_grid_point = value;
            }
        }
        double left = lower + (upper - lower) * std::max(0, best_index - 1) /
                                  (scan_points - 1.0);
        double right = lower + (upper - lower) * std::min(scan_points - 1, best_index + 1) /
                                   (scan_points - 1.0);
        double c = right - golden * (right - left);
        double d = left + golden * (right - left);
        double fc = objective(c), fd = objective(d);
        for (int iter = 0; iter < 60; ++iter) {
            if (fc > fd) {
                right = d; d = c; fd = fc;
                c = right - golden * (right - left); fc = objective(c);
            } else {
                left = c; c = d; fc = fd;
                d = left + golden * (right - left); fd = objective(d);
            }
        }
        const double estimate = 0.5 * (left + right);
        const double estimate_value = objective(estimate);
        return best_value >= estimate_value
            ? std::make_pair(best_grid_point, best_value)
            : std::make_pair(estimate, estimate_value);
    };

    const double tau_bound = 8.0 * pleiotropy_prior_sd;
    double max_se_y = 0.0;
    for (double value : se_y) max_se_y = std::max(max_se_y, value);
    const double point_tau_bound = std::max(tau_bound, 10.0 * max_se_y);
    auto point_beta_objective = [&](double beta, double tau) {
        return eiv_adjusted_profile_log_likelihood(
            beta, x, se_x, y, se_y, ld, sampling_corr, tau);
    };
    auto point_tau_objective = [&](double beta, double tau) {
        if (!independent_errors) return likelihood(beta, tau);
        double value = 0.0;
        const double tau2 = tau * tau;
        for (size_t i = 0; i < x.size(); ++i) {
            const double residual = y[i] - beta * x[i];
            const double variance = se_y[i] * se_y[i] +
                                    beta * beta * se_x[i] * se_x[i] + tau2;
            const double weight = se_x[i] * se_x[i];
            value -= 0.5 * weight *
                     (std::log(variance) + residual * residual / variance);
        }
        return value;
    };
    auto optimize_beta = [&](double tau) {
        return maximize_interval(
            [&](double beta) { return point_beta_objective(beta, tau); },
            -bound, bound, 81);
    };
    auto optimize_tau = [&](double beta) {
        return maximize_interval(
            [&](double tau) { return point_tau_objective(beta, tau); },
            0.0, point_tau_bound, 51);
    };

    double best_tau = 0.0;
    double best_beta = initial;
    double best_ll = -std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < 12; ++iter) {
        const auto beta_fit = optimize_beta(best_tau);
        best_beta = beta_fit.first;
        const auto tau_fit = optimize_tau(best_beta);
        const double previous_ll = best_ll;
        best_tau = tau_fit.first;
        best_ll = tau_fit.second;
        if (std::isfinite(previous_ll) && std::fabs(best_ll - previous_ll) < 1e-9) break;
    }
    const auto final_beta_fit = optimize_beta(best_tau);
    best_beta = final_beta_fit.first;
    const auto final_tau_fit = optimize_tau(best_beta);
    best_tau = final_tau_fit.first;
    best_ll = final_tau_fit.second;
    result.tau = best_tau;

    // Joint adjusted-profile slope and oriented intercept. This aligns the
    // reported point estimate with the nuisance-projected score test.
    bool joint_fit = false;
    double joint_eta = 0.0;
    double joint_beta_se = std::numeric_limits<double>::quiet_NaN();
    double joint_eta_se = std::numeric_limits<double>::quiet_NaN();
    if (result.n >= 3 && std::isfinite(result.directional_collinearity) &&
        result.directional_collinearity < 1.0 - 1e-10) {
        double max_abs_y = 0.0;
        for (double value : y) max_abs_y = std::max(max_abs_y, std::fabs(value));
        const double eta_bound = std::max(
            1.0, 4.0 * max_abs_y + 8.0 * std::sqrt(directional_variance));
        auto joint_adjusted = [&](double beta, double eta, double tau) {
            return eiv_fixed_directional_log_likelihood(
                beta, eta, x, se_x, y, se_y, ld, sampling_corr, tau, false);
        };
        auto joint_full = [&](double beta, double eta, double tau) {
            return eiv_fixed_directional_log_likelihood(
                beta, eta, x, se_x, y, se_y, ld, sampling_corr, tau, true);
        };
        double jb = best_beta;
        double je = 0.0;
        double jt = best_tau;
        auto optimize_joint_beta = [&](double eta, double tau) {
            return maximize_interval(
                [&](double beta) { return joint_adjusted(beta, eta, tau); },
                -bound, bound, 81);
        };
        auto optimize_joint_eta = [&](double beta, double tau) {
            return maximize_interval(
                [&](double eta) { return joint_adjusted(beta, eta, tau); },
                -eta_bound, eta_bound, 81);
        };
        auto optimize_joint_tau = [&](double beta, double eta) {
            return maximize_interval(
                [&](double tau) { return joint_full(beta, eta, tau); },
                0.0, point_tau_bound, 51);
        };
        for (int iter = 0; iter < 20; ++iter) {
            const double previous_beta = jb;
            const double previous_eta = je;
            jb = optimize_joint_beta(je, jt).first;
            je = optimize_joint_eta(jb, jt).first;
            jt = optimize_joint_tau(jb, je).first;
            if (std::fabs(jb - previous_beta) < 1e-9 &&
                std::fabs(je - previous_eta) < 1e-9) break;
        }
        jb = optimize_joint_beta(je, jt).first;
        je = optimize_joint_eta(jb, jt).first;
        jt = optimize_joint_tau(jb, je).first;

        const double hb = std::max(1e-5, std::max(std::fabs(jb), 1.0) * 1e-4);
        const double he = std::max(1e-5, std::max(std::fabs(je), 1.0) * 1e-4);
        auto joint_profile = [&](double beta, double eta) {
            const double tau = optimize_joint_tau(beta, eta).first;
            return joint_adjusted(beta, eta, tau);
        };
        const double center = joint_profile(jb, je);
        const double hbb = -(joint_profile(jb + hb, je) - 2.0 * center +
                             joint_profile(jb - hb, je)) / (hb * hb);
        const double hee = -(joint_profile(jb, je + he) - 2.0 * center +
                             joint_profile(jb, je - he)) / (he * he);
        const double hbe = -(joint_profile(jb + hb, je + he) -
                             joint_profile(jb + hb, je - he) -
                             joint_profile(jb - hb, je + he) +
                             joint_profile(jb - hb, je - he)) /
                           (4.0 * hb * he);
        const double determinant = hbb * hee - hbe * hbe;
        if (hbb > 0.0 && hee > 0.0 && determinant > 0.0 &&
            std::isfinite(determinant)) {
            joint_beta_se = std::sqrt(hee / determinant);
            joint_eta_se = std::sqrt(hbb / determinant);
            best_beta = jb;
            best_tau = jt;
            joint_eta = je;
            result.tau = jt;
            joint_fit = std::isfinite(joint_beta_se) &&
                        std::isfinite(joint_eta_se);
        }
    }

    if (joint_fit) {
        result.se = joint_beta_se;
        result.directional_intercept = joint_eta;
        result.directional_intercept_se = joint_eta_se;
    } else if (independent_errors) {
        double score_variance = 0.0;
        double sensitivity = 0.0;
        const double tau2 = best_tau * best_tau;
        for (size_t i = 0; i < x.size(); ++i) {
            const double sx2 = se_x[i] * se_x[i];
            const double outcome_variance = tau2 + se_y[i] * se_y[i];
            const double variance = outcome_variance + best_beta * best_beta * sx2;
            const double corrected_x2 = x[i] * x[i] - sx2;
            const double corrected_y2 = y[i] * y[i] - outcome_variance;
            score_variance +=
                (corrected_x2 * outcome_variance + corrected_y2 * sx2 +
                 sx2 * outcome_variance) / (variance * variance);
            sensitivity +=
                (corrected_x2 * outcome_variance + corrected_y2 * sx2) /
                (variance * variance);
        }
        if (score_variance > 0.0 && std::fabs(sensitivity) > 1e-12 &&
            std::isfinite(score_variance) && std::isfinite(sensitivity)) {
            result.se = std::sqrt(score_variance) / std::fabs(sensitivity);
        }
    } else {
        const double h = std::max(1e-5,
            std::max(std::fabs(best_beta), 1.0) * 1e-4);
        const auto adjusted_profile_likelihood = [&](double beta) {
            const double tau = optimize_tau(beta).first;
            return point_beta_objective(beta, tau);
        };
        const double center = adjusted_profile_likelihood(best_beta);
        const double curvature = -(
            adjusted_profile_likelihood(best_beta + h) - 2.0 * center +
            adjusted_profile_likelihood(best_beta - h)) / (h * h);
        if (curvature > 0.0 && std::isfinite(curvature)) {
            result.se = 1.0 / std::sqrt(curvature);
        }
    }

    // The same half-normal pleiotropy prior is integrated under H0 and H1.
    // This prevents residual heterogeneity from becoming evidence for beta.
    int points = std::max(41, quadrature_points);
    if (points % 2 == 0) points++;
    int tau_points = std::max(21, pleiotropy_quadrature_points);
    if (tau_points % 2 == 0) tau_points++;
    const double q_bound = 8.0 * prior_sd;
    const double beta_step = 2.0 * q_bound / (points - 1.0);
    const double tau_step = tau_bound / (tau_points - 1.0);
    const double log_half_normal_constant =
        0.5 * std::log(2.0 / M_PI) - std::log(pleiotropy_prior_sd);
    std::vector<double> beta_log_weights(points);
    std::vector<double> tau_log_weights(tau_points);
    for (int i = 0; i < points; ++i) {
        const double beta = -q_bound + i * beta_step;
        const double simpson_weight = (i == 0 || i == points - 1)
            ? 1.0 : (i % 2 ? 4.0 : 2.0);
        beta_log_weights[i] = std::log(simpson_weight) -
            0.5 * beta * beta / prior_variance -
            0.5 * std::log(2.0 * M_PI * prior_variance);
    }
    for (int j = 0; j < tau_points; ++j) {
        const double tau = j * tau_step;
        const double simpson_weight = (j == 0 || j == tau_points - 1)
            ? 1.0 : (j % 2 ? 4.0 : 2.0);
        tau_log_weights[j] = std::log(simpson_weight) + log_half_normal_constant -
            0.5 * tau * tau / (pleiotropy_prior_sd * pleiotropy_prior_sd);
    }
    const double beta_weight_normalizer = log_sum_exp(beta_log_weights);
    const double tau_weight_normalizer = log_sum_exp(tau_log_weights);
    for (double& value : beta_log_weights) value -= beta_weight_normalizer;
    for (double& value : tau_log_weights) value -= tau_weight_normalizer;

    std::vector<double> alternative_terms(points * tau_points,
        -std::numeric_limits<double>::infinity());
    std::vector<double> null_terms(tau_points,
        -std::numeric_limits<double>::infinity());
    std::vector<double> directional_terms(tau_points,
        -std::numeric_limits<double>::infinity());
    std::vector<double> slope_directional_terms(points * tau_points,
        -std::numeric_limits<double>::infinity());
    for (int i = 0; i < points; ++i) {
        const double beta = -q_bound + i * beta_step;
        for (int j = 0; j < tau_points; ++j) {
            const double tau = j * tau_step;
            const double value = likelihood(beta, tau) + beta_log_weights[i] +
                                 tau_log_weights[j];
            alternative_terms[i * tau_points + j] = value;
            slope_directional_terms[i * tau_points + j] =
                directional_likelihood(beta, tau) + beta_log_weights[i] +
                tau_log_weights[j];
            if (i == 0) {
                null_terms[j] = likelihood(0.0, tau) + tau_log_weights[j];
                directional_terms[j] = directional_likelihood(0.0, tau) +
                    tau_log_weights[j];
            }
        }
    }
    const double log_alternative = log_sum_exp(alternative_terms);
    const double log_null = log_sum_exp(null_terms);
    const double log_directional = log_sum_exp(directional_terms);
    const double log_slope_directional = log_sum_exp(slope_directional_terms);
    std::vector<double> no_heterogeneity_terms(points,
        -std::numeric_limits<double>::infinity());
    for (int i = 0; i < points; ++i) {
        const double beta = -q_bound + i * beta_step;
        no_heterogeneity_terms[i] = likelihood(beta, 0.0) + beta_log_weights[i];
    }
    const double log_no_heterogeneity = log_sum_exp(no_heterogeneity_terms);
    result.beta = best_beta;
    const ScoreFit score = gls_score(x, y, se_y, ld);
    result.p = score.p_robust;
    result.p_strict = score.p_strict;
    const double log_slope_marginal = log_weighted_pair(
        log_alternative, 1.0 - directional_prior_probability,
        log_slope_directional, directional_prior_probability);
    const double log_no_slope_marginal = log_weighted_pair(
        log_null, 1.0 - directional_prior_probability,
        log_directional, directional_prior_probability);
    const double log_directional_marginal = log_weighted_pair(
        log_directional, 1.0 - slope_prior_probability,
        log_slope_directional, slope_prior_probability);
    const double log_no_directional_marginal = log_weighted_pair(
        log_null, 1.0 - slope_prior_probability,
        log_alternative, slope_prior_probability);
    result.log_bf = log_slope_marginal - log_no_slope_marginal;
    result.log_bf_directional =
        log_directional_marginal - log_no_directional_marginal;
    result.log_bf_slope_only = log_alternative - log_null;
    result.log_bf_directional_only = log_directional - log_null;
    result.log_bf_slope_directional = log_slope_directional - log_null;
    const std::vector<double> model_log_weights = {
        log_null + std::log1p(-slope_prior_probability) +
            std::log1p(-directional_prior_probability),
        log_alternative + std::log(slope_prior_probability) +
            std::log1p(-directional_prior_probability),
        log_directional + std::log1p(-slope_prior_probability) +
            std::log(directional_prior_probability),
        log_slope_directional + std::log(slope_prior_probability) +
            std::log(directional_prior_probability)};
    const double log_model_normalizer = log_sum_exp(model_log_weights);
    result.pp_slope = std::exp(log_sum_exp({model_log_weights[1],
                                            model_log_weights[3]}) -
                               log_model_normalizer);
    result.pp_directional = std::exp(log_sum_exp({model_log_weights[2],
                                                  model_log_weights[3]}) -
                                     log_model_normalizer);
    result.log_bf_heterogeneity = log_alternative - log_no_heterogeneity;
    result.log_e = score.log_e;
    result.effect_valid = joint_fit && std::isfinite(result.beta) &&
                          std::isfinite(result.se);
    result.valid = std::isfinite(result.p) && std::isfinite(result.p_strict) &&
                   std::isfinite(result.log_bf) &&
                   std::isfinite(result.log_bf_heterogeneity) &&
                   std::isfinite(result.log_e) && std::isfinite(result.tau) &&
                   std::isfinite(result.pp_slope);
    return result;
}

std::string factor_pattern(bool xm, bool my, bool xy, bool pleiotropy) {
    std::string result = xm ? "XM+" : "XM-";
    result += my ? "_MY+" : "_MY-";
    result += xy ? "_XY+" : "_XY-";
    result += pleiotropy ? "_PLEIO+" : "_PLEIO-";
    return result;
}

std::string identification_gate(const ProteinResult& result,
                                const Options& opts,
                                const std::string& evidence_status,
                                const std::string& supported_status) {
    if (evidence_status != "TWO_STAGE_EVIDENCE") return evidence_status;
    if (result.factor_ld_source != "reference") {
        return "UNRESOLVED_NO_LD_REFERENCE";
    }
    if (!std::isfinite(result.factor_cross_set_max_r2) ||
        result.factor_cross_set_max_r2 > opts.r2_thresh + 1e-12) {
        return "UNRESOLVED_CROSS_SET_LD";
    }
    if (result.mediation_identifiability == "LD_DISTINCT_SUPPORTED") {
        return "REJECTED_DISTINCT_REGIONAL_SIGNALS";
    }
    if (result.mediation_identifiability == "UNRESOLVED_SINGLE_SHARED_SIGNAL") {
        return "UNRESOLVED_SINGLE_SHARED_SIGNAL";
    }
    if (result.mediation_identifiability !=
        "OVERIDENTIFIED_SHARED_SIGNALS_ASSUMPTION_CONDITIONAL") {
        return "UNRESOLVED_REGIONAL_CONFIGURATION";
    }
    if (!std::isfinite(result.factor_beta1) ||
        !std::isfinite(result.factor_beta1_se) ||
        !std::isfinite(result.factor_beta2) ||
        !std::isfinite(result.factor_beta2_se)) {
        return "UNRESOLVED_EFFECT_ESTIMATION";
    }
    const bool xm_heterogeneity =
        std::isfinite(result.factor_log_bf_heterogeneity_xm) &&
        result.factor_log_bf_heterogeneity_xm >=
            std::log(opts.factor_bf_threshold);
    const bool my_heterogeneity =
        std::isfinite(result.factor_log_bf_heterogeneity_my) &&
        result.factor_log_bf_heterogeneity_my >=
            std::log(opts.factor_bf_threshold);
    const bool xm_directional =
        std::isfinite(result.factor_log_bf_directional_xm) &&
        result.factor_log_bf_directional_xm >=
            std::log(opts.factor_bf_threshold);
    const bool my_directional =
        std::isfinite(result.factor_log_bf_directional_my) &&
        result.factor_log_bf_directional_my >=
            std::log(opts.factor_bf_threshold);
    if (xm_directional && my_directional) {
        return "UNRESOLVED_XM_AND_MY_DIRECTIONAL_PLEIOTROPY";
    }
    if (xm_directional) return "UNRESOLVED_XM_DIRECTIONAL_PLEIOTROPY";
    if (my_directional) return "UNRESOLVED_MY_DIRECTIONAL_PLEIOTROPY";
    if (xm_heterogeneity && my_heterogeneity) {
        return "UNRESOLVED_XM_AND_MY_HETEROGENEITY";
    }
    if (xm_heterogeneity) return "UNRESOLVED_XM_HETEROGENEITY";
    if (my_heterogeneity) return "UNRESOLVED_MY_HETEROGENEITY";
    return supported_status;
}

} // namespace

void run_factorized_inference(const ProteinData& prot,
                              ProteinResult& result,
                              const Options& opts) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    result.factor_beta1 = result.factor_beta1_se = result.factor_p_xm = nan;
    result.factor_beta1_ci_lower = result.factor_beta1_ci_upper = nan;
    result.factor_p_xm_strict = nan;
    result.factor_log_bf_xm = nan;
    result.factor_log_bf_heterogeneity_xm = nan;
    result.factor_log_bf_directional_xm = nan;
    result.factor_log_bf_slope_only_xm = nan;
    result.factor_log_bf_directional_only_xm = nan;
    result.factor_log_bf_slope_directional_xm = nan;
    result.factor_pp_directional_xm = nan;
    result.factor_directional_intercept_xm = nan;
    result.factor_directional_intercept_xm_se = nan;
    result.factor_directional_collinearity_xm = nan;
    result.factor_beta2 = result.factor_beta2_se = result.factor_p_my = nan;
    result.factor_beta2_ci_lower = result.factor_beta2_ci_upper = nan;
    result.factor_p_my_strict = nan;
    result.factor_log_bf_my = nan;
    result.factor_log_bf_heterogeneity_my = nan;
    result.factor_log_bf_directional_my = nan;
    result.factor_log_bf_slope_only_my = nan;
    result.factor_log_bf_directional_only_my = nan;
    result.factor_log_bf_slope_directional_my = nan;
    result.factor_pp_directional_my = nan;
    result.factor_directional_intercept_my = nan;
    result.factor_directional_intercept_my_se = nan;
    result.factor_directional_collinearity_my = nan;
    result.factor_cross_set_max_r2 = prot.setAB_max_r2;
    result.factor_beta3 = result.factor_beta3_se = result.factor_p_xy = nan;
    result.factor_p_xy_strict = nan;
    result.factor_log_bf_xy = nan;
    result.factor_log_bf_heterogeneity_xy = nan;
    result.factor_log_e_xm = result.factor_log_e_my = result.factor_log_e_xy = nan;
    result.factor_log_e_mediation = result.factor_e_q_ebh = nan;
    result.factor_tau_xm = result.factor_tau_my = result.factor_tau_xy = nan;
    result.factor_indirect = result.factor_indirect_se = nan;
    result.factor_indirect_ci_lower = result.factor_indirect_ci_upper = nan;
    result.factor_conjunction_p = result.factor_conjunction_q_by = nan;
    result.factor_strict_conjunction_p = result.factor_strict_conjunction_q_by = nan;
    result.factor_min_log_bf = nan;
    result.factor_log_e_p2e_mediation = result.factor_e_q_p2e_ebh = nan;
    result.factor_pp_xm = result.factor_pp_my = result.factor_pp_two_stage = nan;
    result.factor_posterior_local_fdr = result.factor_posterior_cum_fdr = nan;
    result.factor_posterior_rank = -1;
    result.factor_pleiotropy_rho = result.factor_pleiotropy_p = nan;
    result.factor_nA = 0;
    result.factor_nB = prot.nB();
    result.factor_ld_source = prot.ld_reference_used ? "reference" : "identity";
    result.factor_selection_design = opts.factor_independent_selection
        ? "independent-discovery" : "same-sample";
    result.factor_effect_estimator =
        "joint-directional-generalized-adjusted-profile-score";
    result.factor_pattern = "NOT_RUN";
    result.factor_two_stage_status = "NOT_RUN";
    result.factor_mediation_status = "NOT_RUN";
    result.factor_frequentist_status = "NOT_RUN";
    result.factor_strict_status = "NOT_RUN";
    result.factor_p2e_status = "NOT_RUN";
    result.factor_ebh_status = "NOT_RUN";
    result.factor_posterior_status = "NOT_RUN";
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
        opts.factor_quadrature_points, opts.factor_pleio_sd_xm,
        opts.factor_pleio_quadrature_points,
        opts.factor_prior_xm, opts.factor_prior_directional,
        opts.factor_directional_variance);
    const FactorFit stage2 = fit_eiv(
        prot.setB_alpha_cis, prot.setB_se_alpha_cis,
        prot.setB_Gamma_cis, prot.setB_se_Gamma_cis, ld_b,
        opts.sampling_corr_pqtl_outcome, opts.prior_sigma2_beta2,
        opts.factor_quadrature_points, opts.factor_pleio_sd_my,
        opts.factor_pleio_quadrature_points,
        opts.factor_prior_my, opts.factor_prior_directional,
        opts.factor_directional_variance);
    if (!stage1.effect_valid || !stage2.effect_valid) {
        result.factor_effect_estimator =
            "unresolved-joint-directional-curvature";
    }

    if (stage1.valid) {
        result.factor_p_xm = stage1.p;
        result.factor_p_xm_strict = stage1.p_strict;
        result.factor_log_bf_xm = stage1.log_bf;
        result.factor_log_bf_heterogeneity_xm = stage1.log_bf_heterogeneity;
        result.factor_log_bf_directional_xm = stage1.log_bf_directional;
        result.factor_log_bf_slope_only_xm = stage1.log_bf_slope_only;
        result.factor_log_bf_directional_only_xm = stage1.log_bf_directional_only;
        result.factor_log_bf_slope_directional_xm = stage1.log_bf_slope_directional;
        result.factor_pp_directional_xm = stage1.pp_directional;
        result.factor_directional_collinearity_xm = stage1.directional_collinearity;
        result.factor_log_e_xm = stage1.log_e;
        result.factor_tau_xm = stage1.tau;
    }
    if (stage1.effect_valid) {
        result.factor_beta1 = stage1.beta;
        result.factor_beta1_se = stage1.se;
        result.factor_directional_intercept_xm = stage1.directional_intercept;
        result.factor_directional_intercept_xm_se = stage1.directional_intercept_se;
        const double critical = student_t_critical(0.05, result.factor_nA - 2.0);
        result.factor_beta1_ci_lower = stage1.beta - critical * stage1.se;
        result.factor_beta1_ci_upper = stage1.beta + critical * stage1.se;
    }
    if (stage2.valid) {
        result.factor_p_my = stage2.p;
        result.factor_p_my_strict = stage2.p_strict;
        result.factor_log_bf_my = stage2.log_bf;
        result.factor_log_bf_heterogeneity_my = stage2.log_bf_heterogeneity;
        result.factor_log_bf_directional_my = stage2.log_bf_directional;
        result.factor_log_bf_slope_only_my = stage2.log_bf_slope_only;
        result.factor_log_bf_directional_only_my = stage2.log_bf_directional_only;
        result.factor_log_bf_slope_directional_my = stage2.log_bf_slope_directional;
        result.factor_pp_directional_my = stage2.pp_directional;
        result.factor_directional_collinearity_my = stage2.directional_collinearity;
        result.factor_log_e_my = stage2.log_e;
        result.factor_tau_my = stage2.tau;
    }
    if (stage2.effect_valid) {
        result.factor_beta2 = stage2.beta;
        result.factor_beta2_se = stage2.se;
        result.factor_directional_intercept_my = stage2.directional_intercept;
        result.factor_directional_intercept_my_se = stage2.directional_intercept_se;
        const double critical = student_t_critical(0.05, result.factor_nB - 2.0);
        result.factor_beta2_ci_lower = stage2.beta - critical * stage2.se;
        result.factor_beta2_ci_upper = stage2.beta + critical * stage2.se;
    }
    if (stage1.valid && stage2.valid) {
        result.factor_pp_xm = stage1.pp_slope;
        result.factor_pp_my = stage2.pp_slope;
        result.factor_pp_two_stage = result.factor_pp_xm * result.factor_pp_my;
        result.factor_posterior_local_fdr = 1.0 - result.factor_pp_two_stage;
    }

    FactorFit direct;
    if (stage2.effect_valid && !gamma.empty()) {
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
                         opts.factor_quadrature_points, opts.factor_pleio_sd_xy,
                         opts.factor_pleio_quadrature_points, 0.5,
                         opts.factor_prior_directional,
                         opts.factor_directional_variance);
        if (direct.valid) {
            result.factor_p_xy = direct.p;
            result.factor_p_xy_strict = direct.p_strict;
            result.factor_log_bf_xy = direct.log_bf;
            result.factor_log_bf_heterogeneity_xy = direct.log_bf_heterogeneity;
            result.factor_log_e_xy = direct.log_e;
            result.factor_tau_xy = direct.tau;
        }
        if (direct.effect_valid) {
            result.factor_beta3 = direct.beta;
            result.factor_beta3_se = direct.se;
        }
    }

    if (stage1.effect_valid && stage2.effect_valid) {
        result.factor_indirect = stage1.beta * stage2.beta;
        result.factor_indirect_se = std::sqrt(
            stage2.beta * stage2.beta * stage1.se * stage1.se +
            stage1.beta * stage1.beta * stage2.se * stage2.se +
            stage1.se * stage1.se * stage2.se * stage2.se);
        const double critical_a = student_t_critical(
            0.025, result.factor_nA - 2.0);
        const double critical_b = student_t_critical(
            0.025, result.factor_nB - 2.0);
        const double a_lower = stage1.beta - critical_a * stage1.se;
        const double a_upper = stage1.beta + critical_a * stage1.se;
        const double b_lower = stage2.beta - critical_b * stage2.se;
        const double b_upper = stage2.beta + critical_b * stage2.se;
        const std::vector<double> products = {
            a_lower * b_lower, a_lower * b_upper,
            a_upper * b_lower, a_upper * b_upper};
        result.factor_indirect_ci_lower =
            *std::min_element(products.begin(), products.end());
        result.factor_indirect_ci_upper =
            *std::max_element(products.begin(), products.end());
    }
    if (stage1.valid && stage2.valid) {
        result.factor_conjunction_p = std::max(stage1.p, stage2.p);
        result.factor_strict_conjunction_p =
            std::max(stage1.p_strict, stage2.p_strict);
        result.factor_log_e_p2e_mediation =
            p_to_e_log_mixture(result.factor_strict_conjunction_p);
        result.factor_min_log_bf = std::min(stage1.log_bf, stage2.log_bf);
        result.factor_log_e_mediation = std::min(stage1.log_e, stage2.log_e);
    }

    // A de-noised residual-correlation diagnostic for coexisting pleiotropy.
    // It is deliberately not part of the confirmatory mediation p-value.
    if (stage1.effect_valid && stage2.effect_valid && direct.effect_valid &&
        gamma.size() >= 4) {
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

    auto evidence_gate = [&](bool evidence, const std::string& no_evidence,
                             const std::string& supported) {
        if (result.factor_nA < opts.factor_min_set_a && result.factor_nB < opts.factor_min_set_b)
            return std::string("INSUFFICIENT_SET_A_AND_SET_B");
        if (result.factor_nA < opts.factor_min_set_a) return std::string("INSUFFICIENT_SET_A");
        if (result.factor_nB < opts.factor_min_set_b) return std::string("INSUFFICIENT_SET_B");
        if (!stage1.valid || !stage2.valid) return std::string("NUMERICAL_FAILURE");
        if (std::fabs(opts.sampling_corr_rf_pqtl) > EPS ||
            std::fabs(opts.sampling_corr_pqtl_outcome) > EPS)
            return std::string("UNRESOLVED_SAMPLE_OVERLAP");
        if (!evidence) return no_evidence;
        return supported;
    };
    const bool bayes_evidence = std::isfinite(result.factor_min_log_bf) &&
        result.factor_min_log_bf >= std::log(opts.factor_bf_threshold);
    result.factor_two_stage_status = evidence_gate(
        bayes_evidence, "NO_TWO_STAGE_BAYES_EVIDENCE", "TWO_STAGE_EVIDENCE");
    result.factor_mediation_status = identification_gate(
        result, opts, result.factor_two_stage_status,
        "SUPPORTED_EXCLUSION_RESTRICTION_CONDITIONAL");
    const bool frequentist_evidence = std::isfinite(result.factor_conjunction_p) &&
        result.factor_conjunction_p <= opts.factor_alpha;
    result.factor_frequentist_status = evidence_gate(
        frequentist_evidence, "NO_TWO_STAGE_FREQUENTIST_EVIDENCE",
        "PENDING_MULTIPLE_TESTING");
    const bool strict_evidence =
        std::isfinite(result.factor_strict_conjunction_p) &&
        result.factor_strict_conjunction_p <= opts.factor_alpha;
    result.factor_strict_status = evidence_gate(
        strict_evidence, "NO_TWO_STAGE_STRICT_EVIDENCE",
        "PENDING_MULTIPLE_TESTING");
    if (result.factor_nA < opts.factor_min_set_a && result.factor_nB < opts.factor_min_set_b)
        result.factor_ebh_status = "INSUFFICIENT_SET_A_AND_SET_B";
    else if (result.factor_nA < opts.factor_min_set_a)
        result.factor_ebh_status = "INSUFFICIENT_SET_A";
    else if (result.factor_nB < opts.factor_min_set_b)
        result.factor_ebh_status = "INSUFFICIENT_SET_B";
    else if (!stage1.valid || !stage2.valid)
        result.factor_ebh_status = "NUMERICAL_FAILURE";
    else
        result.factor_ebh_status = "PENDING_MULTIPLE_TESTING";
    if (result.factor_nA < opts.factor_min_set_a && result.factor_nB < opts.factor_min_set_b)
        result.factor_p2e_status = "INSUFFICIENT_SET_A_AND_SET_B";
    else if (result.factor_nA < opts.factor_min_set_a)
        result.factor_p2e_status = "INSUFFICIENT_SET_A";
    else if (result.factor_nB < opts.factor_min_set_b)
        result.factor_p2e_status = "INSUFFICIENT_SET_B";
    else if (!stage1.valid || !stage2.valid)
        result.factor_p2e_status = "NUMERICAL_FAILURE";
    else
        result.factor_p2e_status = "PENDING_MULTIPLE_TESTING";
    if (result.factor_nA < opts.factor_min_set_a && result.factor_nB < opts.factor_min_set_b)
        result.factor_posterior_status = "INSUFFICIENT_SET_A_AND_SET_B";
    else if (result.factor_nA < opts.factor_min_set_a)
        result.factor_posterior_status = "INSUFFICIENT_SET_A";
    else if (result.factor_nB < opts.factor_min_set_b)
        result.factor_posterior_status = "INSUFFICIENT_SET_B";
    else if (!stage1.valid || !stage2.valid)
        result.factor_posterior_status = "NUMERICAL_FAILURE";
    else
        result.factor_posterior_status = "PENDING_MULTIPLE_TESTING";
}

void finalize_factorized_multiple_testing(std::vector<ProteinResult>& results,
                                          const Options& opts) {
    if (opts.structural_method != "factorized") return;
    const bool causal_leg_overlap =
        std::fabs(opts.sampling_corr_rf_pqtl) > EPS ||
        std::fabs(opts.sampling_corr_pqtl_outcome) > EPS;
    if (causal_leg_overlap) {
        // Selection occurs in the exposure GWAS for each leg. With correlated
        // estimation errors, the conditional p/e-value arguments used below
        // no longer apply. Preserve the exploratory leg statistics, but do
        // not expose confirmatory multiple-testing values.
        for (auto& result : results) {
            result.factor_conjunction_q_by =
                std::numeric_limits<double>::quiet_NaN();
            result.factor_strict_conjunction_q_by =
                std::numeric_limits<double>::quiet_NaN();
            result.factor_e_q_ebh = std::numeric_limits<double>::quiet_NaN();
            result.factor_e_q_p2e_ebh = std::numeric_limits<double>::quiet_NaN();
            if (result.factor_ebh_status == "PENDING_MULTIPLE_TESTING") {
                result.factor_ebh_status = "UNRESOLVED_SAMPLE_OVERLAP";
            }
            if (result.factor_p2e_status == "PENDING_MULTIPLE_TESTING") {
                result.factor_p2e_status = "UNRESOLVED_SAMPLE_OVERLAP";
            }
        }
        return;
    }
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
        if (result.factor_conjunction_q_by > opts.factor_alpha) {
            result.factor_frequentist_status = "NOT_SELECTED_BY_FDR";
            continue;
        }
        result.factor_frequentist_status = identification_gate(
            result, opts, "TWO_STAGE_EVIDENCE",
            "SUPPORTED_SCALAR_DISPERSION_EXCLUSION_CONDITIONAL");
    }

    std::vector<int> strict_order;
    for (int i = 0; i < static_cast<int>(results.size()); ++i) {
        if (std::isfinite(results[i].factor_strict_conjunction_p)) {
            strict_order.push_back(i);
        }
    }
    std::sort(strict_order.begin(), strict_order.end(), [&](int a, int b) {
        return results[a].factor_strict_conjunction_p <
               results[b].factor_strict_conjunction_p;
    });
    const int strict_m = static_cast<int>(strict_order.size());
    double strict_harmonic = 0.0;
    for (int i = 1; i <= strict_m; ++i) strict_harmonic += 1.0 / i;
    double strict_running = 1.0;
    for (int rank = strict_m; rank >= 1; --rank) {
        ProteinResult& result = results[strict_order[rank - 1]];
        const double raw = result.factor_strict_conjunction_p *
                           strict_m * strict_harmonic / rank;
        strict_running = std::min(strict_running, std::min(1.0, raw));
        result.factor_strict_conjunction_q_by = strict_running;
    }
    for (auto& result : results) {
        if (result.factor_strict_status != "PENDING_MULTIPLE_TESTING") continue;
        if (result.factor_strict_conjunction_q_by > opts.factor_alpha) {
            result.factor_strict_status = "NOT_SELECTED_BY_STRICT_BY";
            continue;
        }
        result.factor_strict_status = identification_gate(
            result, opts, "TWO_STAGE_EVIDENCE",
            "SUPPORTED_GAUSSIAN_COVARIANCE_EXCLUSION_CONDITIONAL");
    }

    std::vector<int> e_order;
    for (int i = 0; i < static_cast<int>(results.size()); ++i) {
        if (std::isfinite(results[i].factor_log_e_mediation)) e_order.push_back(i);
    }
    std::sort(e_order.begin(), e_order.end(), [&](int a, int b) {
        return results[a].factor_log_e_mediation > results[b].factor_log_e_mediation;
    });
    const int e_m = static_cast<int>(e_order.size());
    double running_log_q = std::numeric_limits<double>::infinity();
    for (int rank = e_m; rank >= 1; --rank) {
        ProteinResult& result = results[e_order[rank - 1]];
        const double log_raw_q = std::log(static_cast<double>(e_m)) -
            std::log(static_cast<double>(rank)) - result.factor_log_e_mediation;
        running_log_q = std::min(running_log_q, log_raw_q);
        result.factor_e_q_ebh = running_log_q >= 0.0 ? 1.0 : std::exp(running_log_q);
    }
    for (auto& result : results) {
        if (result.factor_ebh_status != "PENDING_MULTIPLE_TESTING") continue;
        if (!(std::isfinite(result.factor_e_q_ebh) &&
              result.factor_e_q_ebh <= opts.factor_alpha)) {
            result.factor_ebh_status = "NOT_SELECTED_BY_EBH";
            continue;
        }
        result.factor_ebh_status = identification_gate(
            result, opts, "TWO_STAGE_EVIDENCE",
            "SUPPORTED_EXCLUSION_RESTRICTION_CONDITIONAL");
    }

    std::vector<int> p2e_order;
    for (int i = 0; i < static_cast<int>(results.size()); ++i) {
        if (std::isfinite(results[i].factor_log_e_p2e_mediation)) {
            p2e_order.push_back(i);
        }
    }
    std::sort(p2e_order.begin(), p2e_order.end(), [&](int a, int b) {
        return results[a].factor_log_e_p2e_mediation >
               results[b].factor_log_e_p2e_mediation;
    });
    const int p2e_m = static_cast<int>(p2e_order.size());
    double p2e_running_log_q = std::numeric_limits<double>::infinity();
    for (int rank = p2e_m; rank >= 1; --rank) {
        ProteinResult& result = results[p2e_order[rank - 1]];
        const double log_raw_q = std::log(static_cast<double>(p2e_m)) -
            std::log(static_cast<double>(rank)) -
            result.factor_log_e_p2e_mediation;
        p2e_running_log_q = std::min(p2e_running_log_q, log_raw_q);
        result.factor_e_q_p2e_ebh = p2e_running_log_q >= 0.0
            ? 1.0 : std::exp(p2e_running_log_q);
    }
    for (auto& result : results) {
        if (result.factor_p2e_status != "PENDING_MULTIPLE_TESTING") continue;
        if (!(std::isfinite(result.factor_e_q_p2e_ebh) &&
              result.factor_e_q_p2e_ebh <= opts.factor_alpha)) {
            result.factor_p2e_status = "NOT_SELECTED_BY_P2E_EBH";
            continue;
        }
        result.factor_p2e_status = identification_gate(
            result, opts, "TWO_STAGE_EVIDENCE",
            "SUPPORTED_GAUSSIAN_COVARIANCE_EXCLUSION_CONDITIONAL");
    }

    std::vector<int> posterior_order;
    for (int i = 0; i < static_cast<int>(results.size()); ++i) {
        if (std::isfinite(results[i].factor_posterior_local_fdr)) {
            posterior_order.push_back(i);
        }
    }
    std::sort(posterior_order.begin(), posterior_order.end(), [&](int a, int b) {
        return results[a].factor_posterior_local_fdr <
               results[b].factor_posterior_local_fdr;
    });
    double posterior_fdr_sum = 0.0;
    for (int rank = 1; rank <= static_cast<int>(posterior_order.size()); ++rank) {
        ProteinResult& result = results[posterior_order[rank - 1]];
        posterior_fdr_sum += result.factor_posterior_local_fdr;
        result.factor_posterior_rank = rank;
        result.factor_posterior_cum_fdr = posterior_fdr_sum / rank;
    }
    for (auto& result : results) {
        if (result.factor_posterior_status != "PENDING_MULTIPLE_TESTING") continue;
        if (!(std::isfinite(result.factor_posterior_cum_fdr) &&
              result.factor_posterior_cum_fdr <= opts.factor_alpha)) {
            result.factor_posterior_status = "NOT_SELECTED_BY_POSTERIOR_FDR";
            continue;
        }
        result.factor_posterior_status = identification_gate(
            result, opts, "TWO_STAGE_EVIDENCE",
            "SUPPORTED_BAYES_POSTERIOR_FDR_ASSUMPTION_CONDITIONAL");
    }
}

} // namespace bmediator
