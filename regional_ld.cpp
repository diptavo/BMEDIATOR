#include "regional_ld.h"

#include <array>
#include <unordered_set>

namespace bmediator {

namespace {

struct FineMapSignal {
    int lead = -1;
    std::vector<double> log_bf;
    std::vector<int> credible_set;
};

struct ConditionalScores {
    std::vector<double> z;
    std::vector<double> residual_variance;
};

using Matrix = std::vector<std::vector<double>>;

struct JointModelFit {
    double log_evidence = -std::numeric_limits<double>::infinity();
    std::array<double, 3> beta{{0.0, 0.0, 0.0}};
    std::array<double, 3> beta_se{{0.0, 0.0, 0.0}};
};

double log_sum_exp(const std::vector<double>& values) {
    if (values.empty()) return -std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (double value : values) {
        if (std::isfinite(value)) maximum = std::max(maximum, value);
    }
    if (!std::isfinite(maximum)) return maximum;
    double total = 0.0;
    for (double value : values) {
        if (std::isfinite(value)) total += std::exp(value - maximum);
    }
    return maximum + std::log(total);
}

double wakefield_log_bf_from_z(double z, double se, double prior_variance) {
    if (!std::isfinite(z) || !std::isfinite(se) || se <= 0.0 ||
        !std::isfinite(prior_variance) || prior_variance <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    const double variance = se * se;
    const double shrinkage = prior_variance / (variance + prior_variance);
    return 0.5 * (std::log1p(-shrinkage) + shrinkage * z * z);
}

bool invert_matrix(std::vector<std::vector<double>> matrix,
                   std::vector<std::vector<double>>& inverse) {
    const int n = static_cast<int>(matrix.size());
    inverse.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        matrix[i][i] += 1e-8;
        inverse[i][i] = 1.0;
    }
    for (int column = 0; column < n; ++column) {
        int pivot = column;
        for (int row = column + 1; row < n; ++row) {
            if (std::fabs(matrix[row][column]) > std::fabs(matrix[pivot][column])) {
                pivot = row;
            }
        }
        if (std::fabs(matrix[pivot][column]) < 1e-10) return false;
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            std::swap(inverse[pivot], inverse[column]);
        }
        const double scale = matrix[column][column];
        for (int j = 0; j < n; ++j) {
            matrix[column][j] /= scale;
            inverse[column][j] /= scale;
        }
        for (int row = 0; row < n; ++row) {
            if (row == column) continue;
            const double factor = matrix[row][column];
            for (int j = 0; j < n; ++j) {
                matrix[row][j] -= factor * matrix[column][j];
                inverse[row][j] -= factor * inverse[column][j];
            }
        }
    }
    return true;
}

bool cholesky_decompose(const Matrix& matrix, Matrix& lower) {
    const int n = static_cast<int>(matrix.size());
    lower.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double value = matrix[i][j];
            for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
            if (i == j) {
                if (!(value > 1e-12) || !std::isfinite(value)) return false;
                lower[i][j] = std::sqrt(value);
            } else {
                lower[i][j] = value / lower[j][j];
            }
        }
    }
    return true;
}

std::vector<double> solve_cholesky(const Matrix& lower,
                                   const std::vector<double>& rhs) {
    const int n = static_cast<int>(lower.size());
    std::vector<double> intermediate(n, 0.0), solution(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double value = rhs[i];
        for (int j = 0; j < i; ++j) value -= lower[i][j] * intermediate[j];
        intermediate[i] = value / lower[i][i];
    }
    for (int i = n - 1; i >= 0; --i) {
        double value = intermediate[i];
        for (int j = i + 1; j < n; ++j) value -= lower[j][i] * solution[j];
        solution[i] = value / lower[i][i];
    }
    return solution;
}

double log_determinant_cholesky(const Matrix& lower) {
    double value = 0.0;
    for (int i = 0; i < static_cast<int>(lower.size()); ++i) {
        value += 2.0 * std::log(lower[i][i]);
    }
    return value;
}

double log_mvn_zero(const std::vector<double>& observed, Matrix covariance) {
    Matrix lower;
    double jitter = 1e-10;
    for (int attempt = 0; attempt < 8; ++attempt) {
        if (cholesky_decompose(covariance, lower)) {
            const std::vector<double> solved = solve_cholesky(lower, observed);
            double quadratic = 0.0;
            for (size_t i = 0; i < observed.size(); ++i) {
                quadratic += observed[i] * solved[i];
            }
            return -0.5 * (static_cast<double>(observed.size()) * LOG2PI +
                           log_determinant_cholesky(lower) + quadratic);
        }
        for (int i = 0; i < static_cast<int>(covariance.size()); ++i) {
            covariance[i][i] += jitter;
        }
        jitter *= 10.0;
    }
    return -std::numeric_limits<double>::infinity();
}

class RegionalLD {
public:
    RegionalLD(const PlinkData& plink, const std::vector<int>& indices) {
        values_ = plink.extract_genotypes_double(indices);
        n_samples_ = plink.n_samples;
        for (auto& variant : values_) {
            double mean = std::accumulate(variant.begin(), variant.end(), 0.0) /
                          static_cast<double>(n_samples_);
            double sum_squares = 0.0;
            for (double value : variant) {
                const double centered = value - mean;
                sum_squares += centered * centered;
            }
            const double sd = std::sqrt(sum_squares / static_cast<double>(n_samples_));
            if (sd <= 1e-10) {
                std::fill(variant.begin(), variant.end(), 0.0);
            } else {
                for (double& value : variant) value = (value - mean) / sd;
            }
        }
    }

    void cache_rows(const std::vector<int>& indices) const {
        for (int index : indices) {
            if (row_cache_.count(index)) continue;
            std::vector<double> row(values_.size(), 0.0);
            for (int other = 0; other < static_cast<int>(values_.size()); ++other) {
                row[other] = direct_correlation(index, other);
            }
            row_cache_[index] = std::move(row);
        }
    }

    double correlation(int left, int right) const {
        auto left_row = row_cache_.find(left);
        if (left_row != row_cache_.end()) return left_row->second[right];
        auto right_row = row_cache_.find(right);
        if (right_row != row_cache_.end()) return right_row->second[left];
        return direct_correlation(left, right);
    }

private:
    double direct_correlation(int left, int right) const {
        if (left == right) return 1.0;
        double value = 0.0;
        for (int sample = 0; sample < n_samples_; ++sample) {
            value += values_[left][sample] * values_[right][sample];
        }
        value /= static_cast<double>(n_samples_);
        return std::max(-1.0, std::min(1.0, value));
    }

    std::vector<std::vector<double>> values_;
    int n_samples_ = 0;
    mutable std::map<int, std::vector<double>> row_cache_;
};

ConditionalScores conditional_z_scores(const std::vector<double>& z,
                                       const std::vector<int>& conditioning,
                                       const RegionalLD& ld) {
    if (conditioning.empty()) return {z, std::vector<double>(z.size(), 1.0)};
    ld.cache_rows(conditioning);
    const int k = static_cast<int>(conditioning.size());
    std::vector<std::vector<double>> r_ss(k, std::vector<double>(k));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            r_ss[i][j] = ld.correlation(conditioning[i], conditioning[j]);
        }
    }
    std::vector<std::vector<double>> inverse;
    ConditionalScores result{
        std::vector<double>(z.size(), std::numeric_limits<double>::quiet_NaN()),
        std::vector<double>(z.size(), std::numeric_limits<double>::quiet_NaN())
    };
    if (!invert_matrix(r_ss, inverse)) return result;

    std::vector<double> inverse_z(k, 0.0);
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            inverse_z[i] += inverse[i][j] * z[conditioning[j]];
        }
    }
    std::unordered_set<int> conditioned(conditioning.begin(), conditioning.end());
    for (int variant = 0; variant < static_cast<int>(z.size()); ++variant) {
        if (conditioned.count(variant)) continue;
        std::vector<double> r(k, 0.0);
        for (int i = 0; i < k; ++i) {
            r[i] = ld.correlation(variant, conditioning[i]);
        }
        double conditional_mean = 0.0;
        double explained_variance = 0.0;
        for (int i = 0; i < k; ++i) {
            conditional_mean += r[i] * inverse_z[i];
            for (int j = 0; j < k; ++j) {
                explained_variance += r[i] * inverse[i][j] * r[j];
            }
        }
        const double residual_variance = 1.0 - explained_variance;
        if (residual_variance <= 1e-4) continue;
        result.residual_variance[variant] = residual_variance;
        result.z[variant] = (z[variant] - conditional_mean) /
                            std::sqrt(residual_variance);
    }
    return result;
}

std::vector<int> select_independent_signals(const std::vector<double>& z,
                                            const RegionalLD& ld,
                                            const Options& opts) {
    std::vector<int> selected;
    while (static_cast<int>(selected.size()) < opts.regional_max_signals) {
        ConditionalScores conditional = conditional_z_scores(z, selected, ld);
        int best = -1;
        double best_abs_z = 0.0;
        for (int variant = 0; variant < static_cast<int>(conditional.z.size()); ++variant) {
            if (!std::isfinite(conditional.z[variant])) continue;
            const double abs_z = std::fabs(conditional.z[variant]);
            if (abs_z > best_abs_z) {
                best_abs_z = abs_z;
                best = variant;
            }
        }
        if (best < 0) break;
        const double p = std::erfc(best_abs_z / std::sqrt(2.0));
        if (!std::isfinite(p) || p > opts.regional_signal_p) break;
        selected.push_back(best);
    }
    return selected;
}

std::vector<FineMapSignal> fine_map_signals(const std::vector<double>& beta,
                                            const std::vector<double>& se,
                                            const RegionalLD& ld,
                                            double prior_variance,
                                            const Options& opts) {
    std::vector<double> z(beta.size());
    for (size_t i = 0; i < beta.size(); ++i) z[i] = beta[i] / se[i];
    const std::vector<int> selected = select_independent_signals(z, ld, opts);
    std::vector<FineMapSignal> signals;
    for (size_t signal_index = 0; signal_index < selected.size(); ++signal_index) {
        std::vector<int> other_signals;
        for (size_t other = 0; other < selected.size(); ++other) {
            if (other != signal_index) other_signals.push_back(selected[other]);
        }
        const ConditionalScores conditional = conditional_z_scores(z, other_signals, ld);
        FineMapSignal signal;
        signal.lead = selected[signal_index];
        signal.log_bf.resize(beta.size());
        for (size_t variant = 0; variant < beta.size(); ++variant) {
            signal.log_bf[variant] = wakefield_log_bf_from_z(
                conditional.z[variant],
                se[variant] / std::sqrt(conditional.residual_variance[variant]),
                prior_variance
            );
        }

        const double normalizer = log_sum_exp(signal.log_bf);
        std::vector<std::pair<double, int>> posterior;
        posterior.reserve(beta.size());
        for (int variant = 0; variant < static_cast<int>(beta.size()); ++variant) {
            if (!std::isfinite(signal.log_bf[variant])) continue;
            posterior.push_back({std::exp(signal.log_bf[variant] - normalizer), variant});
        }
        std::sort(posterior.begin(), posterior.end(),
                  [](const auto& left, const auto& right) { return left.first > right.first; });
        double coverage = 0.0;
        for (const auto& item : posterior) {
            signal.credible_set.push_back(item.second);
            coverage += item.first;
            if (coverage >= opts.regional_coverage) break;
        }
        signals.push_back(std::move(signal));
    }
    return signals;
}

Matrix transformed_summary_covariance(const std::vector<double>& left_se,
                                      const std::vector<double>& right_se,
                                      const Matrix& ld,
                                      const Matrix& inverse_ld,
                                      double correlation) {
    const int n = static_cast<int>(ld.size());
    Matrix covariance(n, std::vector<double>(n, 0.0));
    if (correlation == 0.0) return covariance;
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < n; ++column) {
            double value = 0.0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    value += inverse_ld[row][i] * left_se[i] * ld[i][j] *
                             right_se[j] * inverse_ld[column][j];
                }
            }
            covariance[row][column] = correlation * value;
        }
    }
    return covariance;
}

std::vector<double> transform_marginal_effects(const std::vector<double>& effects,
                                               const Matrix& inverse_ld) {
    const int n = static_cast<int>(effects.size());
    std::vector<double> transformed(n, 0.0);
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < n; ++column) {
            transformed[row] += inverse_ld[row][column] * effects[column];
        }
    }
    return transformed;
}

std::array<bool, 3> active_betas(int scenario) {
    return {{scenario == 1 || scenario == 2 || scenario == 5,
             scenario == 1 || scenario == 4,
             scenario == 1 || scenario == 3 || scenario == 5}};
}

std::array<double, 3> unpack_betas(int scenario,
                                   const std::vector<double>& parameters) {
    std::array<double, 3> beta{{0.0, 0.0, 0.0}};
    const auto active = active_betas(scenario);
    int position = 0;
    for (int index = 0; index < 3; ++index) {
        if (active[index]) beta[index] = parameters[position++];
    }
    return beta;
}

Matrix joint_structural_covariance(const Matrix& measurement_covariance,
                                   int components,
                                   int scenario,
                                   const std::vector<double>& parameters,
                                   const Options& opts) {
    Matrix covariance = measurement_covariance;
    const auto beta = unpack_betas(scenario, parameters);
    const double beta1 = beta[0], beta2 = beta[1], beta3 = beta[2];
    const Matrix mapping = {
        {1.0, 0.0, 0.0},
        {beta1, 1.0, 0.0},
        {beta3 + beta1 * beta2, beta2, 1.0}
    };
    Matrix latent = {
        {opts.regional_prior_var_rf, 0.0, 0.0},
        {0.0, opts.regional_prior_var_pp, 0.0},
        {0.0, 0.0, opts.regional_prior_var_outcome}
    };
    if (scenario == 5) {
        const double covariance_dh = opts.regional_pleiotropy_rho *
            std::sqrt(opts.regional_prior_var_pp * opts.regional_prior_var_outcome);
        latent[1][2] = covariance_dh;
        latent[2][1] = covariance_dh;
    }

    Matrix trait_covariance(3, std::vector<double>(3, 0.0));
    for (int left = 0; left < 3; ++left) {
        for (int right = 0; right < 3; ++right) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    trait_covariance[left][right] +=
                        mapping[left][i] * latent[i][j] * mapping[right][j];
                }
            }
        }
    }
    for (int component = 0; component < components; ++component) {
        for (int left = 0; left < 3; ++left) {
            for (int right = 0; right < 3; ++right) {
                covariance[left * components + component]
                          [right * components + component] += trait_covariance[left][right];
            }
        }
    }
    return covariance;
}

double joint_log_posterior(const std::vector<double>& observed,
                           const Matrix& measurement_covariance,
                           int components,
                           int scenario,
                           const std::vector<double>& parameters,
                           const Options& opts) {
    Matrix covariance = joint_structural_covariance(
        measurement_covariance, components, scenario, parameters, opts
    );
    double value = log_mvn_zero(observed, std::move(covariance));
    if (!std::isfinite(value)) return value;
    const auto active = active_betas(scenario);
    const std::array<double, 3> prior_variance{{
        opts.prior_sigma2_beta1,
        opts.prior_sigma2_beta2,
        opts.prior_sigma2_beta3
    }};
    int position = 0;
    for (int index = 0; index < 3; ++index) {
        if (!active[index]) continue;
        const double variance = prior_variance[index];
        const double parameter = parameters[position++];
        value += -0.5 * (LOG2PI + std::log(variance) + parameter * parameter / variance);
    }
    return value;
}

std::array<double, 3> initial_joint_betas(const std::vector<double>& rf,
                                          const std::vector<double>& protein,
                                          const std::vector<double>& outcome) {
    std::array<double, 3> result{{0.0, 0.0, 0.0}};
    double xx = 0.0, mm = 0.0, xm = 0.0, xy = 0.0, my = 0.0;
    for (size_t i = 0; i < rf.size(); ++i) {
        xx += rf[i] * rf[i];
        mm += protein[i] * protein[i];
        xm += rf[i] * protein[i];
        xy += rf[i] * outcome[i];
        my += protein[i] * outcome[i];
    }
    if (xx > 1e-12) result[0] = xm / xx;
    const double determinant = xx * mm - xm * xm;
    if (determinant > 1e-12) {
        result[1] = (xx * my - xm * xy) / determinant;
        result[2] = (mm * xy - xm * my) / determinant;
    } else {
        if (mm > 1e-12) result[1] = my / mm;
        if (xx > 1e-12) result[2] = xy / xx;
    }
    for (double& value : result) value = std::max(-3.0, std::min(3.0, value));
    return result;
}

std::vector<double> pack_active_betas(int scenario,
                                      const std::array<double, 3>& beta) {
    const auto active = active_betas(scenario);
    std::vector<double> result;
    for (int index = 0; index < 3; ++index) {
        if (active[index]) result.push_back(beta[index]);
    }
    return result;
}

JointModelFit fit_joint_scenario(const std::vector<double>& observed,
                                 const Matrix& measurement_covariance,
                                 int components,
                                 int scenario,
                                 const std::array<double, 3>& initial,
                                 const Options& opts) {
    JointModelFit fit;
    std::vector<double> parameters = pack_active_betas(scenario, initial);
    if (parameters.empty()) {
        fit.log_evidence = joint_log_posterior(
            observed, measurement_covariance, components, scenario, parameters, opts
        );
        return fit;
    }

    auto objective = [&](const std::vector<double>& values) {
        return joint_log_posterior(
            observed, measurement_covariance, components, scenario, values, opts
        );
    };
    std::vector<std::vector<double>> starts;
    starts.push_back(parameters);
    starts.push_back(std::vector<double>(parameters.size(), 0.0));
    std::vector<double> reversed = parameters;
    for (double& value : reversed) value = -value;
    starts.push_back(reversed);

    double best_value = -std::numeric_limits<double>::infinity();
    std::vector<double> best_parameters = parameters;
    for (auto current : starts) {
        double current_value = objective(current);
        double step = 0.5;
        for (int iteration = 0; iteration < 160 && step >= 1e-4; ++iteration) {
            bool improved = false;
            for (int coordinate = 0; coordinate < static_cast<int>(current.size()); ++coordinate) {
                for (double direction : {-1.0, 1.0}) {
                    std::vector<double> candidate = current;
                    candidate[coordinate] = std::max(
                        -5.0, std::min(5.0, candidate[coordinate] + direction * step)
                    );
                    const double candidate_value = objective(candidate);
                    if (candidate_value > current_value) {
                        current = std::move(candidate);
                        current_value = candidate_value;
                        improved = true;
                    }
                }
            }
            if (!improved) step *= 0.5;
        }
        if (current_value > best_value) {
            best_value = current_value;
            best_parameters = std::move(current);
        }
    }

    const int dimension = static_cast<int>(best_parameters.size());
    Matrix precision(dimension, std::vector<double>(dimension, 0.0));
    std::vector<double> increments(dimension);
    for (int i = 0; i < dimension; ++i) {
        increments[i] = 1e-3 * (1.0 + std::fabs(best_parameters[i]));
        std::vector<double> plus = best_parameters, minus = best_parameters;
        plus[i] += increments[i];
        minus[i] -= increments[i];
        precision[i][i] = -(objective(plus) - 2.0 * best_value + objective(minus)) /
                           (increments[i] * increments[i]);
        precision[i][i] = std::max(precision[i][i], 1e-6);
    }
    for (int i = 0; i < dimension; ++i) {
        for (int j = i + 1; j < dimension; ++j) {
            std::vector<double> pp = best_parameters, pm = best_parameters;
            std::vector<double> mp = best_parameters, mm = best_parameters;
            pp[i] += increments[i]; pp[j] += increments[j];
            pm[i] += increments[i]; pm[j] -= increments[j];
            mp[i] -= increments[i]; mp[j] += increments[j];
            mm[i] -= increments[i]; mm[j] -= increments[j];
            const double hessian = (objective(pp) - objective(pm) - objective(mp) + objective(mm)) /
                (4.0 * increments[i] * increments[j]);
            precision[i][j] = -hessian;
            precision[j][i] = -hessian;
        }
    }

    Matrix lower;
    double jitter = 1e-8;
    bool precision_pd = false;
    for (int attempt = 0; attempt < 8; ++attempt) {
        if (cholesky_decompose(precision, lower)) {
            precision_pd = true;
            break;
        }
        for (int i = 0; i < dimension; ++i) precision[i][i] += jitter;
        jitter *= 10.0;
    }
    Matrix posterior_covariance;
    if (!precision_pd || !invert_matrix(precision, posterior_covariance)) {
        posterior_covariance.assign(dimension, std::vector<double>(dimension, 0.0));
        for (int i = 0; i < dimension; ++i) posterior_covariance[i][i] = 1.0 / precision[i][i];
        double log_det = 0.0;
        for (int i = 0; i < dimension; ++i) log_det += std::log(precision[i][i]);
        fit.log_evidence = best_value + 0.5 * dimension * LOG2PI - 0.5 * log_det;
    } else {
        fit.log_evidence = best_value + 0.5 * dimension * LOG2PI -
                           0.5 * log_determinant_cholesky(lower);
    }
    fit.beta = unpack_betas(scenario, best_parameters);
    const auto active = active_betas(scenario);
    int position = 0;
    for (int index = 0; index < 3; ++index) {
        if (active[index]) {
            fit.beta_se[index] = std::sqrt(std::max(
                posterior_covariance[position][position], 0.0
            ));
            position++;
        }
    }
    return fit;
}

double infinity_condition_number(const Matrix& matrix, const Matrix& inverse) {
    double matrix_norm = 0.0, inverse_norm = 0.0;
    for (size_t row = 0; row < matrix.size(); ++row) {
        double matrix_sum = 0.0, inverse_sum = 0.0;
        for (size_t column = 0; column < matrix.size(); ++column) {
            matrix_sum += std::fabs(matrix[row][column]);
            inverse_sum += std::fabs(inverse[row][column]);
        }
        matrix_norm = std::max(matrix_norm, matrix_sum);
        inverse_norm = std::max(inverse_norm, inverse_sum);
    }
    return matrix_norm * inverse_norm;
}

void compute_joint_regional_model(ProteinData& protein,
                                  const std::vector<FineMapSignal>& rf_signals,
                                  const std::vector<FineMapSignal>& protein_signals,
                                  const std::vector<FineMapSignal>& outcome_signals,
                                  const RegionalLD& ld,
                                  const Options& opts) {
    protein.regional_joint_evaluated = true;
    std::set<int> component_set;
    for (const auto& signal : rf_signals) component_set.insert(signal.lead);
    for (const auto& signal : protein_signals) component_set.insert(signal.lead);
    for (const auto& signal : outcome_signals) component_set.insert(signal.lead);
    const std::vector<int> components(component_set.begin(), component_set.end());
    protein.regional_joint_components = static_cast<int>(components.size());
    if (components.empty()) {
        protein.regional_joint_status = "NO_REGIONAL_SIGNAL";
        return;
    }

    const int n = static_cast<int>(components.size());
    Matrix component_ld(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double raw = ld.correlation(components[i], components[j]);
            component_ld[i][j] = (1.0 - opts.regional_ld_shrinkage) * raw;
            if (i == j) component_ld[i][j] += opts.regional_ld_shrinkage;
        }
    }
    Matrix inverse_ld;
    if (!invert_matrix(component_ld, inverse_ld)) {
        protein.regional_joint_status = "SINGULAR_COMPONENT_LD";
        return;
    }
    protein.regional_joint_condition_number =
        infinity_condition_number(component_ld, inverse_ld);

    std::array<std::vector<double>, 3> marginal_effects;
    std::array<std::vector<double>, 3> standard_errors;
    for (int component : components) {
        marginal_effects[0].push_back(protein.regional_joint_rf_beta[component]);
        marginal_effects[1].push_back(protein.regional_joint_pp_beta[component]);
        marginal_effects[2].push_back(protein.regional_joint_outcome_beta[component]);
        standard_errors[0].push_back(protein.regional_joint_rf_se[component]);
        standard_errors[1].push_back(protein.regional_joint_pp_se[component]);
        standard_errors[2].push_back(protein.regional_joint_outcome_se[component]);
    }
    std::array<std::vector<double>, 3> joint_effects;
    for (int trait = 0; trait < 3; ++trait) {
        joint_effects[trait] = transform_marginal_effects(marginal_effects[trait], inverse_ld);
    }
    std::vector<double> observed;
    for (int trait = 0; trait < 3; ++trait) {
        observed.insert(observed.end(), joint_effects[trait].begin(), joint_effects[trait].end());
    }

    const std::array<std::array<double, 3>, 3> overlap{{
        {{1.0, opts.overlap_rf_protein, opts.overlap_rf_outcome}},
        {{opts.overlap_rf_protein, 1.0, opts.overlap_protein_outcome}},
        {{opts.overlap_rf_outcome, opts.overlap_protein_outcome, 1.0}}
    }};
    Matrix measurement(3 * n, std::vector<double>(3 * n, 0.0));
    for (int left = 0; left < 3; ++left) {
        for (int right = 0; right < 3; ++right) {
            Matrix block = transformed_summary_covariance(
                standard_errors[left], standard_errors[right],
                component_ld, inverse_ld, overlap[left][right]
            );
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    measurement[left * n + i][right * n + j] = block[i][j];
                }
            }
        }
    }

    const auto initial = initial_joint_betas(
        joint_effects[0], joint_effects[1], joint_effects[2]
    );
    std::array<JointModelFit, 6> fits;
    for (int scenario = 0; scenario < 6; ++scenario) {
        fits[scenario] = fit_joint_scenario(
            observed, measurement, n, scenario, initial, opts
        );
        if (!std::isfinite(fits[scenario].log_evidence)) {
            protein.regional_joint_status = "NUMERICAL_FAILURE";
            return;
        }
    }
    const double baseline = fits[0].log_evidence;
    for (int scenario = 0; scenario < 6; ++scenario) {
        protein.regional_joint_log_bf[scenario] = fits[scenario].log_evidence - baseline;
        for (int index = 0; index < 3; ++index) {
            protein.regional_joint_scenario_beta[scenario * 3 + index] =
                fits[scenario].beta[index];
            protein.regional_joint_scenario_beta_se[scenario * 3 + index] =
                fits[scenario].beta_se[index];
        }
    }
    for (int index = 0; index < 3; ++index) {
        protein.regional_joint_beta[index] = fits[1].beta[index];
        protein.regional_joint_beta_se[index] = fits[1].beta_se[index];
    }
    protein.regional_joint_status = n < 2
        ? "UNRESOLVED_SINGLE_COMPONENT" : "EVALUATED";
}

RegionalSignalPairResult coloc_signal_pair(const FineMapSignal& protein_signal,
                                           const FineMapSignal& outcome_signal,
                                           int protein_index,
                                           int outcome_index,
                                           const ProteinData& protein,
                                           const RegionalLD& ld,
                                           const Options& opts) {
    RegionalSignalPairResult result;
    result.protein_signal = protein_index + 1;
    result.outcome_signal = outcome_index + 1;
    result.protein_lead = protein.regional_cis_rsid[protein_signal.lead];
    result.outcome_lead = protein.regional_cis_rsid[outcome_signal.lead];

    const double log_sum_protein = log_sum_exp(protein_signal.log_bf);
    const double log_sum_outcome = log_sum_exp(outcome_signal.log_bf);
    std::vector<double> shared_terms;
    shared_terms.reserve(protein_signal.log_bf.size());
    for (size_t i = 0; i < protein_signal.log_bf.size(); ++i) {
        shared_terms.push_back(protein_signal.log_bf[i] + outcome_signal.log_bf[i]);
    }
    const double log_sum_shared = log_sum_exp(shared_terms);
    const double log_all_pairs = log_sum_protein + log_sum_outcome;
    const double shared_fraction = std::exp(std::min(0.0, log_sum_shared - log_all_pairs));
    const double log_sum_distinct = shared_fraction >= 1.0
        ? -std::numeric_limits<double>::infinity()
        : log_all_pairs + std::log1p(-shared_fraction);

    const std::vector<double> hypotheses = {
        0.0,
        std::log(opts.regional_prior_pp) + log_sum_protein,
        std::log(opts.regional_prior_outcome) + log_sum_outcome,
        std::log(opts.regional_prior_pp) + std::log(opts.regional_prior_outcome) +
            log_sum_distinct,
        std::log(opts.regional_prior_shared) + log_sum_shared
    };
    const double normalizer = log_sum_exp(hypotheses);
    result.pp_h0 = std::exp(hypotheses[0] - normalizer);
    result.pp_h1 = std::exp(hypotheses[1] - normalizer);
    result.pp_h2 = std::exp(hypotheses[2] - normalizer);
    result.pp_h3 = std::exp(hypotheses[3] - normalizer);
    result.pp_h4 = std::exp(hypotheses[4] - normalizer);
    const double both = result.pp_h3 + result.pp_h4;
    result.shared_given_both = both > 0.0 ? result.pp_h4 / both : 0.0;
    const double lead_r = ld.correlation(protein_signal.lead, outcome_signal.lead);
    result.lead_pair_r2 = lead_r * lead_r;
    for (int left : protein_signal.credible_set) {
        for (int right : outcome_signal.credible_set) {
            const double r = ld.correlation(left, right);
            result.max_credible_set_pair_r2 = std::max(
                result.max_credible_set_pair_r2, r * r
            );
        }
    }

    if (both >= opts.regional_min_both &&
        result.shared_given_both >= opts.regional_min_shared) {
        result.interpretation = "SHARED_SIGNAL_SUPPORTED";
    } else if (both >= opts.regional_min_both &&
               result.shared_given_both <= 1.0 - opts.regional_min_shared) {
        result.interpretation = result.max_credible_set_pair_r2 >= opts.regional_high_ld_r2
            ? "DISTINCT_SIGNALS_HIGH_LD"
            : "DISTINCT_SIGNALS_LOW_MODERATE_LD";
    } else {
        result.interpretation = "SIGNAL_PAIR_AMBIGUOUS";
    }
    return result;
}

} // namespace

void compute_multisignal_regional_evidence(ProteinData& protein,
                                           const PlinkData& plink,
                                           const Options& opts) {
    protein.regional_method = opts.regional_method;
    protein.regional_multisignal_evaluated = true;
    protein.regional_signal_pairs.clear();
    const size_t variant_count = protein.regional_bim_index.size();
    if (!protein.regional_data_complete || variant_count < 2 ||
        protein.regional_pp_beta.size() != variant_count ||
        protein.regional_pp_se.size() != variant_count ||
        protein.regional_outcome_beta.size() != variant_count ||
        protein.regional_outcome_se.size() != variant_count) {
        protein.regional_multisignal_interpretation = "NO_REGIONAL_DATA";
        if (opts.regional_method == "joint-ld") {
            protein.regional_joint_evaluated = true;
            protein.regional_joint_status = "NO_REGIONAL_DATA";
        }
        return;
    }

    RegionalLD ld(plink, protein.regional_bim_index);
    const auto protein_signals = fine_map_signals(
        protein.regional_pp_beta, protein.regional_pp_se, ld,
        opts.regional_prior_var_pp, opts
    );
    const auto outcome_signals = fine_map_signals(
        protein.regional_outcome_beta, protein.regional_outcome_se, ld,
        opts.regional_prior_var_outcome, opts
    );
    protein.regional_protein_signals = static_cast<int>(protein_signals.size());
    protein.regional_outcome_signals = static_cast<int>(outcome_signals.size());

    if (opts.regional_method == "joint-ld") {
        protein.regional_joint_evaluated = true;
        const size_t joint_count = protein.regional_joint_bim_index.size();
        protein.regional_joint_n_variants = static_cast<int>(joint_count);
        const bool joint_complete = joint_count >= 2 &&
            protein.regional_joint_rsid.size() == joint_count &&
            protein.regional_joint_rf_beta.size() == joint_count &&
            protein.regional_joint_rf_se.size() == joint_count &&
            protein.regional_joint_pp_beta.size() == joint_count &&
            protein.regional_joint_pp_se.size() == joint_count &&
            protein.regional_joint_outcome_beta.size() == joint_count &&
            protein.regional_joint_outcome_se.size() == joint_count;
        if (!joint_complete) {
            protein.regional_joint_status = "NO_JOINT_REGIONAL_DATA";
        } else {
            RegionalLD joint_ld(plink, protein.regional_joint_bim_index);
            const auto joint_rf_signals = fine_map_signals(
                protein.regional_joint_rf_beta, protein.regional_joint_rf_se,
                joint_ld, opts.regional_prior_var_rf, opts
            );
            const auto joint_protein_signals = fine_map_signals(
                protein.regional_joint_pp_beta, protein.regional_joint_pp_se,
                joint_ld, opts.regional_prior_var_pp, opts
            );
            const auto joint_outcome_signals = fine_map_signals(
                protein.regional_joint_outcome_beta, protein.regional_joint_outcome_se,
                joint_ld, opts.regional_prior_var_outcome, opts
            );
            protein.regional_rf_signals = static_cast<int>(joint_rf_signals.size());
            compute_joint_regional_model(
                protein, joint_rf_signals, joint_protein_signals,
                joint_outcome_signals, joint_ld, opts
            );
        }
    }

    if (protein_signals.empty()) {
        protein.regional_multisignal_interpretation = "NO_PROTEIN_CREDIBLE_SET";
        return;
    }
    if (outcome_signals.empty()) {
        protein.regional_multisignal_interpretation = "NO_OUTCOME_CREDIBLE_SET";
        return;
    }

    for (int i = 0; i < static_cast<int>(protein_signals.size()); ++i) {
        for (int j = 0; j < static_cast<int>(outcome_signals.size()); ++j) {
            protein.regional_signal_pairs.push_back(coloc_signal_pair(
                protein_signals[i], outcome_signals[j], i, j, protein, ld, opts
            ));
        }
    }

    const RegionalSignalPairResult* chosen = nullptr;
    for (const auto& pair : protein.regional_signal_pairs) {
        if (pair.interpretation != "SHARED_SIGNAL_SUPPORTED") continue;
        if (chosen == nullptr || chosen->interpretation != "SHARED_SIGNAL_SUPPORTED" ||
            pair.pp_h4 > chosen->pp_h4) {
            chosen = &pair;
        }
    }
    if (chosen == nullptr) {
        for (const auto& pair : protein.regional_signal_pairs) {
            const bool is_distinct = pair.interpretation == "DISTINCT_SIGNALS_HIGH_LD" ||
                                     pair.interpretation == "DISTINCT_SIGNALS_LOW_MODERATE_LD";
            if (!is_distinct) continue;
            const bool chosen_distinct = chosen != nullptr &&
                (chosen->interpretation == "DISTINCT_SIGNALS_HIGH_LD" ||
                 chosen->interpretation == "DISTINCT_SIGNALS_LOW_MODERATE_LD");
            if (!chosen_distinct || pair.pp_h3 > chosen->pp_h3) chosen = &pair;
        }
    }
    if (chosen == nullptr) {
        for (const auto& pair : protein.regional_signal_pairs) {
            if (chosen == nullptr || pair.pp_h3 + pair.pp_h4 > chosen->pp_h3 + chosen->pp_h4) {
                chosen = &pair;
            }
        }
    }

    if (chosen != nullptr) {
        protein.regional_best_pp_shared = chosen->pp_h4;
        protein.regional_best_pp_distinct = chosen->pp_h3;
        protein.regional_best_shared_given_both = chosen->shared_given_both;
        protein.regional_best_cs_pair_r2 = chosen->max_credible_set_pair_r2;
        protein.regional_multisignal_interpretation = chosen->interpretation;
    } else {
        protein.regional_multisignal_interpretation = "SIGNAL_PAIRS_AMBIGUOUS";
    }
}

} // namespace bmediator
