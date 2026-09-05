#include "joint_graph_v02.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace bmediator {
namespace {

using Matrix = std::vector<std::vector<double>>;

constexpr double LOG_2PI = 1.8378770664093454836;

struct StateDefinition {
    int z_xm;
    int z_my;
    int z_sparse;
    int z_directional;
};

struct PreparedBlock {
    std::vector<int> source_index;
    Matrix ld;
    Matrix kg;
    Matrix kd;
    Matrix ky;
    Matrix sampling;
    std::vector<double> observed;
    std::vector<double> directional_mean;
};

struct PreparedData {
    std::vector<PreparedBlock> blocks;
    double max_ignored_ld = 0.0;
};

struct DecodedParameters {
    double a = 0.0;
    double b = 0.0;
    double c_path = 0.0;
    double lambda = 0.0;
    double q = 0.0;
    double eta = 0.0;
};

struct StateFit {
    double log_evidence = -std::numeric_limits<double>::infinity();
    bool converged = false;
    bool regularized = false;
    double adaptive_laplace_difference = 0.0;
};

struct OptimizerResult {
    std::vector<double> point;
    double value = std::numeric_limits<double>::infinity();
    bool converged = false;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        const size_t tab = line.find('\t', start);
        fields.push_back(line.substr(start, tab - start));
        if (tab == std::string::npos) break;
        start = tab + 1;
    }
    return fields;
}

double parse_double(const std::string& value, const std::string& column) {
    size_t used = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid value in " + column + ": " + value);
    }
    if (used != value.size() || !std::isfinite(parsed)) {
        throw std::runtime_error("invalid value in " + column + ": " + value);
    }
    return parsed;
}

double log_sum_exp(const std::vector<double>& values) {
    if (values.empty()) return -std::numeric_limits<double>::infinity();
    const double maximum = *std::max_element(values.begin(), values.end());
    if (!std::isfinite(maximum)) return maximum;
    double total = 0.0;
    for (double value : values) total += std::exp(value - maximum);
    return maximum + std::log(total);
}

double log_add_exp(double first, double second) {
    const double maximum = std::max(first, second);
    if (!std::isfinite(maximum)) return maximum;
    return maximum + std::log(std::exp(first - maximum) +
                              std::exp(second - maximum));
}

bool cholesky(const Matrix& input, Matrix& lower, double& log_determinant,
              double ridge = 0.0) {
    const int n = static_cast<int>(input.size());
    lower.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double value = input[i][j];
            if (i == j) value += ridge;
            for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
            if (i == j) {
                if (!(value > 0.0) || !std::isfinite(value)) return false;
                lower[i][j] = std::sqrt(value);
            } else {
                lower[i][j] = value / lower[j][j];
            }
        }
    }
    log_determinant = 0.0;
    for (int i = 0; i < n; ++i) {
        log_determinant += 2.0 * std::log(lower[i][i]);
    }
    return true;
}

double zero_mean_mvn(const std::vector<double>& observed,
                     const std::vector<double>& mean,
                     const Matrix& covariance) {
    Matrix lower;
    double log_determinant = 0.0;
    if (!cholesky(covariance, lower, log_determinant)) {
        return -std::numeric_limits<double>::infinity();
    }
    const int n = static_cast<int>(observed.size());
    std::vector<double> solved(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double value = observed[i] - mean[i];
        for (int j = 0; j < i; ++j) value -= lower[i][j] * solved[j];
        solved[i] = value / lower[i][i];
    }
    const double quadratic = std::inner_product(
        solved.begin(), solved.end(), solved.begin(), 0.0);
    return -0.5 * (n * LOG_2PI + log_determinant + quadratic);
}

Matrix random_effect_covariance(const Matrix& ld,
                                const std::vector<double>& variance) {
    const int n = static_cast<int>(ld.size());
    Matrix result(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double value = 0.0;
            for (int k = 0; k < n; ++k) {
                value += ld[i][k] * variance[k] * ld[j][k];
            }
            result[i][j] = value;
        }
    }
    return result;
}

PreparedData prepare_data(
    const std::vector<JointGraphV02Observation>& observations,
    const Matrix& ld,
    const JointGraphV02Options& options) {
    const int n = static_cast<int>(observations.size());
    std::map<std::string, std::vector<int>> block_indices;
    for (int i = 0; i < n; ++i) {
        block_indices[observations[i].ld_block].push_back(i);
    }

    PreparedData prepared;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (observations[i].ld_block == observations[j].ld_block) continue;
            prepared.max_ignored_ld = std::max(prepared.max_ignored_ld,
                                               std::fabs(ld[i][j]));
        }
    }
    if (prepared.max_ignored_ld > options.max_cross_block_ld + 1e-12) {
        throw std::invalid_argument(
            "cross-block LD exceeds max_cross_block_ld; merge the LD blocks");
    }

    const double rho_xm = observations.front().rho_xm;
    const double rho_xy = observations.front().rho_xy;
    const double rho_my = observations.front().rho_my;
    for (const auto& item : block_indices) {
        const auto& indices = item.second;
        const int m = static_cast<int>(indices.size());
        PreparedBlock block;
        block.source_index = indices;
        block.ld.assign(m, std::vector<double>(m, 0.0));
        std::vector<double> vx(m), vm(m), vy(m);
        for (int i = 0; i < m; ++i) {
            vx[i] = observations[indices[i]].v_x;
            vm[i] = observations[indices[i]].v_m;
            vy[i] = observations[indices[i]].v_y;
            for (int j = 0; j < m; ++j) {
                block.ld[i][j] = ld[indices[i]][indices[j]];
            }
        }
        block.kg = random_effect_covariance(block.ld, vx);
        block.kd = random_effect_covariance(block.ld, vm);
        block.ky = random_effect_covariance(block.ld, vy);
        block.sampling.assign(3 * m, std::vector<double>(3 * m, 0.0));
        block.observed.assign(3 * m, 0.0);
        block.directional_mean.assign(m, 0.0);

        for (int i = 0; i < m; ++i) {
            const auto& oi = observations[indices[i]];
            block.observed[i] = oi.beta_x;
            block.observed[m + i] = oi.beta_m;
            block.observed[2 * m + i] = oi.beta_y;
            for (int k = 0; k < m; ++k) {
                block.directional_mean[i] +=
                    block.ld[i][k] * observations[indices[k]].orientation;
            }
            for (int j = 0; j < m; ++j) {
                const auto& oj = observations[indices[j]];
                const double r = block.ld[i][j];
                block.sampling[i][j] = r * oi.se_x * oj.se_x;
                block.sampling[m + i][m + j] = r * oi.se_m * oj.se_m;
                block.sampling[2 * m + i][2 * m + j] =
                    r * oi.se_y * oj.se_y;
                block.sampling[i][m + j] = rho_xm * r * oi.se_x * oj.se_m;
                block.sampling[m + j][i] = block.sampling[i][m + j];
                block.sampling[i][2 * m + j] = rho_xy * r * oi.se_x * oj.se_y;
                block.sampling[2 * m + j][i] = block.sampling[i][2 * m + j];
                block.sampling[m + i][2 * m + j] =
                    rho_my * r * oi.se_m * oj.se_y;
                block.sampling[2 * m + j][m + i] =
                    block.sampling[m + i][2 * m + j];
            }
        }
        prepared.blocks.push_back(std::move(block));
    }
    return prepared;
}

double block_component_log_likelihood(const PreparedBlock& block,
                                      const DecodedParameters& parameters,
                                      int pleiotropic) {
    const int n = static_cast<int>(block.ld.size());
    Matrix covariance = block.sampling;
    std::vector<double> mean(3 * n, 0.0);
    const double bg = parameters.c_path + parameters.a * parameters.b;
    const double bd = parameters.b + pleiotropic * parameters.lambda;
    for (int i = 0; i < n; ++i) {
        mean[2 * n + i] = parameters.eta * block.directional_mean[i];
        for (int j = 0; j < n; ++j) {
            const double kg = block.kg[i][j];
            const double kd = block.kd[i][j];
            covariance[i][j] += kg;
            covariance[i][n + j] += parameters.a * kg;
            covariance[n + j][i] = covariance[i][n + j];
            covariance[i][2 * n + j] += bg * kg;
            covariance[2 * n + j][i] = covariance[i][2 * n + j];
            covariance[n + i][n + j] += parameters.a * parameters.a * kg + kd;
            covariance[n + i][2 * n + j] +=
                parameters.a * bg * kg + bd * kd;
            covariance[2 * n + j][n + i] = covariance[n + i][2 * n + j];
            covariance[2 * n + i][2 * n + j] +=
                bg * bg * kg + bd * bd * kd + block.ky[i][j];
        }
    }
    return zero_mean_mvn(block.observed, mean, covariance);
}

double model_log_likelihood(const PreparedData& data,
                            const DecodedParameters& parameters,
                            bool sparse_active) {
    double result = 0.0;
    for (const auto& block : data.blocks) {
        const double clean = block_component_log_likelihood(block, parameters, 0);
        if (!sparse_active || parameters.lambda == 0.0) {
            result += clean;
            continue;
        }
        const double contaminated =
            block_component_log_likelihood(block, parameters, 1);
        result += log_add_exp(std::log1p(-parameters.q) + clean,
                              std::log(parameters.q) + contaminated);
    }
    return result;
}

double normal_log_density(double value, double sd) {
    return -0.5 * LOG_2PI - std::log(sd) -
           0.5 * value * value / (sd * sd);
}

DecodedParameters decode_parameters(const std::vector<double>& value,
                                    const StateDefinition& state) {
    DecodedParameters result;
    int position = 0;
    if (state.z_xm) result.a = value[position++];
    if (state.z_my) result.b = value[position++];
    result.c_path = value[position++];
    if (state.z_sparse) {
        result.lambda = value[position++];
        const double logit_q = value[position++];
        if (logit_q >= 0.0) {
            result.q = 1.0 / (1.0 + std::exp(-logit_q));
        } else {
            const double exponential = std::exp(logit_q);
            result.q = exponential / (1.0 + exponential);
        }
    }
    if (state.z_directional) result.eta = value[position++];
    return result;
}

std::vector<double> parameter_scales(const StateDefinition& state,
                                     const JointGraphV02Options& options) {
    std::vector<double> result;
    if (state.z_xm) result.push_back(options.prior_sd_a);
    if (state.z_my) result.push_back(options.prior_sd_b);
    result.push_back(options.prior_sd_c);
    if (state.z_sparse) {
        result.push_back(options.prior_sd_lambda);
        result.push_back(1.0);
    }
    if (state.z_directional) result.push_back(options.prior_sd_eta);
    return result;
}

double log_kernel(const std::vector<double>& value,
                  const StateDefinition& state,
                  const PreparedData& data,
                  const JointGraphV02Options& options) {
    for (double current : value) {
        if (!std::isfinite(current) || std::fabs(current) > 12.0) {
            return -std::numeric_limits<double>::infinity();
        }
    }
    const DecodedParameters parameters = decode_parameters(value, state);
    double result = model_log_likelihood(data, parameters, state.z_sparse != 0);
    if (!std::isfinite(result)) return result;
    if (state.z_xm) result += normal_log_density(parameters.a, options.prior_sd_a);
    if (state.z_my) result += normal_log_density(parameters.b, options.prior_sd_b);
    result += normal_log_density(parameters.c_path, options.prior_sd_c);
    if (state.z_sparse) {
        result += normal_log_density(parameters.lambda, options.prior_sd_lambda);
        const double q = parameters.q;
        result += (options.q_alpha - 1.0) * std::log(q) +
                  (options.q_beta - 1.0) * std::log1p(-q) -
                  (std::lgamma(options.q_alpha) + std::lgamma(options.q_beta) -
                   std::lgamma(options.q_alpha + options.q_beta));
        result += std::log(q) + std::log1p(-q);
    }
    if (state.z_directional) {
        result += normal_log_density(parameters.eta, options.prior_sd_eta);
    }
    return result;
}

OptimizerResult nelder_mead(const std::function<double(const std::vector<double>&)>& f,
                            const std::vector<double>& initial,
                            const std::vector<double>& scales,
                            int max_iterations, double tolerance) {
    const int dimension = static_cast<int>(initial.size());
    std::vector<std::vector<double>> simplex(dimension + 1, initial);
    for (int i = 0; i < dimension; ++i) {
        simplex[i + 1][i] += std::max(0.025, 0.20 * scales[i]);
    }
    std::vector<double> values(dimension + 1);
    for (int i = 0; i <= dimension; ++i) values[i] = f(simplex[i]);
    bool converged = false;

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::vector<int> order(dimension + 1);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int left, int right) { return values[left] < values[right]; });
        std::vector<std::vector<double>> sorted_simplex;
        std::vector<double> sorted_values;
        for (int index : order) {
            sorted_simplex.push_back(simplex[index]);
            sorted_values.push_back(values[index]);
        }
        simplex.swap(sorted_simplex);
        values.swap(sorted_values);

        double diameter = 0.0;
        for (int i = 1; i <= dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                diameter = std::max(diameter,
                                    std::fabs(simplex[i][j] - simplex[0][j]));
            }
        }
        if (diameter < tolerance &&
            std::fabs(values.back() - values.front()) < tolerance) {
            converged = true;
            break;
        }

        std::vector<double> centroid(dimension, 0.0);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) centroid[j] += simplex[i][j];
        }
        for (double& item : centroid) item /= dimension;
        auto move = [&](const std::vector<double>& from,
                        const std::vector<double>& toward, double multiplier) {
            std::vector<double> result(dimension);
            for (int j = 0; j < dimension; ++j) {
                result[j] = from[j] + multiplier * (toward[j] - from[j]);
            }
            return result;
        };

        const std::vector<double> reflected = move(centroid, simplex.back(), -1.0);
        const double reflected_value = f(reflected);
        if (reflected_value < values.front()) {
            const std::vector<double> expanded = move(centroid, reflected, 2.0);
            const double expanded_value = f(expanded);
            if (expanded_value < reflected_value) {
                simplex.back() = expanded;
                values.back() = expanded_value;
            } else {
                simplex.back() = reflected;
                values.back() = reflected_value;
            }
        } else if (reflected_value < values[dimension - 1]) {
            simplex.back() = reflected;
            values.back() = reflected_value;
        } else {
            const bool outside = reflected_value < values.back();
            const std::vector<double> contracted = outside
                ? move(centroid, reflected, 0.5)
                : move(centroid, simplex.back(), 0.5);
            const double contracted_value = f(contracted);
            if (contracted_value < (outside ? reflected_value : values.back())) {
                simplex.back() = contracted;
                values.back() = contracted_value;
            } else {
                for (int i = 1; i <= dimension; ++i) {
                    simplex[i] = move(simplex[0], simplex[i], 0.5);
                    values[i] = f(simplex[i]);
                }
            }
        }
    }

    int best = static_cast<int>(std::min_element(values.begin(), values.end()) -
                                values.begin());
    return {simplex[best], values[best], converged};
}

std::vector<double> regression_initial(
    const std::vector<JointGraphV02Observation>& observations,
    const StateDefinition& state) {
    double xx = 0.0, xm = 0.0, mm = 0.0, xy = 0.0, my = 0.0;
    for (const auto& observation : observations) {
        xx += observation.beta_x * observation.beta_x;
        xm += observation.beta_x * observation.beta_m;
        mm += observation.beta_m * observation.beta_m;
        xy += observation.beta_x * observation.beta_y;
        my += observation.beta_m * observation.beta_y;
    }
    const double a = xx > 1e-12 ? xm / xx : 0.0;
    const double determinant = xx * mm - xm * xm;
    double c_path = xx > 1e-12 ? xy / xx : 0.0;
    double b = 0.0;
    if (determinant > 1e-12) {
        c_path = (xy * mm - my * xm) / determinant;
        b = (my * xx - xy * xm) / determinant;
    }
    auto limit = [](double value) { return std::max(-2.0, std::min(2.0, value)); };
    std::vector<double> result;
    if (state.z_xm) result.push_back(limit(a));
    if (state.z_my) result.push_back(limit(b));
    result.push_back(limit(c_path));
    if (state.z_sparse) {
        result.push_back(0.25);
        result.push_back(std::log(0.35 / 0.65));
    }
    if (state.z_directional) result.push_back(0.0);
    return result;
}

OptimizerResult optimize_state(
    const std::vector<JointGraphV02Observation>& observations,
    const StateDefinition& state, const PreparedData& data,
    const JointGraphV02Options& options) {
    const std::vector<double> scales = parameter_scales(state, options);
    std::vector<std::vector<double>> starts;
    starts.push_back(regression_initial(observations, state));
    starts.push_back(std::vector<double>(scales.size(), 0.0));
    if (state.z_sparse) {
        int q_index = static_cast<int>(scales.size()) - 1 - state.z_directional;
        starts.back()[q_index] = std::log(0.35 / 0.65);
    }
    for (int i = 0; i < static_cast<int>(scales.size()); ++i) {
        std::vector<double> start = starts.front();
        start[i] += 0.5 * scales[i];
        starts.push_back(start);
        start = starts.front();
        start[i] -= 0.5 * scales[i];
        starts.push_back(start);
    }

    auto objective = [&](const std::vector<double>& value) {
        const double kernel = log_kernel(value, state, data, options);
        return std::isfinite(kernel) ? -kernel : 1e100;
    };
    OptimizerResult best;
    for (const auto& start : starts) {
        OptimizerResult current = nelder_mead(
            objective, start, scales, options.optimizer_iterations,
            options.optimizer_tolerance);
        if (current.value < best.value) best = current;
    }
    return best;
}

Matrix numerical_hessian(const std::function<double(const std::vector<double>&)>& f,
                         const std::vector<double>& point) {
    const int n = static_cast<int>(point.size());
    Matrix result(n, std::vector<double>(n, 0.0));
    std::vector<double> step(n);
    for (int i = 0; i < n; ++i) step[i] = 2e-4 * (1.0 + std::fabs(point[i]));
    const double center = f(point);
    for (int i = 0; i < n; ++i) {
        std::vector<double> plus = point, minus = point;
        plus[i] += step[i];
        minus[i] -= step[i];
        result[i][i] = (f(plus) - 2.0 * center + f(minus)) /
                       (step[i] * step[i]);
        for (int j = 0; j < i; ++j) {
            std::vector<double> pp = point, pm = point, mp = point, mm = point;
            pp[i] += step[i]; pp[j] += step[j];
            pm[i] += step[i]; pm[j] -= step[j];
            mp[i] -= step[i]; mp[j] += step[j];
            mm[i] -= step[i]; mm[j] -= step[j];
            result[i][j] = (f(pp) - f(pm) - f(mp) + f(mm)) /
                           (4.0 * step[i] * step[j]);
            result[j][i] = result[i][j];
        }
    }
    return result;
}

StateFit integrate_state(
    const std::vector<JointGraphV02Observation>& observations,
    const StateDefinition& state, const PreparedData& data,
    const JointGraphV02Options& options) {
    const OptimizerResult mode = optimize_state(observations, state, data, options);
    auto negative_kernel = [&](const std::vector<double>& value) {
        const double kernel = log_kernel(value, state, data, options);
        return std::isfinite(kernel) ? -kernel : 1e100;
    };
    Matrix hessian = numerical_hessian(negative_kernel, mode.point);
    Matrix precision_lower;
    double log_precision_determinant = 0.0;
    bool regularized = false;
    bool factorized = false;
    for (int attempt = 0; attempt < 10; ++attempt) {
        const double ridge = attempt == 0 ? 0.0 : std::pow(10.0, attempt - 7);
        if (cholesky(hessian, precision_lower, log_precision_determinant, ridge)) {
            factorized = true;
            regularized = attempt > 0;
            break;
        }
    }
    if (!factorized) {
        throw std::runtime_error("adaptive integration Hessian is not positive definite");
    }

    constexpr std::array<double, 3> nodes{{
        -1.224744871391589, 0.0, 1.224744871391589
    }};
    constexpr std::array<double, 3> weights{{
        0.2954089751509193, 1.1816359006036772, 0.2954089751509193
    }};
    const int dimension = static_cast<int>(mode.point.size());
    std::vector<double> log_terms;
    int count = 1;
    for (int i = 0; i < dimension; ++i) count *= 3;
    log_terms.reserve(count);
    std::vector<double> z(dimension, 0.0);

    std::function<void(int, double, double)> visit = [&](int position,
                                                         double log_weight,
                                                         double z_squared) {
        if (position < dimension) {
            for (int i = 0; i < 3; ++i) {
                z[position] = nodes[i];
                visit(position + 1, log_weight + std::log(weights[i]),
                      z_squared + nodes[i] * nodes[i]);
            }
            return;
        }
        std::vector<double> delta(dimension, 0.0);
        for (int i = dimension - 1; i >= 0; --i) {
            double value = std::sqrt(2.0) * z[i];
            for (int j = i + 1; j < dimension; ++j) {
                value -= precision_lower[j][i] * delta[j];
            }
            delta[i] = value / precision_lower[i][i];
        }
        std::vector<double> point = mode.point;
        for (int i = 0; i < dimension; ++i) point[i] += delta[i];
        log_terms.push_back(log_kernel(point, state, data, options) +
                            log_weight + z_squared);
    };
    visit(0, 0.0, 0.0);

    StateFit result;
    result.log_evidence = log_sum_exp(log_terms) +
        0.5 * dimension * std::log(2.0) - 0.5 * log_precision_determinant;
    const double laplace_evidence = -mode.value +
        0.5 * dimension * LOG_2PI - 0.5 * log_precision_determinant;
    result.adaptive_laplace_difference =
        std::fabs(result.log_evidence - laplace_evidence);
    result.converged = mode.converged;
    result.regularized = regularized;
    return result;
}

double bernoulli_log_probability(int indicator, double probability) {
    return indicator ? std::log(probability) : std::log1p(-probability);
}

double state_log_prior(const StateDefinition& state,
                       const JointGraphV02Options& options) {
    return bernoulli_log_probability(state.z_xm, options.pi_xm) +
           bernoulli_log_probability(state.z_my, options.pi_my) +
           bernoulli_log_probability(state.z_sparse, options.pi_sparse) +
           bernoulli_log_probability(state.z_directional,
                                     options.pi_directional);
}

std::array<StateDefinition, 16> state_definitions() {
    std::array<StateDefinition, 16> result{};
    for (int index = 0; index < 16; ++index) {
        result[index] = {index & 1, (index >> 1) & 1,
                         (index >> 2) & 1, (index >> 3) & 1};
    }
    return result;
}

void validate_inputs(const std::vector<JointGraphV02Observation>& observations,
                     const Matrix& ld,
                     const JointGraphV02Options& options) {
    if (observations.empty()) throw std::invalid_argument("input contains no variants");
    const int n = static_cast<int>(observations.size());
    if (static_cast<int>(ld.size()) != n) {
        throw std::invalid_argument("LD matrix dimension does not match input");
    }
    std::set<std::string> variants;
    for (int i = 0; i < n; ++i) {
        if (static_cast<int>(ld[i].size()) != n) {
            throw std::invalid_argument("LD matrix must be square");
        }
        if (!variants.insert(observations[i].variant).second) {
            throw std::invalid_argument("duplicate variant: " + observations[i].variant);
        }
        const auto& observation = observations[i];
        if (observation.role != 'A' && observation.role != 'B' &&
            observation.role != 'C') {
            throw std::invalid_argument("role must be A, B, or C");
        }
        if (observation.ld_block.empty()) {
            throw std::invalid_argument("ld_block must not be empty");
        }
        if (!(observation.se_x > 0.0 && observation.se_m > 0.0 &&
              observation.se_y > 0.0 && observation.v_x >= 0.0 &&
              observation.v_m >= 0.0 && observation.v_y > 0.0)) {
            throw std::invalid_argument("SEs and external variances are invalid");
        }
        if (std::fabs(std::fabs(observation.orientation) - 1.0) > 1e-12) {
            throw std::invalid_argument(
                "orientation must be -1 or 1 and derived independently");
        }
        if (std::fabs(ld[i][i] - 1.0) > 1e-6) {
            throw std::invalid_argument("LD diagonal must equal one");
        }
        for (int j = 0; j < n; ++j) {
            if (!std::isfinite(ld[i][j]) || std::fabs(ld[i][j]) > 1.0 + 1e-8 ||
                std::fabs(ld[i][j] - ld[j][i]) > 1e-6) {
                throw std::invalid_argument("LD matrix is invalid or asymmetric");
            }
        }
        if (i > 0 &&
            (std::fabs(observation.rho_xm - observations[0].rho_xm) > 1e-12 ||
             std::fabs(observation.rho_xy - observations[0].rho_xy) > 1e-12 ||
             std::fabs(observation.rho_my - observations[0].rho_my) > 1e-12)) {
            throw std::invalid_argument(
                "sampling correlations must be constant within an analysis");
        }
    }
    const double r12 = observations[0].rho_xm;
    const double r13 = observations[0].rho_xy;
    const double r23 = observations[0].rho_my;
    const double correlation_determinant =
        1.0 + 2.0 * r12 * r13 * r23 - r12 * r12 - r13 * r13 - r23 * r23;
    if (!(correlation_determinant > 1e-8) || std::fabs(r12) >= 1.0 ||
        std::fabs(r13) >= 1.0 || std::fabs(r23) >= 1.0) {
        throw std::invalid_argument("sampling correlation matrix is not positive definite");
    }
    const std::array<double, 4> probabilities{{
        options.pi_xm, options.pi_my, options.pi_sparse, options.pi_directional
    }};
    for (double probability : probabilities) {
        if (!(probability > 0.0 && probability < 1.0)) {
            throw std::invalid_argument("state priors must be in (0,1)");
        }
    }
    if (!(options.q_alpha > 1.0 && options.q_beta > 1.0)) {
        throw std::invalid_argument("q beta-prior shapes must exceed one");
    }
    if (!(options.prior_sd_a > 0.0 && options.prior_sd_b > 0.0 &&
          options.prior_sd_c > 0.0 && options.prior_sd_lambda > 0.0 &&
          options.prior_sd_eta > 0.0 && options.optimizer_iterations >= 100 &&
          options.optimizer_tolerance > 0.0)) {
        throw std::invalid_argument("continuous priors or optimizer options are invalid");
    }
}

}  // namespace

const std::array<std::string, 16>& joint_graph_v02_state_names() {
    static const std::array<std::string, 16> names{{
        "S0000", "S1000", "S0100", "S1100",
        "S0010", "S1010", "S0110", "S1110",
        "S0001", "S1001", "S0101", "S1101",
        "S0011", "S1011", "S0111", "S1111"
    }};
    return names;
}

std::vector<JointGraphV02Observation> read_joint_graph_v02_tsv(
    const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open input: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("input is empty");
    const auto header = split_tab(line);
    std::unordered_map<std::string, int> column;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) column[header[i]] = i;
    const std::array<std::string, 16> required{{
        "variant", "ld_block", "role", "beta_x", "se_x", "beta_m", "se_m",
        "beta_y", "se_y", "v_x", "v_m", "v_y", "orientation", "rho_xm",
        "rho_xy", "rho_my"
    }};
    for (const auto& name : required) {
        if (!column.count(name)) throw std::runtime_error("missing column: " + name);
    }

    std::vector<JointGraphV02Observation> result;
    int line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() != header.size()) {
            throw std::runtime_error("wrong field count at line " +
                                     std::to_string(line_number));
        }
        JointGraphV02Observation observation;
        observation.variant = fields[column["variant"]];
        observation.ld_block = fields[column["ld_block"]];
        const std::string role = fields[column["role"]];
        if (role.size() != 1) throw std::runtime_error("invalid role");
        observation.role = role[0];
        observation.beta_x = parse_double(fields[column["beta_x"]], "beta_x");
        observation.se_x = parse_double(fields[column["se_x"]], "se_x");
        observation.beta_m = parse_double(fields[column["beta_m"]], "beta_m");
        observation.se_m = parse_double(fields[column["se_m"]], "se_m");
        observation.beta_y = parse_double(fields[column["beta_y"]], "beta_y");
        observation.se_y = parse_double(fields[column["se_y"]], "se_y");
        observation.v_x = parse_double(fields[column["v_x"]], "v_x");
        observation.v_m = parse_double(fields[column["v_m"]], "v_m");
        observation.v_y = parse_double(fields[column["v_y"]], "v_y");
        observation.orientation =
            parse_double(fields[column["orientation"]], "orientation");
        observation.rho_xm = parse_double(fields[column["rho_xm"]], "rho_xm");
        observation.rho_xy = parse_double(fields[column["rho_xy"]], "rho_xy");
        observation.rho_my = parse_double(fields[column["rho_my"]], "rho_my");
        result.push_back(std::move(observation));
    }
    return result;
}

Matrix read_joint_graph_v02_ld(
    const std::string& path,
    const std::vector<JointGraphV02Observation>& observations) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open LD matrix: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("LD matrix is empty");
    const auto header = split_tab(line);
    if (header.size() != observations.size() + 1) {
        throw std::runtime_error("LD header dimension does not match input");
    }
    std::unordered_map<std::string, int> source_column;
    for (int i = 1; i < static_cast<int>(header.size()); ++i) {
        source_column[header[i]] = i - 1;
    }
    const int n = static_cast<int>(observations.size());
    Matrix source(n, std::vector<double>(n, 0.0));
    std::unordered_map<std::string, int> source_row;
    int row = 0;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (static_cast<int>(fields.size()) != n + 1 || row >= n) {
            throw std::runtime_error("LD row dimension does not match input");
        }
        source_row[fields[0]] = row;
        for (int j = 0; j < n; ++j) {
            source[row][j] = parse_double(fields[j + 1], "LD");
        }
        ++row;
    }
    if (row != n) throw std::runtime_error("LD matrix has wrong number of rows");
    Matrix result(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        if (!source_row.count(observations[i].variant) ||
            !source_column.count(observations[i].variant)) {
            throw std::runtime_error("variant is missing from LD matrix: " +
                                     observations[i].variant);
        }
        for (int j = 0; j < n; ++j) {
            result[i][j] = source[source_row[observations[i].variant]]
                                 [source_column[observations[j].variant]];
        }
    }
    return result;
}

JointGraphV02Result fit_joint_graph_v02(
    const std::vector<JointGraphV02Observation>& observations,
    const Matrix& ld,
    const JointGraphV02Options& options) {
    validate_inputs(observations, ld, options);
    const PreparedData data = prepare_data(observations, ld, options);
    const auto states = state_definitions();
    std::array<StateFit, 16> fits;
    std::vector<double> log_joint(16);
    for (int i = 0; i < 16; ++i) {
        fits[i] = integrate_state(observations, states[i], data, options);
        log_joint[i] = fits[i].log_evidence + state_log_prior(states[i], options);
    }

    JointGraphV02Result result;
    result.log_evidence = log_sum_exp(log_joint);
    result.options = options;
    result.max_ignored_ld = data.max_ignored_ld;
    result.n_blocks = static_cast<int>(data.blocks.size());
    for (const auto& block : data.blocks) {
        result.max_block_size = std::max(
            result.max_block_size, static_cast<int>(block.ld.size()));
    }
    for (int i = 0; i < 16; ++i) {
        result.state_log_evidence[i] = fits[i].log_evidence;
        result.state_pp[i] = std::exp(log_joint[i] - result.log_evidence);
        result.states_converged += fits[i].converged;
        result.states_regularized += fits[i].regularized;
        result.max_adaptive_laplace_difference = std::max(
            result.max_adaptive_laplace_difference,
            fits[i].adaptive_laplace_difference);
        if (states[i].z_xm) result.pp_xm += result.state_pp[i];
        if (states[i].z_my) result.pp_global_my += result.state_pp[i];
        if (states[i].z_sparse) result.pp_sparse_pleio += result.state_pp[i];
        if (states[i].z_directional) {
            result.pp_directional_pleio += result.state_pp[i];
        }
        if (states[i].z_xm && states[i].z_my) {
            result.pp_two_path += result.state_pp[i];
        }
        if (states[i].z_sparse || states[i].z_directional) {
            result.pp_any_pleio += result.state_pp[i];
        }
    }
    if (result.states_converged != 16) {
        throw std::runtime_error("one or more graph-state optimizations did not converge");
    }
    if (result.states_regularized != 0 ||
        result.max_adaptive_laplace_difference > 1.0) {
        throw std::runtime_error(
            "adaptive evidence diagnostic failed; no posterior is reportable");
    }
    return result;
}

double joint_graph_v02_log_likelihood(
    const std::vector<JointGraphV02Observation>& observations,
    const Matrix& ld,
    double a, double b, double c_path, double lambda, double q, double eta,
    bool sparse_active,
    const JointGraphV02Options& options) {
    validate_inputs(observations, ld, options);
    if (sparse_active && !(q > 0.0 && q < 1.0)) {
        throw std::invalid_argument("q must be in (0,1) for sparse likelihood");
    }
    const PreparedData data = prepare_data(observations, ld, options);
    DecodedParameters parameters;
    parameters.a = a;
    parameters.b = b;
    parameters.c_path = c_path;
    parameters.lambda = lambda;
    parameters.q = sparse_active ? q : 0.0;
    parameters.eta = eta;
    return model_log_likelihood(data, parameters, sparse_active);
}

void write_joint_graph_v02_result_tsv(const JointGraphV02Result& result,
                                      const std::string& path) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open output: " + path);
    output << "model_version\tidentification_scope";
    for (const auto& state : joint_graph_v02_state_names()) {
        output << "\tPP_" << state;
    }
    output << "\tPP_XM\tPP_global_MY\tPP_sparse_P\tPP_directional_P"
           << "\tPP_any_P\tPP_two_path\tlog_evidence\tn_blocks"
           << "\tmax_block_size\tmax_ignored_ld\tstates_converged"
           << "\tstates_regularized\tmax_adaptive_laplace_difference"
           << "\tpi_xm\tpi_my\tpi_sparse\tpi_directional\n";
    output << "JG-0.2\tCONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY"
           << std::setprecision(17);
    for (double value : result.state_pp) output << '\t' << value;
    output << '\t' << result.pp_xm
           << '\t' << result.pp_global_my
           << '\t' << result.pp_sparse_pleio
           << '\t' << result.pp_directional_pleio
           << '\t' << result.pp_any_pleio
           << '\t' << result.pp_two_path
           << '\t' << result.log_evidence
           << '\t' << result.n_blocks
           << '\t' << result.max_block_size
           << '\t' << result.max_ignored_ld
           << '\t' << result.states_converged
           << '\t' << result.states_regularized
           << '\t' << result.max_adaptive_laplace_difference
           << '\t' << result.options.pi_xm
           << '\t' << result.options.pi_my
           << '\t' << result.options.pi_sparse
           << '\t' << result.options.pi_directional << '\n';
}

}  // namespace bmediator
