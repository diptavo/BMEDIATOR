#include "joint_graph_model.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace bmediator {
namespace {

constexpr double LOG_2PI = 1.8378770664093454836;

struct StateDefinition {
    int z_xm;
    int z_my;
    int z_pleio;
};

constexpr std::array<StateDefinition, 8> STATES{{
    {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}
}};

struct Quadrature {
    std::array<double, 7> value{};
    std::array<double, 7> log_weight{};
};

struct GridView {
    const double* value;
    const double* log_weight;
    int size;
};

struct IntegrationTerm {
    double log_value;
    double a;
    double b;
    double c_path;
    double lambda;
    double q;
};

struct StateFit {
    double log_evidence = -std::numeric_limits<double>::infinity();
    double mean_a = 0.0;
    double mean_b = 0.0;
    double mean_c = 0.0;
    double mean_lambda = 0.0;
    double mean_q = 0.0;
    double mean_indirect = 0.0;
};

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

Quadrature gauss_hermite_7(double sd) {
    constexpr std::array<double, 7> nodes{{
        -2.6519613568352335, -1.6735516287674714, -0.8162878828589647,
         0.0,
         0.8162878828589647,  1.6735516287674714,  2.6519613568352335
    }};
    constexpr std::array<double, 7> weights{{
        0.0009717812450995192, 0.05451558281912703,
        0.4256072526101278, 0.8102646175568073,
        0.4256072526101278, 0.05451558281912703,
        0.0009717812450995192
    }};
    Quadrature result;
    const double scale = std::sqrt(2.0) * sd;
    const double log_sqrt_pi = 0.5 * std::log(std::acos(-1.0));
    for (int index = 0; index < 7; ++index) {
        result.value[index] = scale * nodes[index];
        result.log_weight[index] = std::log(weights[index]) - log_sqrt_pi;
    }
    return result;
}

int role_index(char role) {
    if (role == 'A') return 0;
    if (role == 'B') return 1;
    if (role == 'C') return 2;
    throw std::invalid_argument("joint-graph role must be A, B, or C");
}

double log_mvn3_zero(const JointGraphObservation& observation,
                     double s11, double s12, double s13,
                     double s22, double s23, double s33) {
    const double determinant =
        s11 * (s22 * s33 - s23 * s23) -
        s12 * (s12 * s33 - s13 * s23) +
        s13 * (s12 * s23 - s13 * s22);
    if (!(determinant > 0.0) || !std::isfinite(determinant)) {
        return -std::numeric_limits<double>::infinity();
    }

    const double i11 = (s22 * s33 - s23 * s23) / determinant;
    const double i12 = (s13 * s23 - s12 * s33) / determinant;
    const double i13 = (s12 * s23 - s13 * s22) / determinant;
    const double i22 = (s11 * s33 - s13 * s13) / determinant;
    const double i23 = (s12 * s13 - s11 * s23) / determinant;
    const double i33 = (s11 * s22 - s12 * s12) / determinant;

    const double x = observation.beta_x;
    const double m = observation.beta_m;
    const double y = observation.beta_y;
    const double quadratic =
        i11 * x * x + i22 * m * m + i33 * y * y +
        2.0 * i12 * x * m + 2.0 * i13 * x * y + 2.0 * i23 * m * y;
    return -0.5 * (3.0 * LOG_2PI + std::log(determinant) + quadratic);
}

double component_log_likelihood(const JointGraphObservation& observation,
                                double a, double b, double c_path,
                                double lambda, int h,
                                const JointGraphOptions& options) {
    const int role = role_index(observation.role);
    const double vx = options.vx[role];
    const double vm = options.vm[role];
    const double bg = c_path + a * b;
    const double bd = b + h * lambda;
    const double s11 = vx + observation.se_x * observation.se_x;
    const double s12 = a * vx;
    const double s13 = bg * vx;
    const double s22 = a * a * vx + vm +
                       observation.se_m * observation.se_m;
    const double s23 = a * bg * vx + bd * vm;
    const double s33 = bg * bg * vx + bd * bd * vm + options.vy +
                       observation.se_y * observation.se_y;
    return log_mvn3_zero(observation, s11, s12, s13, s22, s23, s33);
}

double log_likelihood(const std::vector<JointGraphObservation>& observations,
                      double a, double b, double c_path,
                      double lambda, double q,
                      const JointGraphOptions& options) {
    double total = 0.0;
    for (const auto& observation : observations) {
        const double clean = component_log_likelihood(
            observation, a, b, c_path, lambda, 0, options);
        if (q <= 0.0 || lambda == 0.0) {
            total += clean;
            continue;
        }
        const double shared = component_log_likelihood(
            observation, a, b, c_path, lambda, 1, options);
        total += log_add_exp(std::log1p(-q) + clean, std::log(q) + shared);
    }
    return total;
}

double bernoulli_log_probability(int indicator, double probability) {
    return indicator ? std::log(probability) : std::log1p(-probability);
}

double state_log_prior(const StateDefinition& state,
                       const JointGraphOptions& options) {
    return bernoulli_log_probability(state.z_xm, options.pi_xm) +
           bernoulli_log_probability(state.z_my, options.pi_my) +
           bernoulli_log_probability(state.z_pleio, options.pi_pleio);
}

StateFit integrate_state(const std::vector<JointGraphObservation>& observations,
                         const StateDefinition& state,
                         const JointGraphOptions& options,
                         const Quadrature& qa,
                         const Quadrature& qb,
                         const Quadrature& qc,
                         const Quadrature& qlambda) {
    const double zero = 0.0;
    const double zero_log_weight = 0.0;
    const GridView a_grid = state.z_xm
        ? GridView{qa.value.data(), qa.log_weight.data(), 7}
        : GridView{&zero, &zero_log_weight, 1};
    const GridView b_grid = state.z_my
        ? GridView{qb.value.data(), qb.log_weight.data(), 7}
        : GridView{&zero, &zero_log_weight, 1};
    const GridView c_grid{qc.value.data(), qc.log_weight.data(), 7};
    const GridView lambda_grid = state.z_pleio
        ? GridView{qlambda.value.data(), qlambda.log_weight.data(), 7}
        : GridView{&zero, &zero_log_weight, 1};
    const int q_size = state.z_pleio ? 3 : 1;
    const int count = a_grid.size * b_grid.size * c_grid.size *
                      lambda_grid.size * q_size;
    std::vector<IntegrationTerm> terms;
    terms.reserve(count);

    for (int ia = 0; ia < a_grid.size; ++ia) {
        for (int ib = 0; ib < b_grid.size; ++ib) {
            for (int ic = 0; ic < c_grid.size; ++ic) {
                for (int il = 0; il < lambda_grid.size; ++il) {
                    for (int iq = 0; iq < q_size; ++iq) {
                        const double q = state.z_pleio ? options.q_values[iq] : 0.0;
                        const double q_log_weight = state.z_pleio
                            ? -std::log(3.0) : 0.0;
                        const double a = a_grid.value[ia];
                        const double b = b_grid.value[ib];
                        const double c_path = c_grid.value[ic];
                        const double lambda = lambda_grid.value[il];
                        const double log_value =
                            a_grid.log_weight[ia] + b_grid.log_weight[ib] +
                            c_grid.log_weight[ic] + lambda_grid.log_weight[il] +
                            q_log_weight + log_likelihood(
                                observations, a, b, c_path, lambda, q, options);
                        terms.push_back({log_value, a, b, c_path, lambda, q});
                    }
                }
            }
        }
    }

    std::vector<double> log_values;
    log_values.reserve(terms.size());
    for (const auto& term : terms) log_values.push_back(term.log_value);
    StateFit fit;
    fit.log_evidence = log_sum_exp(log_values);
    for (const auto& term : terms) {
        const double weight = std::exp(term.log_value - fit.log_evidence);
        fit.mean_a += weight * term.a;
        fit.mean_b += weight * term.b;
        fit.mean_c += weight * term.c_path;
        fit.mean_lambda += weight * term.lambda;
        fit.mean_q += weight * term.q;
        fit.mean_indirect += weight * term.a * term.b;
    }
    return fit;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

double parse_double(const std::string& value, const std::string& column) {
    size_t used = 0;
    const double parsed = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(parsed)) {
        throw std::runtime_error("invalid value in " + column + ": " + value);
    }
    return parsed;
}

}  // namespace

const std::array<std::string, 8>& joint_graph_state_names() {
    static const std::array<std::string, 8> names{{
        "S000", "S100", "S010", "S110",
        "S001", "S101", "S011", "S111"
    }};
    return names;
}

JointGraphResult fit_joint_graph(
    const std::vector<JointGraphObservation>& observations,
    const JointGraphOptions& options) {
    if (observations.empty()) {
        throw std::invalid_argument("joint-graph input contains no variants");
    }
    if (!(options.pi_xm > 0.0 && options.pi_xm < 1.0) ||
        !(options.pi_my > 0.0 && options.pi_my < 1.0) ||
        !(options.pi_pleio > 0.0 && options.pi_pleio < 1.0)) {
        throw std::invalid_argument("joint-graph inclusion priors must be in (0,1)");
    }
    if (!(options.prior_sd_a > 0.0) || !(options.prior_sd_b > 0.0) ||
        !(options.prior_sd_c > 0.0) || !(options.prior_sd_lambda > 0.0) ||
        !(options.vy > 0.0)) {
        throw std::invalid_argument("joint-graph prior scales must be positive");
    }
    for (int role = 0; role < 3; ++role) {
        if (!(options.vx[role] > 0.0) || !(options.vm[role] > 0.0)) {
            throw std::invalid_argument("joint-graph role variances must be positive");
        }
        if (!(options.q_values[role] > 0.0 && options.q_values[role] < 1.0)) {
            throw std::invalid_argument("joint-graph q values must be in (0,1)");
        }
    }
    for (const auto& observation : observations) {
        role_index(observation.role);
        if (!(observation.se_x > 0.0) || !(observation.se_m > 0.0) ||
            !(observation.se_y > 0.0)) {
            throw std::invalid_argument("joint-graph standard errors must be positive");
        }
        if (!std::isfinite(observation.beta_x) ||
            !std::isfinite(observation.beta_m) ||
            !std::isfinite(observation.beta_y) ||
            !std::isfinite(observation.se_x) ||
            !std::isfinite(observation.se_m) ||
            !std::isfinite(observation.se_y)) {
            throw std::invalid_argument("joint-graph observations must be finite");
        }
    }

    const Quadrature qa = gauss_hermite_7(options.prior_sd_a);
    const Quadrature qb = gauss_hermite_7(options.prior_sd_b);
    const Quadrature qc = gauss_hermite_7(options.prior_sd_c);
    const Quadrature qlambda = gauss_hermite_7(options.prior_sd_lambda);
    std::array<StateFit, 8> fits;
    std::vector<double> log_joint(8);
    for (int index = 0; index < 8; ++index) {
        fits[index] = integrate_state(
            observations, STATES[index], options, qa, qb, qc, qlambda);
        log_joint[index] = fits[index].log_evidence +
                           state_log_prior(STATES[index], options);
    }

    JointGraphResult result;
    result.log_evidence = log_sum_exp(log_joint);
    for (int index = 0; index < 8; ++index) {
        result.state_log_evidence[index] = fits[index].log_evidence;
        result.state_pp[index] = std::exp(log_joint[index] - result.log_evidence);
        if (STATES[index].z_xm) result.pp_xm += result.state_pp[index];
        if (STATES[index].z_my) result.pp_global_my += result.state_pp[index];
        if (STATES[index].z_pleio) {
            result.pp_nonaligned_pleio += result.state_pp[index];
        }
        if (STATES[index].z_xm && STATES[index].z_my) {
            result.pp_two_path += result.state_pp[index];
        }
        result.mean_c += result.state_pp[index] * fits[index].mean_c;
    }
    result.pp_two_path_plus_pleio = result.state_pp[7];

    const double nan = std::numeric_limits<double>::quiet_NaN();
    result.mean_a_given_xm = result.pp_xm > 0.0 ? 0.0 : nan;
    result.mean_b_given_my = result.pp_global_my > 0.0 ? 0.0 : nan;
    result.mean_lambda_given_pleio = result.pp_nonaligned_pleio > 0.0 ? 0.0 : nan;
    result.mean_q_given_pleio = result.pp_nonaligned_pleio > 0.0 ? 0.0 : nan;
    result.mean_indirect_given_two_path = result.pp_two_path > 0.0 ? 0.0 : nan;
    for (int index = 0; index < 8; ++index) {
        if (STATES[index].z_xm && result.pp_xm > 0.0) {
            result.mean_a_given_xm += result.state_pp[index] * fits[index].mean_a /
                                      result.pp_xm;
        }
        if (STATES[index].z_my && result.pp_global_my > 0.0) {
            result.mean_b_given_my += result.state_pp[index] * fits[index].mean_b /
                                      result.pp_global_my;
        }
        if (STATES[index].z_pleio && result.pp_nonaligned_pleio > 0.0) {
            result.mean_lambda_given_pleio +=
                result.state_pp[index] * fits[index].mean_lambda /
                result.pp_nonaligned_pleio;
            result.mean_q_given_pleio += result.state_pp[index] * fits[index].mean_q /
                                         result.pp_nonaligned_pleio;
        }
        if (STATES[index].z_xm && STATES[index].z_my && result.pp_two_path > 0.0) {
            result.mean_indirect_given_two_path +=
                result.state_pp[index] * fits[index].mean_indirect /
                result.pp_two_path;
        }
    }
    return result;
}

std::vector<JointGraphObservation> read_joint_graph_tsv(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open joint-graph input: " + path);
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("joint-graph input is empty: " + path);
    }
    const std::vector<std::string> header = split_tab(line);
    std::unordered_map<std::string, int> column;
    for (int index = 0; index < static_cast<int>(header.size()); ++index) {
        column[header[index]] = index;
    }
    const std::array<std::string, 8> required{{
        "variant", "role", "beta_x", "se_x", "beta_m", "se_m", "beta_y", "se_y"
    }};
    for (const auto& name : required) {
        if (!column.count(name)) {
            throw std::runtime_error("missing joint-graph column: " + name);
        }
    }

    std::vector<JointGraphObservation> observations;
    int line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_tab(line);
        if (fields.size() != header.size()) {
            throw std::runtime_error("wrong field count at line " +
                                     std::to_string(line_number));
        }
        JointGraphObservation observation;
        observation.variant = fields[column["variant"]];
        const std::string& role = fields[column["role"]];
        if (role.size() != 1) {
            throw std::runtime_error("invalid role at line " +
                                     std::to_string(line_number));
        }
        observation.role = role[0];
        observation.beta_x = parse_double(fields[column["beta_x"]], "beta_x");
        observation.se_x = parse_double(fields[column["se_x"]], "se_x");
        observation.beta_m = parse_double(fields[column["beta_m"]], "beta_m");
        observation.se_m = parse_double(fields[column["se_m"]], "se_m");
        observation.beta_y = parse_double(fields[column["beta_y"]], "beta_y");
        observation.se_y = parse_double(fields[column["se_y"]], "se_y");
        observations.push_back(observation);
    }
    return observations;
}

void write_joint_graph_result_tsv(const JointGraphResult& result,
                                  const std::string& path) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open joint-graph output: " + path);
    output << "model_version\tidentification_scope";
    for (const auto& name : joint_graph_state_names()) output << "\tPP_" << name;
    output << "\tPP_XM\tPP_global_MY\tPP_nonaligned_P\tPP_two_path"
           << "\tPP_two_path_plus_P\tmean_a_given_XM\tmean_b_given_MY"
           << "\tmean_c\tmean_lambda_given_P\tmean_q_given_P"
           << "\tmean_indirect_given_two_path\tlog_evidence\n";
    output << "JG-0.1\tCONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY"
           << std::setprecision(17);
    for (double value : result.state_pp) output << '\t' << value;
    output << '\t' << result.pp_xm
           << '\t' << result.pp_global_my
           << '\t' << result.pp_nonaligned_pleio
           << '\t' << result.pp_two_path
           << '\t' << result.pp_two_path_plus_pleio
           << '\t' << result.mean_a_given_xm
           << '\t' << result.mean_b_given_my
           << '\t' << result.mean_c
           << '\t' << result.mean_lambda_given_pleio
           << '\t' << result.mean_q_given_pleio
           << '\t' << result.mean_indirect_given_two_path
           << '\t' << result.log_evidence << '\n';
}

}  // namespace bmediator
