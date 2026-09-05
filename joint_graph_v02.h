#ifndef BMEDIATOR_JOINT_GRAPH_V02_H
#define BMEDIATOR_JOINT_GRAPH_V02_H

#include <array>
#include <string>
#include <vector>

namespace bmediator {

struct JointGraphV02Observation {
    std::string variant;
    std::string ld_block;
    char role = 'A';
    double beta_x = 0.0;
    double se_x = 0.0;
    double beta_m = 0.0;
    double se_m = 0.0;
    double beta_y = 0.0;
    double se_y = 0.0;
    double v_x = 0.0;
    double v_m = 0.0;
    double v_y = 0.0;
    double orientation = 0.0;
    double orientation_probability = 1.0;
    double rho_xm = 0.0;
    double rho_xy = 0.0;
    double rho_my = 0.0;
};

struct JointGraphV02Options {
    double pi_xm = 0.25;
    double pi_my = 0.10;
    double pi_sparse = 0.10;
    double pi_directional = 0.10;
    double prior_sd_a = 0.70;
    double prior_sd_b = 0.70;
    double prior_sd_c = 0.175;
    double prior_sd_lambda = 0.70;
    double prior_sd_eta = 0.70;
    double q_alpha = 2.0;
    double q_beta = 2.0;
    double max_cross_block_ld = 0.05;
    double max_evidence_discrepancy = 1.0;
    double quadrature_escalation_threshold = 0.05;
    double max_quadrature_posterior_error = 0.01;
    int max_sparse_grid_level = 12;
    int min_role_blocks = 3;
    int optimizer_iterations = 1500;
    double optimizer_tolerance = 1e-6;
};

struct JointGraphV02Result {
    std::array<double, 16> state_pp{};
    std::array<double, 16> state_log_evidence{};
    std::array<double, 16> state_quadrature_difference{};
    std::array<int, 16> state_quadrature_order{};
    std::array<int, 16> state_sparse_grid_level{};
    std::array<double, 16> state_sparse_grid_cancellation{};
    std::array<double, 16> state_tensor_sparse_difference{};
    double pp_xm = 0.0;
    double pp_global_my = 0.0;
    double pp_sparse_pleio = 0.0;
    double pp_directional_pleio = 0.0;
    double pp_two_path = 0.0;
    double pp_any_pleio = 0.0;
    double log_evidence = 0.0;
    double max_ignored_ld = 0.0;
    int n_blocks = 0;
    int max_block_size = 0;
    int n_role_a = 0;
    int n_role_b = 0;
    int n_role_c = 0;
    int n_role_a_blocks = 0;
    int n_role_b_blocks = 0;
    int n_role_c_blocks = 0;
    int states_converged = 0;
    int states_regularized = 0;
    double max_adaptive_laplace_difference = 0.0;
    double max_relevant_evidence_difference = 0.0;
    double max_relevant_quadrature_difference = 0.0;
    double estimated_quadrature_posterior_error = 0.0;
    int max_quadrature_order = 3;
    int posterior_aware_refinements = 0;
    int sparse_grid_states = 0;
    int max_sparse_grid_level = 0;
    double max_sparse_grid_cancellation = 0.0;
    double max_tensor_sparse_difference = 0.0;
    JointGraphV02Options options;
};

const std::array<std::string, 16>& joint_graph_v02_state_names();

std::vector<JointGraphV02Observation> read_joint_graph_v02_tsv(
    const std::string& path);
std::vector<std::vector<double>> read_joint_graph_v02_ld(
    const std::string& path,
    const std::vector<JointGraphV02Observation>& observations);

JointGraphV02Result fit_joint_graph_v02(
    const std::vector<JointGraphV02Observation>& observations,
    const std::vector<std::vector<double>>& ld,
    const JointGraphV02Options& options = JointGraphV02Options());

double joint_graph_v02_log_likelihood(
    const std::vector<JointGraphV02Observation>& observations,
    const std::vector<std::vector<double>>& ld,
    double a, double b, double c_path, double lambda, double q, double eta,
    bool sparse_active = true,
    const JointGraphV02Options& options = JointGraphV02Options());

double joint_graph_v02_log_likelihood_integrated_q(
    const std::vector<JointGraphV02Observation>& observations,
    const std::vector<std::vector<double>>& ld,
    double a, double b, double c_path, double lambda, double eta,
    const JointGraphV02Options& options = JointGraphV02Options());

double joint_graph_v02_sparse_grid_normalized_moment(
    const std::vector<int>& powers, int level_increment);

void write_joint_graph_v02_result_tsv(const JointGraphV02Result& result,
                                      const std::string& path);

}  // namespace bmediator

#endif
