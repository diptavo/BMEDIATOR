#ifndef BMEDIATOR_JOINT_GRAPH_MODEL_H
#define BMEDIATOR_JOINT_GRAPH_MODEL_H

#include <array>
#include <string>
#include <vector>

namespace bmediator {

struct JointGraphObservation {
    std::string variant;
    char role = 'A';
    double beta_x = 0.0;
    double se_x = 0.0;
    double beta_m = 0.0;
    double se_m = 0.0;
    double beta_y = 0.0;
    double se_y = 0.0;
};

struct JointGraphOptions {
    double pi_xm = 0.25;
    double pi_my = 0.10;
    double pi_pleio = 0.10;
    double prior_sd_a = 0.70;
    double prior_sd_b = 0.70;
    double prior_sd_c = 0.175;
    double prior_sd_lambda = 0.70;
    std::array<double, 3> q_values{{0.15, 0.35, 0.60}};
    std::array<double, 3> vx{{0.0400, 0.0004, 0.0400}};
    std::array<double, 3> vm{{0.0025, 0.0400, 0.0400}};
    double vy = 0.0025;
};

struct JointGraphResult {
    std::array<double, 8> state_pp{};
    std::array<double, 8> state_log_evidence{};
    double pp_xm = 0.0;
    double pp_global_my = 0.0;
    double pp_nonaligned_pleio = 0.0;
    double pp_two_path = 0.0;
    double pp_two_path_plus_pleio = 0.0;
    double mean_a_given_xm = 0.0;
    double mean_b_given_my = 0.0;
    double mean_c = 0.0;
    double mean_lambda_given_pleio = 0.0;
    double mean_q_given_pleio = 0.0;
    double mean_indirect_given_two_path = 0.0;
    double log_evidence = 0.0;
};

const std::array<std::string, 8>& joint_graph_state_names();

JointGraphResult fit_joint_graph(
    const std::vector<JointGraphObservation>& observations,
    const JointGraphOptions& options = JointGraphOptions());

std::vector<JointGraphObservation> read_joint_graph_tsv(const std::string& path);
void write_joint_graph_result_tsv(const JointGraphResult& result,
                                  const std::string& path);

}  // namespace bmediator

#endif
