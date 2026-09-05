#include "joint_graph_v02.h"

#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

bmediator::JointGraphV02Options read_options(const std::string& path) {
    bmediator::JointGraphV02Options result;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open options file: " + path);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        const size_t tab = line.find('\t');
        if (tab == std::string::npos) {
            throw std::runtime_error("options must be tab-delimited key/value pairs");
        }
        const std::string key = line.substr(0, tab);
        const double value = std::stod(line.substr(tab + 1));
        if (key == "pi_xm") result.pi_xm = value;
        else if (key == "pi_my") result.pi_my = value;
        else if (key == "pi_sparse") result.pi_sparse = value;
        else if (key == "pi_directional") result.pi_directional = value;
        else if (key == "prior_sd_a") result.prior_sd_a = value;
        else if (key == "prior_sd_b") result.prior_sd_b = value;
        else if (key == "prior_sd_c") result.prior_sd_c = value;
        else if (key == "prior_sd_lambda") result.prior_sd_lambda = value;
        else if (key == "prior_sd_eta") result.prior_sd_eta = value;
        else if (key == "q_alpha") result.q_alpha = value;
        else if (key == "q_beta") result.q_beta = value;
        else if (key == "max_cross_block_ld") result.max_cross_block_ld = value;
        else if (key == "optimizer_iterations") {
            result.optimizer_iterations = static_cast<int>(value);
        } else if (key == "optimizer_tolerance") {
            result.optimizer_tolerance = value;
        } else {
            throw std::runtime_error("unknown option: " + key);
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    const bool likelihood_mode = argc == 10 && std::string(argv[1]) == "--loglik";
    if ((argc != 4 && argc != 5) && !likelihood_mode) {
        std::cerr << "usage: joint_graph_v02_cli INPUT.tsv LD.tsv OUTPUT.tsv [OPTIONS.tsv]\n"
                  << "       joint_graph_v02_cli --loglik INPUT.tsv LD.tsv"
                  << " A B C LAMBDA Q ETA\n";
        return 2;
    }
    try {
        const int offset = likelihood_mode ? 1 : 0;
        const auto observations = bmediator::read_joint_graph_v02_tsv(argv[1 + offset]);
        const auto ld = bmediator::read_joint_graph_v02_ld(argv[2 + offset], observations);
        if (likelihood_mode) {
            const double value = bmediator::joint_graph_v02_log_likelihood(
                observations, ld, std::stod(argv[4]), std::stod(argv[5]),
                std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]),
                std::stod(argv[9]));
            std::cout << std::setprecision(17) << value << '\n';
            return 0;
        }
        const auto options = argc == 5
            ? read_options(argv[4]) : bmediator::JointGraphV02Options();
        const auto result = bmediator::fit_joint_graph_v02(observations, ld, options);
        bmediator::write_joint_graph_v02_result_tsv(result, argv[3]);
    } catch (const std::exception& error) {
        std::cerr << "joint_graph_v02_cli: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
