#include "joint_graph_v02.h"

#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

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
        else if (key == "max_evidence_discrepancy") {
            result.max_evidence_discrepancy = value;
        }
        else if (key == "max_quadrature_discrepancy") {
            result.max_quadrature_discrepancy = value;
        }
        else if (key == "quadrature_escalation_threshold") {
            result.quadrature_escalation_threshold = value;
        }
        else if (key == "min_role_blocks") {
            if (value != std::floor(value)) {
                throw std::runtime_error("min_role_blocks must be an integer");
            }
            result.min_role_blocks = static_cast<int>(value);
        }
        else if (key == "optimizer_iterations") {
            if (value != std::floor(value)) {
                throw std::runtime_error("optimizer_iterations must be an integer");
            }
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
    const bool integrated_q_mode =
        argc == 9 && std::string(argv[1]) == "--loglik-integrated-q";
    const std::string program = argc > 0 ? argv[0] : "bmediator-joint";
    const auto print_usage = [&]() {
        std::cerr << "usage: " << program
                  << " --input INPUT.tsv --ld LD.tsv --out OUTPUT.tsv"
                  << " [--options OPTIONS.tsv]\n"
                  << "       " << program
                  << " INPUT.tsv LD.tsv OUTPUT.tsv [OPTIONS.tsv]\n"
                  << "       " << program << " --loglik INPUT.tsv LD.tsv"
                  << " A B C LAMBDA Q ETA\n"
                  << "       " << program
                  << " --loglik-integrated-q INPUT.tsv LD.tsv"
                  << " A B C LAMBDA ETA\n";
    };
    if (argc == 2 && std::string(argv[1]) == "--help") {
        print_usage();
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--version") {
        std::cout << "BMEDIATOR joint model JG-0.2.3\n";
        return 0;
    }
    if (argc < 2) {
        print_usage();
        return 2;
    }
    try {
        std::string input_path;
        std::string ld_path;
        std::string output_path;
        std::string options_path;
        if (!likelihood_mode && !integrated_q_mode &&
            std::string(argv[1]).rfind("--", 0) == 0) {
            std::unordered_map<std::string, std::string> named;
            for (int i = 1; i < argc; i += 2) {
                if (i + 1 >= argc) throw std::runtime_error("missing option value");
                const std::string key = argv[i];
                if (key != "--input" && key != "--ld" && key != "--out" &&
                    key != "--options") {
                    throw std::runtime_error("unknown command option: " + key);
                }
                if (named.count(key)) {
                    throw std::runtime_error("duplicate command option: " + key);
                }
                named[key] = argv[i + 1];
            }
            if (!named.count("--input") || !named.count("--ld") ||
                !named.count("--out")) {
                throw std::runtime_error("--input, --ld, and --out are required");
            }
            input_path = named["--input"];
            ld_path = named["--ld"];
            output_path = named["--out"];
            if (named.count("--options")) options_path = named["--options"];
        } else if (!likelihood_mode && !integrated_q_mode) {
            if (argc != 4 && argc != 5) {
                print_usage();
                return 2;
            }
            input_path = argv[1];
            ld_path = argv[2];
            output_path = argv[3];
            if (argc == 5) options_path = argv[4];
        }
        if (likelihood_mode || integrated_q_mode) {
            input_path = argv[2];
            ld_path = argv[3];
        }
        const auto observations = bmediator::read_joint_graph_v02_tsv(input_path);
        const auto ld = bmediator::read_joint_graph_v02_ld(ld_path, observations);
        if (likelihood_mode) {
            const double value = bmediator::joint_graph_v02_log_likelihood(
                observations, ld, std::stod(argv[4]), std::stod(argv[5]),
                std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]),
                std::stod(argv[9]));
            std::cout << std::setprecision(17) << value << '\n';
            return 0;
        }
        if (integrated_q_mode) {
            const double value =
                bmediator::joint_graph_v02_log_likelihood_integrated_q(
                    observations, ld, std::stod(argv[4]), std::stod(argv[5]),
                    std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]));
            std::cout << std::setprecision(17) << value << '\n';
            return 0;
        }
        const auto options = !options_path.empty()
            ? read_options(options_path) : bmediator::JointGraphV02Options();
        const auto result = bmediator::fit_joint_graph_v02(observations, ld, options);
        bmediator::write_joint_graph_v02_result_tsv(result, output_path);
    } catch (const std::exception& error) {
        std::cerr << "joint_graph_v02_cli: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
