#include "joint_graph_model.h"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: joint_graph_cli INPUT.tsv OUTPUT.tsv\n";
        return 2;
    }
    try {
        const auto observations = bmediator::read_joint_graph_tsv(argv[1]);
        const auto result = bmediator::fit_joint_graph(observations);
        bmediator::write_joint_graph_result_tsv(result, argv[2]);
    } catch (const std::exception& error) {
        std::cerr << "joint_graph_cli: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
