#include "bmediator.h"

using namespace bmediator;

int main(int argc, char* argv[]) {
    auto t_start = std::chrono::steady_clock::now();
    Options opts;
    parse_args(argc, argv, opts);

    if (!opts.protein_gwas_list_file.empty()) {
        run_full_pipeline(opts);
    } else {
        run_legacy_pipeline(opts);
    }

    auto t_end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "\nAnalysis completed in " << std::fixed << std::setprecision(1)
              << elapsed << " seconds.\n";

    return 0;
}
