#include "../plink_ld.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: count_clumped_bmi <bmi_sumstats.tsv> <plink_prefix> <p_thresh> <r2_thresh> <window_kb>\n";
        return 1;
    }

    const std::string sumstat_path = argv[1];
    const std::string plink_prefix = argv[2];
    const double p_thresh = std::stod(argv[3]);
    const double r2_thresh = std::stod(argv[4]);
    const int window_kb = std::stoi(argv[5]);

    bmediator::PlinkData plink;
    if (!plink.load(plink_prefix)) {
        std::cerr << "failed to load PLINK reference\n";
        return 1;
    }

    std::ifstream fin(sumstat_path);
    if (!fin.is_open()) {
        std::cerr << "failed to open " << sumstat_path << "\n";
        return 1;
    }

    std::string line;
    std::getline(fin, line); // header

    std::vector<int> bim_indices;
    std::vector<double> pvalues;
    std::vector<std::string> rsids;
    long long total = 0;
    long long passing = 0;
    long long in_ref = 0;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        ++total;
        std::istringstream iss(line);
        std::string snp, a1, a2;
        double freq, beta, se, p;
        int chr, bp;
        if (!(iss >> snp >> a1 >> a2 >> freq >> beta >> se >> p >> chr >> bp)) {
            continue;
        }
        if (p >= p_thresh) continue;
        ++passing;
        auto it = plink.rsid_to_idx.find(snp);
        if (it == plink.rsid_to_idx.end()) continue;
        ++in_ref;
        bim_indices.push_back(it->second);
        pvalues.push_back(p);
        rsids.push_back(snp);
    }

    auto clumped = bmediator::ld_clump(plink, bim_indices, pvalues, r2_thresh, window_kb);

    std::cout << "total_rows\t" << total << "\n";
    std::cout << "passing_p_threshold\t" << passing << "\n";
    std::cout << "passing_and_in_reference\t" << in_ref << "\n";
    std::cout << "clumped_instruments\t" << clumped.index_snps.size() << "\n";
    return 0;
}
