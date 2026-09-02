#include "../plink_ld.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Row {
    std::string snp;
    std::string a1;
    std::string a2;
    double freq = 0.0;
    double beta = 0.0;
    double se = 0.0;
    double p = 1.0;
    int chr = 0;
    int bp = 0;
    int bim_idx = -1;
};

int main(int argc, char* argv[]) {
    if (argc != 7) {
        std::cerr << "Usage: export_clumped_sumstats <sumstats.tsv> <plink_prefix> <p_thresh> <r2_thresh> <window_kb> <out.tsv>\n";
        return 1;
    }

    const std::string sumstat_path = argv[1];
    const std::string plink_prefix = argv[2];
    const double p_thresh = std::stod(argv[3]);
    const double r2_thresh = std::stod(argv[4]);
    const int window_kb = std::stoi(argv[5]);
    const std::string out_path = argv[6];

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
    std::getline(fin, line);

    std::vector<Row> rows;
    std::vector<int> bim_indices;
    std::vector<double> pvalues;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Row row;
        if (!(iss >> row.snp >> row.a1 >> row.a2 >> row.freq >> row.beta >> row.se >> row.p >> row.chr >> row.bp)) {
            continue;
        }
        if (row.p >= p_thresh) continue;
        auto it = plink.rsid_to_idx.find(row.snp);
        if (it == plink.rsid_to_idx.end()) continue;
        row.bim_idx = it->second;
        rows.push_back(row);
        bim_indices.push_back(row.bim_idx);
        pvalues.push_back(row.p);
    }

    auto clumped = bmediator::ld_clump(plink, bim_indices, pvalues, r2_thresh, window_kb);

    std::ofstream fout(out_path);
    if (!fout.is_open()) {
        std::cerr << "failed to open " << out_path << "\n";
        return 1;
    }

    fout << "SNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tCHR\tBP\n";
    for (int idx : clumped.index_snps) {
        const auto& row = rows[idx];
        fout << row.snp << '\t' << row.a1 << '\t' << row.a2 << '\t'
             << row.freq << '\t' << row.beta << '\t' << row.se << '\t'
             << row.p << '\t' << row.chr << '\t' << row.bp << '\n';
    }

    std::cout << "input_candidates\t" << rows.size() << "\n";
    std::cout << "clumped_instruments\t" << clumped.index_snps.size() << "\n";
    std::cout << "output_file\t" << out_path << "\n";
    return 0;
}
