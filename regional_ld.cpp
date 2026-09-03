#include "regional_ld.h"

#include <unordered_set>

namespace bmediator {

namespace {

struct FineMapSignal {
    int lead = -1;
    std::vector<double> log_bf;
    std::vector<int> credible_set;
};

struct ConditionalScores {
    std::vector<double> z;
    std::vector<double> residual_variance;
};

double log_sum_exp(const std::vector<double>& values) {
    if (values.empty()) return -std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (double value : values) {
        if (std::isfinite(value)) maximum = std::max(maximum, value);
    }
    if (!std::isfinite(maximum)) return maximum;
    double total = 0.0;
    for (double value : values) {
        if (std::isfinite(value)) total += std::exp(value - maximum);
    }
    return maximum + std::log(total);
}

double wakefield_log_bf_from_z(double z, double se, double prior_variance) {
    if (!std::isfinite(z) || !std::isfinite(se) || se <= 0.0 ||
        !std::isfinite(prior_variance) || prior_variance <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    const double variance = se * se;
    const double shrinkage = prior_variance / (variance + prior_variance);
    return 0.5 * (std::log1p(-shrinkage) + shrinkage * z * z);
}

bool invert_matrix(std::vector<std::vector<double>> matrix,
                   std::vector<std::vector<double>>& inverse) {
    const int n = static_cast<int>(matrix.size());
    inverse.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        matrix[i][i] += 1e-8;
        inverse[i][i] = 1.0;
    }
    for (int column = 0; column < n; ++column) {
        int pivot = column;
        for (int row = column + 1; row < n; ++row) {
            if (std::fabs(matrix[row][column]) > std::fabs(matrix[pivot][column])) {
                pivot = row;
            }
        }
        if (std::fabs(matrix[pivot][column]) < 1e-10) return false;
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            std::swap(inverse[pivot], inverse[column]);
        }
        const double scale = matrix[column][column];
        for (int j = 0; j < n; ++j) {
            matrix[column][j] /= scale;
            inverse[column][j] /= scale;
        }
        for (int row = 0; row < n; ++row) {
            if (row == column) continue;
            const double factor = matrix[row][column];
            for (int j = 0; j < n; ++j) {
                matrix[row][j] -= factor * matrix[column][j];
                inverse[row][j] -= factor * inverse[column][j];
            }
        }
    }
    return true;
}

class RegionalLD {
public:
    RegionalLD(const PlinkData& plink, const std::vector<int>& indices) {
        values_ = plink.extract_genotypes_double(indices);
        n_samples_ = plink.n_samples;
        for (auto& variant : values_) {
            double mean = std::accumulate(variant.begin(), variant.end(), 0.0) /
                          static_cast<double>(n_samples_);
            double sum_squares = 0.0;
            for (double value : variant) {
                const double centered = value - mean;
                sum_squares += centered * centered;
            }
            const double sd = std::sqrt(sum_squares / static_cast<double>(n_samples_));
            if (sd <= 1e-10) {
                std::fill(variant.begin(), variant.end(), 0.0);
            } else {
                for (double& value : variant) value = (value - mean) / sd;
            }
        }
    }

    void cache_rows(const std::vector<int>& indices) const {
        for (int index : indices) {
            if (row_cache_.count(index)) continue;
            std::vector<double> row(values_.size(), 0.0);
            for (int other = 0; other < static_cast<int>(values_.size()); ++other) {
                row[other] = direct_correlation(index, other);
            }
            row_cache_[index] = std::move(row);
        }
    }

    double correlation(int left, int right) const {
        auto left_row = row_cache_.find(left);
        if (left_row != row_cache_.end()) return left_row->second[right];
        auto right_row = row_cache_.find(right);
        if (right_row != row_cache_.end()) return right_row->second[left];
        return direct_correlation(left, right);
    }

private:
    double direct_correlation(int left, int right) const {
        if (left == right) return 1.0;
        double value = 0.0;
        for (int sample = 0; sample < n_samples_; ++sample) {
            value += values_[left][sample] * values_[right][sample];
        }
        value /= static_cast<double>(n_samples_);
        return std::max(-1.0, std::min(1.0, value));
    }

    std::vector<std::vector<double>> values_;
    int n_samples_ = 0;
    mutable std::map<int, std::vector<double>> row_cache_;
};

ConditionalScores conditional_z_scores(const std::vector<double>& z,
                                       const std::vector<int>& conditioning,
                                       const RegionalLD& ld) {
    if (conditioning.empty()) return {z, std::vector<double>(z.size(), 1.0)};
    ld.cache_rows(conditioning);
    const int k = static_cast<int>(conditioning.size());
    std::vector<std::vector<double>> r_ss(k, std::vector<double>(k));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            r_ss[i][j] = ld.correlation(conditioning[i], conditioning[j]);
        }
    }
    std::vector<std::vector<double>> inverse;
    ConditionalScores result{
        std::vector<double>(z.size(), std::numeric_limits<double>::quiet_NaN()),
        std::vector<double>(z.size(), std::numeric_limits<double>::quiet_NaN())
    };
    if (!invert_matrix(r_ss, inverse)) return result;

    std::vector<double> inverse_z(k, 0.0);
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            inverse_z[i] += inverse[i][j] * z[conditioning[j]];
        }
    }
    std::unordered_set<int> conditioned(conditioning.begin(), conditioning.end());
    for (int variant = 0; variant < static_cast<int>(z.size()); ++variant) {
        if (conditioned.count(variant)) continue;
        std::vector<double> r(k, 0.0);
        for (int i = 0; i < k; ++i) {
            r[i] = ld.correlation(variant, conditioning[i]);
        }
        double conditional_mean = 0.0;
        double explained_variance = 0.0;
        for (int i = 0; i < k; ++i) {
            conditional_mean += r[i] * inverse_z[i];
            for (int j = 0; j < k; ++j) {
                explained_variance += r[i] * inverse[i][j] * r[j];
            }
        }
        const double residual_variance = 1.0 - explained_variance;
        if (residual_variance <= 1e-4) continue;
        result.residual_variance[variant] = residual_variance;
        result.z[variant] = (z[variant] - conditional_mean) /
                            std::sqrt(residual_variance);
    }
    return result;
}

std::vector<int> select_independent_signals(const std::vector<double>& z,
                                            const RegionalLD& ld,
                                            const Options& opts) {
    std::vector<int> selected;
    while (static_cast<int>(selected.size()) < opts.regional_max_signals) {
        ConditionalScores conditional = conditional_z_scores(z, selected, ld);
        int best = -1;
        double best_abs_z = 0.0;
        for (int variant = 0; variant < static_cast<int>(conditional.z.size()); ++variant) {
            if (!std::isfinite(conditional.z[variant])) continue;
            const double abs_z = std::fabs(conditional.z[variant]);
            if (abs_z > best_abs_z) {
                best_abs_z = abs_z;
                best = variant;
            }
        }
        if (best < 0) break;
        const double p = std::erfc(best_abs_z / std::sqrt(2.0));
        if (!std::isfinite(p) || p > opts.regional_signal_p) break;
        selected.push_back(best);
    }
    return selected;
}

std::vector<FineMapSignal> fine_map_signals(const std::vector<double>& beta,
                                            const std::vector<double>& se,
                                            const RegionalLD& ld,
                                            double prior_variance,
                                            const Options& opts) {
    std::vector<double> z(beta.size());
    for (size_t i = 0; i < beta.size(); ++i) z[i] = beta[i] / se[i];
    const std::vector<int> selected = select_independent_signals(z, ld, opts);
    std::vector<FineMapSignal> signals;
    for (size_t signal_index = 0; signal_index < selected.size(); ++signal_index) {
        std::vector<int> other_signals;
        for (size_t other = 0; other < selected.size(); ++other) {
            if (other != signal_index) other_signals.push_back(selected[other]);
        }
        const ConditionalScores conditional = conditional_z_scores(z, other_signals, ld);
        FineMapSignal signal;
        signal.lead = selected[signal_index];
        signal.log_bf.resize(beta.size());
        for (size_t variant = 0; variant < beta.size(); ++variant) {
            signal.log_bf[variant] = wakefield_log_bf_from_z(
                conditional.z[variant],
                se[variant] / std::sqrt(conditional.residual_variance[variant]),
                prior_variance
            );
        }

        const double normalizer = log_sum_exp(signal.log_bf);
        std::vector<std::pair<double, int>> posterior;
        posterior.reserve(beta.size());
        for (int variant = 0; variant < static_cast<int>(beta.size()); ++variant) {
            if (!std::isfinite(signal.log_bf[variant])) continue;
            posterior.push_back({std::exp(signal.log_bf[variant] - normalizer), variant});
        }
        std::sort(posterior.begin(), posterior.end(),
                  [](const auto& left, const auto& right) { return left.first > right.first; });
        double coverage = 0.0;
        for (const auto& item : posterior) {
            signal.credible_set.push_back(item.second);
            coverage += item.first;
            if (coverage >= opts.regional_coverage) break;
        }
        signals.push_back(std::move(signal));
    }
    return signals;
}

RegionalSignalPairResult coloc_signal_pair(const FineMapSignal& protein_signal,
                                           const FineMapSignal& outcome_signal,
                                           int protein_index,
                                           int outcome_index,
                                           const ProteinData& protein,
                                           const RegionalLD& ld,
                                           const Options& opts) {
    RegionalSignalPairResult result;
    result.protein_signal = protein_index + 1;
    result.outcome_signal = outcome_index + 1;
    result.protein_lead = protein.regional_cis_rsid[protein_signal.lead];
    result.outcome_lead = protein.regional_cis_rsid[outcome_signal.lead];

    const double log_sum_protein = log_sum_exp(protein_signal.log_bf);
    const double log_sum_outcome = log_sum_exp(outcome_signal.log_bf);
    std::vector<double> shared_terms;
    shared_terms.reserve(protein_signal.log_bf.size());
    for (size_t i = 0; i < protein_signal.log_bf.size(); ++i) {
        shared_terms.push_back(protein_signal.log_bf[i] + outcome_signal.log_bf[i]);
    }
    const double log_sum_shared = log_sum_exp(shared_terms);
    const double log_all_pairs = log_sum_protein + log_sum_outcome;
    const double shared_fraction = std::exp(std::min(0.0, log_sum_shared - log_all_pairs));
    const double log_sum_distinct = shared_fraction >= 1.0
        ? -std::numeric_limits<double>::infinity()
        : log_all_pairs + std::log1p(-shared_fraction);

    const std::vector<double> hypotheses = {
        0.0,
        std::log(opts.regional_prior_pp) + log_sum_protein,
        std::log(opts.regional_prior_outcome) + log_sum_outcome,
        std::log(opts.regional_prior_pp) + std::log(opts.regional_prior_outcome) +
            log_sum_distinct,
        std::log(opts.regional_prior_shared) + log_sum_shared
    };
    const double normalizer = log_sum_exp(hypotheses);
    result.pp_h0 = std::exp(hypotheses[0] - normalizer);
    result.pp_h1 = std::exp(hypotheses[1] - normalizer);
    result.pp_h2 = std::exp(hypotheses[2] - normalizer);
    result.pp_h3 = std::exp(hypotheses[3] - normalizer);
    result.pp_h4 = std::exp(hypotheses[4] - normalizer);
    const double both = result.pp_h3 + result.pp_h4;
    result.shared_given_both = both > 0.0 ? result.pp_h4 / both : 0.0;
    const double lead_r = ld.correlation(protein_signal.lead, outcome_signal.lead);
    result.lead_pair_r2 = lead_r * lead_r;
    for (int left : protein_signal.credible_set) {
        for (int right : outcome_signal.credible_set) {
            const double r = ld.correlation(left, right);
            result.max_credible_set_pair_r2 = std::max(
                result.max_credible_set_pair_r2, r * r
            );
        }
    }

    if (both >= opts.regional_min_both &&
        result.shared_given_both >= opts.regional_min_shared) {
        result.interpretation = "SHARED_SIGNAL_SUPPORTED";
    } else if (both >= opts.regional_min_both &&
               result.shared_given_both <= 1.0 - opts.regional_min_shared) {
        result.interpretation = result.max_credible_set_pair_r2 >= opts.regional_high_ld_r2
            ? "DISTINCT_SIGNALS_HIGH_LD"
            : "DISTINCT_SIGNALS_LOW_MODERATE_LD";
    } else {
        result.interpretation = "SIGNAL_PAIR_AMBIGUOUS";
    }
    return result;
}

} // namespace

void compute_multisignal_regional_evidence(ProteinData& protein,
                                           const PlinkData& plink,
                                           const Options& opts) {
    protein.regional_method = "ld-multisignal";
    protein.regional_multisignal_evaluated = true;
    protein.regional_signal_pairs.clear();
    if (!protein.regional_data_complete || protein.regional_bim_index.size() < 2) {
        protein.regional_multisignal_interpretation = "NO_REGIONAL_DATA";
        return;
    }

    RegionalLD ld(plink, protein.regional_bim_index);
    const auto protein_signals = fine_map_signals(
        protein.regional_pp_beta, protein.regional_pp_se, ld,
        opts.regional_prior_var_pp, opts
    );
    const auto outcome_signals = fine_map_signals(
        protein.regional_outcome_beta, protein.regional_outcome_se, ld,
        opts.regional_prior_var_outcome, opts
    );
    protein.regional_protein_signals = static_cast<int>(protein_signals.size());
    protein.regional_outcome_signals = static_cast<int>(outcome_signals.size());

    if (protein_signals.empty()) {
        protein.regional_multisignal_interpretation = "NO_PROTEIN_CREDIBLE_SET";
        return;
    }
    if (outcome_signals.empty()) {
        protein.regional_multisignal_interpretation = "NO_OUTCOME_CREDIBLE_SET";
        return;
    }

    for (int i = 0; i < static_cast<int>(protein_signals.size()); ++i) {
        for (int j = 0; j < static_cast<int>(outcome_signals.size()); ++j) {
            protein.regional_signal_pairs.push_back(coloc_signal_pair(
                protein_signals[i], outcome_signals[j], i, j, protein, ld, opts
            ));
        }
    }

    const RegionalSignalPairResult* chosen = nullptr;
    for (const auto& pair : protein.regional_signal_pairs) {
        if (pair.interpretation != "SHARED_SIGNAL_SUPPORTED") continue;
        if (chosen == nullptr || chosen->interpretation != "SHARED_SIGNAL_SUPPORTED" ||
            pair.pp_h4 > chosen->pp_h4) {
            chosen = &pair;
        }
    }
    if (chosen == nullptr) {
        for (const auto& pair : protein.regional_signal_pairs) {
            const bool is_distinct = pair.interpretation == "DISTINCT_SIGNALS_HIGH_LD" ||
                                     pair.interpretation == "DISTINCT_SIGNALS_LOW_MODERATE_LD";
            if (!is_distinct) continue;
            const bool chosen_distinct = chosen != nullptr &&
                (chosen->interpretation == "DISTINCT_SIGNALS_HIGH_LD" ||
                 chosen->interpretation == "DISTINCT_SIGNALS_LOW_MODERATE_LD");
            if (!chosen_distinct || pair.pp_h3 > chosen->pp_h3) chosen = &pair;
        }
    }
    if (chosen == nullptr) {
        for (const auto& pair : protein.regional_signal_pairs) {
            if (chosen == nullptr || pair.pp_h3 + pair.pp_h4 > chosen->pp_h3 + chosen->pp_h4) {
                chosen = &pair;
            }
        }
    }

    if (chosen != nullptr) {
        protein.regional_best_pp_shared = chosen->pp_h4;
        protein.regional_best_pp_distinct = chosen->pp_h3;
        protein.regional_best_shared_given_both = chosen->shared_given_both;
        protein.regional_best_cs_pair_r2 = chosen->max_credible_set_pair_r2;
        protein.regional_multisignal_interpretation = chosen->interpretation;
    } else {
        protein.regional_multisignal_interpretation = "SIGNAL_PAIRS_AMBIGUOUS";
    }
}

} // namespace bmediator
