#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/bmediator"
OUT_DIR="$ROOT/build/test"
OUT_PREFIX="$OUT_DIR/test_run"

if [[ ! -x "$BIN" ]]; then
  echo "error: binary not found at $BIN" >&2
  echo "build first with: make" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
rm -f "${OUT_PREFIX}.mediation" "${OUT_PREFIX}.hyp" \
  "${OUT_PREFIX}.instruments" "${OUT_PREFIX}.regional"

"$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --out "$OUT_PREFIX" \
  --p-thresh-rf 5e-6 \
  --p-thresh-cis 5e-4 \
  --max-eb-iter 3 \
  --max-cavi-iter 50

if [[ ! -f "${OUT_PREFIX}.mediation" ]]; then
  echo "error: mediation output was not created" >&2
  exit 1
fi

if [[ ! -f "${OUT_PREFIX}.hyp" ]]; then
  echo "error: hyperparameter output was not created" >&2
  exit 1
fi

if ! awk -F '\t' '$1=="regional_prior_shared" {found=1; exit !($2==0.00000001)} END {if (!found) exit 1}' "${OUT_PREFIX}.hyp"; then
  echo "error: analytical independence prior was not preserved in hyperparameter output" >&2
  exit 1
fi

if ! awk -F '\t' 'NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next} NR==2 {exit !( $1=="P1" && $3==1 && $4==1 && $5==1 && $(h["P_M1"])>0.5 )}' "${OUT_PREFIX}.mediation"; then
  echo "error: smoke test output did not match expected structure" >&2
  exit 1
fi

if ! awk -F '\t' 'NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next} NR==2 {exit !( $(h["P_mediator_ld_resolved"])==0 && $(h["P_mediator_identified"])==$(h["P_mediator_ld_resolved"]) && $(h["mediation_identifiability"])=="UNRESOLVED_NO_LD_REFERENCE" )}' "${OUT_PREFIX}.mediation"; then
  echo "error: legacy mode must not identify mediation without an LD reference" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --out "$OUT_DIR/invalid_priors" \
  --prior-p0 0 --prior-p1 0 --prior-p2 0 \
  --prior-p3 0 --prior-p4 0 --prior-p5 0 \
  >/dev/null 2>&1; then
  echo "error: invalid zero-sum scenario priors were accepted" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --out "$OUT_DIR/invalid_clump_r2" \
  --clump-r2 0 \
  >/dev/null 2>&1; then
  echo "error: invalid zero LD clumping threshold was accepted" >&2
  exit 1
fi

FULL_DIR="$OUT_DIR/full_mode_fixture"
FULL_PREFIX="$OUT_DIR/full_mode_run"
FULL_FACTOR_PREFIX="$OUT_DIR/full_mode_factorized_run"
FULL_FACTOR_SINGLE_SIGNAL_PREFIX="$OUT_DIR/full_mode_factorized_single_signal_run"
FULL_FACTOR_OVERLAP_PREFIX="$OUT_DIR/full_mode_factorized_overlap_run"
FULL_SINGLE_PREFIX="$OUT_DIR/full_mode_single_run"
FULL_LEGACY_REGIONAL_PREFIX="$OUT_DIR/full_mode_legacy_regional_run"
rm -rf "$FULL_DIR"
python3 "$ROOT/test/generate_full_mode_fixture.py" "$FULL_DIR"

"$BIN" \
  --rf-sumstat "$FULL_DIR/rf.txt" \
  --protein-gwas-list "$FULL_DIR/manifest.txt" \
  --cancer-sumstat "$FULL_DIR/outcome.txt" \
  --protein-info "$FULL_DIR/protein_info.txt" \
  --bfile "$FULL_DIR/ldref" \
  --out "$FULL_PREFIX" \
  --p-thresh-rf 5e-6 \
  --p-thresh-cis 1e-4 \
  --cis-window 50 \
  --clump-r2 0.05 \
  --min-instruments 1 \
  --heidi-off \
  --no-steiger \
  --max-cavi-iter 80

if ! awk -F '\t' '$1=="clump_r2" {found=1; exit !($2==0.05)} END {if (!found) exit 1}' "${FULL_PREFIX}.hyp"; then
  echo "error: configured LD clumping threshold was not preserved in hyperparameter output" >&2
  exit 1
fi

if ! awk -F '\t' '$1=="regional_method" {found=1; exit !($2=="ld-multisignal")} END {if (!found) exit 1}' "${FULL_PREFIX}.hyp"; then
  echo "error: LD-aware multi-signal regional method is not the full-mode default" >&2
  exit 1
fi

if [[ ! -f "${FULL_PREFIX}.regional" ]]; then
  echo "error: regional signal-pair output was not created" >&2
  exit 1
fi

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  $1=="P_SHARED" {
    if ($(h["protein_lead"]) == $(h["outcome_lead"]) &&
        $(h["interpretation"]) == "SHARED_SIGNAL_SUPPORTED") diagonal_shared++
    if ($(h["protein_lead"]) != $(h["outcome_lead"]) &&
        $(h["interpretation"]) ~ /^DISTINCT_SIGNALS_/) off_diagonal_distinct++
  }
  $1=="P_DISTINCT" && $(h["interpretation"]) ~ /^DISTINCT_SIGNALS_/ {distinct=1}
  END {exit !(diagonal_shared == 2 && off_diagonal_distinct == 2 && distinct)}
' "${FULL_PREFIX}.regional"; then
  echo "error: signal-pair output did not resolve the shared and distinct fixtures" >&2
  exit 1
fi

"$BIN" \
  --rf-sumstat "$FULL_DIR/rf.txt" \
  --protein-gwas-list "$FULL_DIR/manifest.txt" \
  --cancer-sumstat "$FULL_DIR/outcome.txt" \
  --protein-info "$FULL_DIR/protein_info.txt" \
  --bfile "$FULL_DIR/ldref" \
  --out "$FULL_FACTOR_PREFIX" \
  --structural-method factorized \
  --factor-independent-selection \
  --p-thresh-rf 5e-6 \
  --p-thresh-cis 1e-4 \
  --cis-window 50 \
  --clump-r2 0.05 \
  --min-instruments 1 \
  --heidi-off \
  --no-steiger \
  --factor-pleio-sd-xm 0.01 \
  --factor-pleio-sd-my 0.01 \
  --factor-pleio-sd-xy 0.01 \
  --max-cavi-iter 80

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  $1=="P_FACTOR" {
    found=1
    if ($(h["factor_nA"]) != 4 || $(h["factor_nB"]) != 4 ||
        $(h["factor_log_BF_XM"]) < 2.302585 ||
        $(h["factor_log_BF_MY"]) < 2.302585 ||
        $(h["factor_log_e_p2e_mediation"]) == "nan" ||
        $(h["factor_e_q_p2e_EBH"]) > 0.05 ||
        $(h["factor_PP_XM"]) < 0.99 ||
        $(h["factor_PP_MY"]) < 0.99 ||
        $(h["factor_PP_two_stage"]) < 0.99 ||
        $(h["factor_log_BF_slope_only_XM"]) == "nan" ||
        $(h["factor_log_BF_directional_only_XM"]) == "nan" ||
        $(h["factor_log_BF_slope_directional_XM"]) == "nan" ||
        $(h["factor_log_BF_slope_only_MY"]) == "nan" ||
        $(h["factor_log_BF_directional_only_MY"]) == "nan" ||
        $(h["factor_log_BF_slope_directional_MY"]) == "nan" ||
        $(h["factor_directional_intercept_XM_se"]) <= 0 ||
        $(h["factor_directional_intercept_MY_se"]) <= 0 ||
        $(h["factor_cross_set_max_r2"]) > 0.05 ||
        $(h["factor_posterior_local_fdr"]) > 0.01 ||
        $(h["factor_posterior_cum_fdr"]) > 0.05 ||
        $(h["factor_posterior_rank"]) != 1 ||
        $(h["factor_log_e_XM"]) == "nan" ||
        $(h["factor_log_e_MY"]) == "nan" ||
        $(h["factor_log_e_mediation"]) == "nan" ||
        $(h["factor_ld_source"]) != "reference" ||
        $(h["factor_selection_design"]) != "independent-discovery" ||
        $(h["factor_effect_estimator"]) != "joint-directional-generalized-adjusted-profile-score" ||
        $(h["factor_beta1_ci_lower"]) > $(h["factor_beta1"]) ||
        $(h["factor_beta1_ci_upper"]) < $(h["factor_beta1"]) ||
        $(h["factor_beta2_ci_lower"]) > $(h["factor_beta2"]) ||
        $(h["factor_beta2_ci_upper"]) < $(h["factor_beta2"]) ||
        $(h["factor_indirect_ci_lower"]) > $(h["factor_indirect"]) ||
        $(h["factor_indirect_ci_upper"]) < $(h["factor_indirect"]) ||
        $(h["factor_two_stage_status"]) != "TWO_STAGE_EVIDENCE" ||
        $(h["factor_mediation_status"]) != "SUPPORTED_EXCLUSION_RESTRICTION_CONDITIONAL" ||
        $(h["factor_frequentist_status"]) != "SUPPORTED_SCALAR_DISPERSION_EXCLUSION_CONDITIONAL" ||
        $(h["factor_strict_status"]) != "SUPPORTED_GAUSSIAN_COVARIANCE_EXCLUSION_CONDITIONAL" ||
        $(h["factor_p2e_status"]) != "SUPPORTED_GAUSSIAN_COVARIANCE_EXCLUSION_CONDITIONAL" ||
        $(h["factor_posterior_status"]) != "SUPPORTED_BAYES_POSTERIOR_FDR_ASSUMPTION_CONDITIONAL" ||
        $(h["regional_independent_shared_signals"]) < 2 ||
        $(h["factor_log_BF_heterogeneity_MY"]) >= 2.302585) exit 1
  }
  END {if (!found) exit 1}
' "${FULL_FACTOR_PREFIX}.mediation"; then
  echo "error: factorized mode did not recover the deterministic mediation fixture" >&2
  exit 1
fi

python3 "$ROOT/test/check_student_t_evalue.py" "${FULL_FACTOR_PREFIX}.mediation"

if ! awk -F '\t' '
  $1=="factor_e_method" && $2=="student-t-density-ratio-mixture" {method=1}
  $1=="factor_e_shift_grid" && $2=="2,4,6" {shifts=1}
  $1=="factor_e_scale_grid" && $2=="2,4,8" {scales=1}
  END {exit !(method && shifts && scales)}
' "${FULL_FACTOR_PREFIX}.hyp"; then
  echo "error: factorized run did not record the frozen Student-t e-value specification" >&2
  exit 1
fi

if ! awk -F '\t' '
  $1=="P_FACTOR" && $2=="A" && $3=="q1" {bad=1}
  $1=="P_FACTOR" && $2=="B" && $3=="f6" {bad=1}
  $1=="P_FACTOR" && $2=="C" && $3=="f5" {cis_overlap=1}
  $1=="P_FACTOR" && $2=="C" && $3=="f6" {ld_overlap=1}
  END {exit bad || !(cis_overlap && ld_overlap)}
' "${FULL_FACTOR_PREFIX}.instruments"; then
  echo "error: factorized Set A/B/C partition did not remove cis or cross-LD overlap" >&2
  exit 1
fi

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  $1=="P_SHARED" {
    found=1
    if ($(h["factor_nA"]) != 2 || $(h["factor_nB"]) != 2 ||
        $(h["factor_two_stage_status"]) != "INSUFFICIENT_SET_A_AND_SET_B" ||
        $(h["factor_directional_collinearity_MY"]) < 0.98) exit 1
  }
  END {if (!found) exit 1}
' "${FULL_FACTOR_PREFIX}.mediation"; then
  echo "error: two-instrument slope/directional ambiguity was not reported" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --out "$OUT_DIR/invalid_factor_prior" \
  --factor-prior-my 1 \
  >/dev/null 2>&1; then
  echo "error: invalid factorized causal-leg prior was accepted" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --out "$OUT_DIR/invalid_factor_directional_prior" \
  --factor-prior-directional 1 \
  >/dev/null 2>&1; then
  echo "error: invalid directional-component prior was accepted" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --out "$OUT_DIR/invalid_factor_directional_variance" \
  --factor-directional-variance 0 \
  >/dev/null 2>&1; then
  echo "error: invalid directional-component variance was accepted" >&2
  exit 1
fi

if ! awk -F '\t' 'NR==1 {n=NF; next} NF!=n {exit 1}' "${FULL_FACTOR_PREFIX}.mediation"; then
  echo "error: mediation output rows do not match the declared columns" >&2
  exit 1
fi

if ! awk -F '\t' 'NR==1 {n=NF; next} NF!=n {exit 1}' "${FULL_FACTOR_PREFIX}.regional"; then
  echo "error: regional output rows do not match the declared columns" >&2
  exit 1
fi

"$BIN" \
  --rf-sumstat "$FULL_DIR/rf.txt" \
  --protein-gwas-list "$FULL_DIR/manifest_factor.txt" \
  --cancer-sumstat "$FULL_DIR/outcome.txt" \
  --protein-info "$FULL_DIR/protein_info.txt" \
  --bfile "$FULL_DIR/ldref" \
  --out "$FULL_FACTOR_SINGLE_SIGNAL_PREFIX" \
  --structural-method factorized \
  --factor-independent-selection \
  --regional-max-signals 1 \
  --p-thresh-rf 5e-6 \
  --p-thresh-cis 1e-4 \
  --cis-window 50 \
  --clump-r2 0.05 \
  --min-instruments 1 \
  --heidi-off \
  --no-steiger \
  --factor-pleio-sd-xm 0.01 \
  --factor-pleio-sd-my 0.01 \
  --factor-pleio-sd-xy 0.01 \
  --max-cavi-iter 80

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  $1=="P_FACTOR" {
    found=1
    if ($(h["factor_two_stage_status"]) != "TWO_STAGE_EVIDENCE" ||
        $(h["regional_independent_shared_signals"]) != 1 ||
        $(h["mediation_identifiability"]) != "UNRESOLVED_SINGLE_SHARED_SIGNAL" ||
        $(h["factor_mediation_status"]) != "UNRESOLVED_SINGLE_SHARED_SIGNAL") exit 1
  }
  END {if (!found) exit 1}
' "${FULL_FACTOR_SINGLE_SIGNAL_PREFIX}.mediation"; then
  echo "error: a single H4 signal was incorrectly promoted to identified mediation" >&2
  exit 1
fi

"$BIN" \
  --rf-sumstat "$FULL_DIR/rf.txt" \
  --protein-gwas-list "$FULL_DIR/manifest.txt" \
  --cancer-sumstat "$FULL_DIR/outcome.txt" \
  --protein-info "$FULL_DIR/protein_info.txt" \
  --bfile "$FULL_DIR/ldref" \
  --out "$FULL_FACTOR_OVERLAP_PREFIX" \
  --structural-method factorized \
  --factor-independent-selection \
  --p-thresh-rf 5e-6 \
  --p-thresh-cis 1e-4 \
  --cis-window 50 \
  --clump-r2 0.05 \
  --min-instruments 1 \
  --heidi-off \
  --no-steiger \
  --factor-pleio-sd-xm 0.01 \
  --factor-pleio-sd-my 0.01 \
  --factor-pleio-sd-xy 0.01 \
  --sampling-corr-pqtl-outcome 0.3 \
  --max-cavi-iter 80

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  $1=="P_FACTOR" {
    found=1
    if ($(h["factor_conjunction_q_BY"]) != "nan" ||
        $(h["factor_strict_conjunction_q_BY"]) != "nan" ||
        $(h["factor_e_q_EBH"]) != "nan" ||
        $(h["factor_e_q_p2e_EBH"]) != "nan" ||
        $(h["factor_frequentist_status"]) != "UNRESOLVED_SAMPLE_OVERLAP" ||
        $(h["factor_strict_status"]) != "UNRESOLVED_SAMPLE_OVERLAP" ||
        $(h["factor_p2e_status"]) != "UNRESOLVED_SAMPLE_OVERLAP" ||
        $(h["factor_ebh_status"]) != "UNRESOLVED_SAMPLE_OVERLAP") exit 1
  }
  END {if (!found) exit 1}
' "${FULL_FACTOR_OVERLAP_PREFIX}.mediation"; then
  echo "error: factorized overlap analysis exposed uncalibrated q-values" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --structural-method factorized \
  --factor-independent-selection \
  --out "$OUT_DIR/missing_selection_p" \
  >/dev/null 2>&1; then
  echo "error: independent selection accepted input without P_SELECT" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --out "$OUT_DIR/invalid_sampling_corr" \
  --sampling-corr-pqtl-outcome 1 \
  >/dev/null 2>&1; then
  echo "error: invalid sampling error correlation was accepted" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$ROOT/testdata/rf_sumstat.txt" \
  --pqtl-sumstat "$ROOT/testdata/pqtl_sumstat.txt" \
  --cancer-sumstat "$ROOT/testdata/cancer_sumstat.txt" \
  --protein-info "$ROOT/testdata/protein_info.txt" \
  --out "$OUT_DIR/invalid_factor_pleio_sd" \
  --factor-pleio-sd-my 0 \
  >/dev/null 2>&1; then
  echo "error: invalid factorized pleiotropy-prior SD was accepted" >&2
  exit 1
fi

"$BIN" \
  --rf-sumstat "$FULL_DIR/rf.txt" \
  --protein-gwas-list "$FULL_DIR/manifest.txt" \
  --cancer-sumstat "$FULL_DIR/outcome.txt" \
  --protein-info "$FULL_DIR/protein_info.txt" \
  --bfile "$FULL_DIR/ldref" \
  --out "$FULL_LEGACY_REGIONAL_PREFIX" \
  --regional-method single \
  --p-thresh-rf 5e-6 \
  --p-thresh-cis 1e-4 \
  --cis-window 50 \
  --clump-r2 0.05 \
  --min-instruments 1 \
  --heidi-off \
  --no-steiger \
  --max-cavi-iter 80

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  {if ($(h["regional_method"]) != "single" || $(h["regional_signal_pairs"]) != 0) exit 1}
  END {if (NR != 4) exit 1}
' "${FULL_LEGACY_REGIONAL_PREFIX}.mediation"; then
  echo "error: single-causal regional compatibility mode was not preserved" >&2
  exit 1
fi

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  $1=="P_SHARED" {
    found=1
    if ($(h["regional_shared_given_both"]) < 0.80 ||
        $(h["P_M1"]) <= $(h["P_M5"]) ||
        $(h["regional_method"]) != "ld-multisignal" ||
        $(h["regional_protein_signals"]) != 2 ||
        $(h["regional_outcome_signals"]) != 2 ||
        $(h["regional_signal_pairs"]) != 4 ||
        $(h["regional_independent_shared_signals"]) != 2 ||
        $(h["P_mediator_identified"]) != $(h["P_mediator_ld_resolved"]) ||
        $(h["mediation_identifiability"]) != "OVERIDENTIFIED_SHARED_SIGNALS_ASSUMPTION_CONDITIONAL") exit 1
  }
  END {if (!found) exit 1}
' "${FULL_PREFIX}.mediation"; then
  echo "error: shared-signal M1 fixture was not distinguished from M5 and conditionally identified" >&2
  exit 1
fi

if "$BIN" \
  --rf-sumstat "$FULL_DIR/rf.txt" \
  --protein-gwas-list "$FULL_DIR/manifest.txt" \
  --cancer-sumstat "$FULL_DIR/outcome.txt" \
  --protein-info "$FULL_DIR/protein_info.txt" \
  --bfile "$FULL_DIR/ldref" \
  --out "$OUT_DIR/invalid_regional_method" \
  --regional-method unsupported \
  >/dev/null 2>&1; then
  echo "error: invalid regional method was accepted" >&2
  exit 1
fi

# A regional model that gives M5 unrestricted protein/outcome effect
# covariance is observationally confounded with beta2 under M1. The withdrawn
# joint-ld implementation must not become accessible through the production CLI.
if "$BIN" \
  --rf-sumstat "$FULL_DIR/rf.txt" \
  --protein-gwas-list "$FULL_DIR/manifest.txt" \
  --cancer-sumstat "$FULL_DIR/outcome.txt" \
  --protein-info "$FULL_DIR/protein_info.txt" \
  --bfile "$FULL_DIR/ldref" \
  --out "$OUT_DIR/invalid_joint_ld" \
  --regional-method joint-ld \
  >/dev/null 2>&1; then
  echo "error: nonidentified joint-ld M1/M5 model was accepted" >&2
  exit 1
fi

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  $1=="P_DISTINCT" {
    found=1
    if ($(h["regional_shared_given_both"]) > 0.20 ||
        $(h["mediation_identifiability"]) != "LD_DISTINCT_SUPPORTED" ||
        $(h["P_mediator_ld_resolved"]) != 0) exit 1
  }
  END {if (!found) exit 1}
' "${FULL_PREFIX}.mediation"; then
  echo "error: distinct LD-confounded fixture was not rejected" >&2
  exit 1
fi

"$BIN" \
  --rf-sumstat "$FULL_DIR/rf.txt" \
  --protein-gwas-list "$FULL_DIR/manifest_single.txt" \
  --cancer-sumstat "$FULL_DIR/outcome.txt" \
  --protein-info "$FULL_DIR/protein_info.txt" \
  --bfile "$FULL_DIR/ldref" \
  --out "$FULL_SINGLE_PREFIX" \
  --p-thresh-rf 5e-6 \
  --p-thresh-cis 1e-4 \
  --cis-window 50 \
  --clump-r2 0.05 \
  --min-instruments 1 \
  --heidi-off \
  --no-steiger \
  --max-cavi-iter 80

if ! awk -F '\t' 'NR==2 {if ($1 != "P_SHARED") exit 1; found=1} NR>2 {exit 1} END {if (!found) exit 1}' \
  "${FULL_SINGLE_PREFIX}.mediation"; then
  echo "error: one-row manifest did not produce exactly one protein result" >&2
  exit 1
fi

M1_M5_PREFIX="$OUT_DIR/m1_m5_independent"
M1_M5_FACTOR_PREFIX="$OUT_DIR/m1_m5_factorized"
python3 "$ROOT/sim/run_m1_m5_task.py" \
  --config "$ROOT/sim/configs/m1_m5_identification_smoke.json" \
  --cell identified_setb2 \
  --replicate 1 \
  --outdir "$M1_M5_PREFIX" \
  --binary "$BIN"

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  {
    rows++
    if ($(h["true_scenario"]) != "M1" && $(h["true_scenario"]) != "M5") exit 1
    if ($(h["nA"]) > 6 || $(h["nB"]) > 2 || $(h["nC"]) != 0) exit 1
  }
  END {if (rows != 16) exit 1}
' "$M1_M5_PREFIX/classification/identified_setb2/rep_0001/task_metrics.tsv"; then
  echo "error: independent M1/M5 simulation runner mixed instrument sets across proteins" >&2
  exit 1
fi

python3 "$ROOT/sim/run_m1_m5_task.py" \
  --config "$ROOT/sim/configs/m1_m5_identification_smoke.json" \
  --cell identified_setb4 \
  --replicate 1 \
  --outdir "$M1_M5_FACTOR_PREFIX" \
  --binary "$BIN" \
  --structural-method factorized

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  {
    both_bf = $(h["factor_log_BF_XM"]) > 2.302585 &&
              $(h["factor_log_BF_MY"]) > 2.302585
    if ($(h["true_scenario"]) == "M1") {m1++; m1_beta2 += $(h["factor_beta2"]); m1_bf += both_bf}
    if ($(h["true_scenario"]) == "M5") {m5++; m5_beta2 += $(h["factor_beta2"]); m5_bf += both_bf}
  }
  END {
    if (m1 != 8 || m5 != 8 || m1_beta2/m1 < 0.20 ||
        (m5_beta2/m5 < -0.15 || m5_beta2/m5 > 0.15) ||
        m5_bf > 1) exit 1
  }
' "$M1_M5_FACTOR_PREFIX/classification/identified_setb4/rep_0001/task_metrics.tsv"; then
  echo "error: factorized engine did not separate the identified M1/M5 smoke fixture" >&2
  exit 1
fi

python3 "$ROOT/sim/build_benchmark_manifest.py" \
  --config "$ROOT/sim/configs/factorized_calibration_smoke.json" \
  --benchmark all \
  --out "$OUT_DIR/factorized_manifest.tsv"
python3 "$ROOT/sim/build_factorized_ld_manifest.py" \
  --replicates 1 \
  --out "$OUT_DIR/factorized_ld_manifest.tsv"
if LC_ALL=C grep -q $'\r' "$OUT_DIR/factorized_manifest.tsv" "$OUT_DIR/factorized_ld_manifest.tsv"; then
  echo "error: SLURM manifest contains CRLF line endings" >&2
  exit 1
fi

python3 "$ROOT/test/check_analytic_calibration.py" \
  "${FULL_FACTOR_PREFIX}.mediation" \
  "${FULL_FACTOR_OVERLAP_PREFIX}.mediation"

echo "BMEDIATOR smoke test passed."
