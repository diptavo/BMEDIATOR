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
  $1=="P_SHARED" {
    found=1
    if ($(h["factor_nA"]) != 2 || $(h["factor_nB"]) != 2 ||
        $(h["factor_log_BF_XM"]) < 2.302585 ||
        $(h["factor_log_BF_MY"]) < 2.302585 ||
        $(h["factor_ld_source"]) != "reference" ||
        $(h["factor_mediation_status"]) != "SUPPORTED_CONDITIONAL") exit 1
  }
  END {if (!found) exit 1}
' "${FULL_FACTOR_PREFIX}.mediation"; then
  echo "error: factorized mode did not recover the deterministic mediation fixture" >&2
  exit 1
fi

if ! awk -F '\t' 'NR==1 {n=NF; next} NF!=n {exit 1}' "${FULL_FACTOR_PREFIX}.regional"; then
  echo "error: regional output rows do not match the declared columns" >&2
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
  END {if (NR != 3) exit 1}
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
        $(h["P_mediator_identified"]) != $(h["P_mediator_ld_resolved"]) ||
        $(h["mediation_identifiability"]) != "LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL") exit 1
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
  --cell identified_setb2 \
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
        m1_bf < 3 || m5_bf > 1) exit 1
  }
' "$M1_M5_FACTOR_PREFIX/classification/identified_setb2/rep_0001/task_metrics.tsv"; then
  echo "error: factorized engine did not separate the identified M1/M5 smoke fixture" >&2
  exit 1
fi

echo "BMEDIATOR smoke test passed."
