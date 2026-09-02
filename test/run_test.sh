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
rm -f "${OUT_PREFIX}.mediation" "${OUT_PREFIX}.hyp" "${OUT_PREFIX}.instruments"

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

if ! awk -F '\t' '
  NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; next}
  $1=="P_SHARED" {
    found=1
    if ($(h["regional_shared_given_both"]) < 0.80 ||
        $(h["P_mediator_identified"]) != $(h["P_mediator_ld_resolved"]) ||
        $(h["mediation_identifiability"]) != "LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL") exit 1
  }
  END {if (!found) exit 1}
' "${FULL_PREFIX}.mediation"; then
  echo "error: shared-signal fixture was not conditionally identified" >&2
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

echo "BMEDIATOR smoke test passed."
