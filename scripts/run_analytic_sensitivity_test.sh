#!/bin/bash
set -euo pipefail

BMEDIATOR_DIR=${BMEDIATOR_DIR:-/data/Dutta_lab/BMEDIATOR/BMEDIATOR}
OUT_ROOT=${OUT_ROOT:-/data/Dutta_lab/BMEDIATOR/sim_benchmark/ng_stress_main}
OUTPUT_NAME=${OUTPUT_NAME:-analytic_sensitivity_test}
BIAS_BOUNDS=${BIAS_BOUNDS:-0,0.05,0.10,0.15,0.20,0.25}
THRESHOLDS=${THRESHOLDS:-0.01,0.05,0.10}

cd "$BMEDIATOR_DIR"
python3 sim/test_analytic_sensitivity.py \
  --outdir "$OUT_ROOT" \
  --output-name "$OUTPUT_NAME" \
  --bias-bounds "$BIAS_BOUNDS" \
  --thresholds "$THRESHOLDS"
