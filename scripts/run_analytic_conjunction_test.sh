#!/bin/bash
set -euo pipefail

: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"

cd "$BMEDIATOR_DIR"
python3 sim/test_analytic_conjunction.py \
  --outdir "$OUT_ROOT" \
  --output-name "${OUTPUT_NAME:-analytic_conjunction_test}"
