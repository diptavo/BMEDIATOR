#!/bin/bash
set -euo pipefail

source /etc/profile.d/modules.sh || true

: "${OUT_ROOT:?}"
: "${BMEDIATOR_DIR:?}"
STRUCTURAL_METHOD="${STRUCTURAL_METHOD:-legacy-six-state}"

cd "$BMEDIATOR_DIR"
if [[ "$STRUCTURAL_METHOD" == "factorized" ]]; then
  python3 sim/summarize_factorized.py \
    --input "$OUT_ROOT" \
    --output "$OUT_ROOT/factorized_summary.tsv"
else
  python3 sim/summarize_benchmark.py --outdir "$OUT_ROOT" --rebuild
  python3 sim/summarize_m1_m5.py --outdir "$OUT_ROOT"
fi
