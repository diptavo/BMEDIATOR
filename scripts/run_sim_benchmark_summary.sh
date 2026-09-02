#!/bin/bash
set -euo pipefail

source /etc/profile.d/modules.sh || true

: "${OUT_ROOT:?}"
: "${BMEDIATOR_DIR:?}"

cd "$BMEDIATOR_DIR"
python3 sim/summarize_benchmark.py --outdir "$OUT_ROOT" --rebuild
python3 sim/run_competitor_benchmark.py --outdir "$OUT_ROOT" --benchmark all
