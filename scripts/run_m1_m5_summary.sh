#!/bin/bash
set -euo pipefail

source /etc/profile.d/modules.sh || true

: "${OUT_ROOT:?}"
: "${BMEDIATOR_DIR:?}"

cd "$BMEDIATOR_DIR"
python3 sim/summarize_benchmark.py --outdir "$OUT_ROOT" --rebuild
python3 sim/summarize_m1_m5.py --outdir "$OUT_ROOT"

