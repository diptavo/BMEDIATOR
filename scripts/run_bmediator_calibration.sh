#!/bin/bash
set -euo pipefail

: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"

cd "$BMEDIATOR_DIR"
python3 sim/calibrate_bmediator.py \
  --outdir "$OUT_ROOT" \
  --train-benchmark "${TRAIN_BENCHMARK:-calibration}" \
  --holdout-mod "${HOLDOUT_MOD:-5}" \
  --ridge "${RIDGE:-1.0}" \
  --model-name "${MODEL_NAME:-calibration_model}"
