#!/bin/bash
#SBATCH --job-name=bmed_effect_summary
#SBATCH --partition=quick
#SBATCH --cpus-per-task=2
#SBATCH --mem=12g
#SBATCH --time=01:00:00

set -euo pipefail

: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"
: "${CONFIG:?}"

python3 "$BMEDIATOR_DIR/sim/summarize_factorized_calibration.py" \
  --input "$OUT_ROOT" \
  --output-dir "$OUT_ROOT/summary" \
  --config "$CONFIG" \
  --require-complete

python3 "$BMEDIATOR_DIR/sim/summarize_effect_validation.py" \
  --input "$OUT_ROOT" \
  --output-dir "$OUT_ROOT/summary" \
  --config "$CONFIG" \
  --require-complete

python3 "$BMEDIATOR_DIR/sim/evaluate_effect_validation.py" \
  --fcr-summary "$OUT_ROOT/summary/effect_fcr_summary.tsv" \
  --effect-summary "$OUT_ROOT/summary/effect_estimation_summary.tsv" \
  --output "$OUT_ROOT/summary/effect_validation_frozen_decision.tsv" \
  --require-pass
