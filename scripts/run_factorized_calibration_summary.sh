#!/bin/bash
#SBATCH --job-name=bmed_factor_cal_summary
#SBATCH --partition=quick
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:15:00

set -euo pipefail

: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"
: "${CONFIG:?}"

python3 "$BMEDIATOR_DIR/sim/summarize_factorized_calibration.py" \
  --input "$OUT_ROOT" \
  --output-dir "$OUT_ROOT/summary" \
  --config "$CONFIG" \
  --require-complete
