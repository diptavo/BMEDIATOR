#!/bin/bash
#SBATCH --job-name=bmed_factor_prior_sensitivity
#SBATCH --partition=quick
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:20:00

set -euo pipefail

: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"

python3 "$BMEDIATOR_DIR/sim/summarize_factorized_bayes_sensitivity.py" \
  --input "$OUT_ROOT" \
  --output "$OUT_ROOT/summary/factorized_bayes_prior_sensitivity.tsv"
