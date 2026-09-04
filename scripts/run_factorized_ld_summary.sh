#!/bin/bash
#SBATCH --job-name=bmed_factor_ld_summary
#SBATCH --partition=quick
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --time=00:10:00

set -euo pipefail

: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"
: "${REPLICATES:?}"

python3 "$BMEDIATOR_DIR/sim/summarize_factorized_ld.py" \
  --input "$OUT_ROOT" \
  --replicates "$REPLICATES" \
  --require-complete
