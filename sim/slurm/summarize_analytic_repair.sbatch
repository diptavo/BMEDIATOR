#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=00:20:00

set -euo pipefail
source /etc/profile.d/modules.sh || true

: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"

cd "$BMEDIATOR_DIR"
python3 sim/summarize_analytic_repair.py \
  --input "$OUT_ROOT" \
  --output "$OUT_ROOT/summary/analytic_repair_fdr_power.tsv"
