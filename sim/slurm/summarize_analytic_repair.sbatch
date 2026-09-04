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

if [[ "${EVALUATE_BH:-0}" == "1" ]]; then
  python3 sim/evaluate_analytic_bh_validation.py \
    --summary "$OUT_ROOT/summary/analytic_repair_fdr_power.tsv" \
    --output "$OUT_ROOT/summary/analytic_bh_frozen_decision.tsv"
fi
