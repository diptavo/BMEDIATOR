#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:10:00

set -euo pipefail
source /etc/profile.d/modules.sh || true
module load R/4.5.2

: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"
: "${REPLICATES_PER_SCENARIO:=50}"

cd "$BMEDIATOR_DIR"
Rscript research/summarize_joint_graph_v021_families.R \
  "$OUT_ROOT" "$OUT_ROOT/summary/joint_graph_family_summary.tsv" \
  "$REPLICATES_PER_SCENARIO"
