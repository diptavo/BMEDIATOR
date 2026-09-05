#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00

set -euo pipefail
source /etc/profile.d/modules.sh || true
module load R/4.5.2

: "${BMEDIATOR_DIR:?}"
: "${OLD_RESULT_ROOT:?}"
: "${OUTPUT:?}"
: "${SEED_BASE:?}"

cd "$BMEDIATOR_DIR"
Rscript research/replay_joint_graph_failures.R \
  "$BMEDIATOR_DIR/build/joint_graph_v02_cli" \
  "$OLD_RESULT_ROOT" "$OUTPUT" "$SEED_BASE"
