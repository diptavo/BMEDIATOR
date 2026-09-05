#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:15:00

set -euo pipefail
source /etc/profile.d/modules.sh || true

: "${BMEDIATOR_DIR:?}"
: "${MANIFEST:?}"
: "${OUT_ROOT:?}"

row=$(awk -F '\t' -v idx="${SLURM_ARRAY_TASK_ID:?}" '
  NR == 1 { next }
  NR - 1 == idx { print; exit }
' "$MANIFEST")
if [[ -z "$row" ]]; then
  echo "No manifest row for array index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi
IFS=$'\t' read -r scenario replicate proteins <<<"$row"

cd "$BMEDIATOR_DIR"
Rscript research/run_joint_graph_v021_family_task.R \
  "$BMEDIATOR_DIR/build/joint_graph_v02_cli" \
  "$scenario" "$replicate" "$OUT_ROOT" "$proteins"
