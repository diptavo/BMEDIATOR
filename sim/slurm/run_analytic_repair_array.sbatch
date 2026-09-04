#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=12g
#SBATCH --time=01:00:00

set -euo pipefail
source /etc/profile.d/modules.sh || true

: "${MANIFEST:?}"
: "${CONFIG:?}"
: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"

row=$(
  awk -F '\t' -v idx="${SLURM_ARRAY_TASK_ID:?}" '
    NR == 1 { next }
    NR - 1 == idx { print; exit }
  ' "$MANIFEST"
)
if [[ -z "$row" ]]; then
  echo "No manifest row for array index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi
IFS=$'\t' read -r benchmark cell replicate <<<"$row"

cd "$BMEDIATOR_DIR"
python3 sim/run_factorized_task.py \
  --config "$CONFIG" \
  --benchmark "$benchmark" \
  --cell "$cell" \
  --replicate "$replicate" \
  --outdir "$OUT_ROOT" \
  --binary "$BMEDIATOR_DIR/bmediator"
