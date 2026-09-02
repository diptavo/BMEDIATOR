#!/bin/bash
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

if [[ -z "${row}" ]]; then
  echo "No manifest row for array index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

IFS=$'\t' read -r benchmark cell replicate <<<"${row}"

echo "[$(date)] benchmark=$benchmark"
echo "[$(date)] cell=$cell"
echo "[$(date)] replicate=$replicate"

cd "$BMEDIATOR_DIR"
python3 sim/run_single_benchmark_task.py \
  --config "$CONFIG" \
  --benchmark "$benchmark" \
  --cell "$cell" \
  --replicate "$replicate" \
  --outdir "$OUT_ROOT" \
  --binary "$BMEDIATOR_DIR/bmediator"
