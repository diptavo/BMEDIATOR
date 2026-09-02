#!/bin/bash
set -euo pipefail

source /etc/profile.d/modules.sh || true

: "${MANIFEST:?}"
: "${CONFIG:?}"
: "${BMEDIATOR_DIR:?}"
: "${OUT_ROOT:?}"
: "${BLOCK_SIZE:?}"

start=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK_SIZE + 1 ))
end=$(( SLURM_ARRAY_TASK_ID * BLOCK_SIZE ))

for idx in $(seq "$start" "$end"); do
  row=$(
    awk -F '\t' -v idx="$idx" '
      NR == 1 { next }
      NR - 1 == idx { print; exit }
    ' "$MANIFEST"
  )

  if [[ -z "${row}" ]]; then
    continue
  fi

  IFS=$'\t' read -r benchmark cell replicate <<<"${row}"

  echo "[$(date)] block_task=${SLURM_ARRAY_TASK_ID} manifest_row=$idx"
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
done
