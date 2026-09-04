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
  row=$(awk -F '\t' -v idx="$idx" 'NR > 1 && NR - 1 == idx {print; exit}' "$MANIFEST")
  if [[ -z "$row" ]]; then
    continue
  fi
  IFS=$'\t' read -r benchmark cell replicate <<<"$row"
  echo "[$(date)] benchmark=$benchmark cell=$cell replicate=$replicate"
  task_metrics="$OUT_ROOT/$benchmark/$cell/$(printf 'rep_%04d' "$replicate")/task_metrics.tsv"
  if [[ -s "$task_metrics" ]]; then
    echo "Skipping completed task: $task_metrics"
    continue
  fi
  cd "$BMEDIATOR_DIR"
  python3 sim/run_factorized_task.py \
    --config "$CONFIG" \
    --benchmark "$benchmark" \
    --cell "$cell" \
    --replicate "$replicate" \
    --outdir "$OUT_ROOT" \
    --binary "$BMEDIATOR_DIR/bmediator"
done
