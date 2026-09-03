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
  echo "[$(date)] cell=$cell replicate=$replicate"
  cd "$BMEDIATOR_DIR"
  python3 sim/run_m1_m5_task.py \
    --config "$CONFIG" \
    --cell "$cell" \
    --replicate "$replicate" \
    --outdir "$OUT_ROOT" \
    --binary "$BMEDIATOR_DIR/bmediator"
done

