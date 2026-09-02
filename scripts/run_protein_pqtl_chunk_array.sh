#!/bin/bash
set -euo pipefail

: "${MANIFEST:?}"
: "${BUILDER_BIN:?}"
: "${LD_PREFIX:?}"
: "${OUT_ROOT:?}"
: "${RESOURCE_MODE:?}"

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

IFS=$'\t' read -r chunk_id chunk_label resource_subset n_rows <<<"${row}"
mkdir -p "$OUT_ROOT"

if [[ "$RESOURCE_MODE" == "aric" ]]; then
  "$BUILDER_BIN" "$resource_subset" "$LD_PREFIX" "$OUT_ROOT" "$chunk_label"
elif [[ "$RESOURCE_MODE" == "ukbb" ]]; then
  : "${VARIANT_DIR:?}"
  "$BUILDER_BIN" "$resource_subset" "$VARIANT_DIR" "$LD_PREFIX" "$OUT_ROOT" "$chunk_label"
else
  echo "Unsupported RESOURCE_MODE: $RESOURCE_MODE" >&2
  exit 1
fi
