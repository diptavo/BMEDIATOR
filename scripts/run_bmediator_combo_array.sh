#!/bin/bash
set -euo pipefail

: "${MANIFEST:?}"
: "${BMEDIATOR_DIR:?}"
: "${BFILE:?}"

RF_P_THRESH="${RF_P_THRESH:-5e-6}"
CIS_P_THRESH="${CIS_P_THRESH:-5e-6}"
CIS_WINDOW_KB="${CIS_WINDOW_KB:-500}"
CLUMP_KB="${CLUMP_KB:-10000}"

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

IFS=$'\t' read -r domain exposure outcome rf_sumstat outcome_sumstat protein_sumstat protein_info out_prefix <<<"${row}"
out_dir=$(dirname "$out_prefix")
mkdir -p "$out_dir"

echo "[$(date)] domain=$domain"
echo "[$(date)] exposure=$exposure"
echo "[$(date)] outcome=$outcome"
echo "[$(date)] rf_sumstat=$rf_sumstat"
echo "[$(date)] outcome_sumstat=$outcome_sumstat"
echo "[$(date)] protein_sumstat=$protein_sumstat"
echo "[$(date)] protein_info=$protein_info"
echo "[$(date)] out_prefix=$out_prefix"

cd "$BMEDIATOR_DIR"
./bmediator \
  --rf-sumstat "$rf_sumstat" \
  --pqtl-sumstat "$protein_sumstat" \
  --cancer-sumstat "$outcome_sumstat" \
  --protein-info "$protein_info" \
  --bfile "$BFILE" \
  --p-thresh-rf "$RF_P_THRESH" \
  --p-thresh-cis "$CIS_P_THRESH" \
  --cis-window "$CIS_WINDOW_KB" \
  --clump-kb "$CLUMP_KB" \
  --threads "${SLURM_CPUS_PER_TASK:-1}" \
  --verbose \
  --out "$out_prefix"
