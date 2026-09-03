#!/bin/bash
set -euo pipefail

repo="${1:-/data/Dutta_lab/BMEDIATOR/BMEDIATOR}"
config="${2:-$repo/sim/configs/factorized_identification_smoke.json}"
out_root="${3:-/data/Dutta_lab/BMEDIATOR/sim_benchmark/factorized_identification}"
array_throttle="${4:-20}"
block_size="${5:-1}"

manifest_dir="$out_root/manifests"
log_dir="$out_root/logs"
mkdir -p "$manifest_dir" "$log_dir"
manifest="$manifest_dir/benchmark_manifest.tsv"

cd "$repo"
python3 sim/build_benchmark_manifest.py \
  --config "$config" \
  --benchmark classification \
  --out "$manifest"

n_tasks=$(( $(wc -l < "$manifest") - 1 ))
if [[ "$n_tasks" -le 0 ]]; then
  echo "Manifest has no tasks: $manifest" >&2
  exit 1
fi
array_tasks=$(( (n_tasks + block_size - 1) / block_size ))

array_job=$(
  sbatch --parsable \
    -J bmed_factor \
    --partition=quick \
    --array="1-${array_tasks}%${array_throttle}" \
    --cpus-per-task=2 \
    --mem=4g \
    --time=00:30:00 \
    -o "$log_dir/bmed_factor.%A_%a.out" \
    -e "$log_dir/bmed_factor.%A_%a.err" \
    --export=ALL,MANIFEST="$manifest",CONFIG="$config",BMEDIATOR_DIR="$repo",OUT_ROOT="$out_root",BLOCK_SIZE="$block_size",STRUCTURAL_METHOD=factorized \
    "$repo/scripts/run_m1_m5_array.sh"
)

summary_job=$(
  sbatch --parsable \
    -J bmed_factor_summary \
    --partition=quick \
    --dependency="afterok:${array_job}" \
    --cpus-per-task=1 \
    --mem=2g \
    --time=00:10:00 \
    -o "$log_dir/bmed_factor_summary.%j.out" \
    -e "$log_dir/bmed_factor_summary.%j.err" \
    --export=ALL,OUT_ROOT="$out_root",BMEDIATOR_DIR="$repo",STRUCTURAL_METHOD=factorized \
    "$repo/scripts/run_m1_m5_summary.sh"
)

echo "ARRAY_JOB=${array_job}"
echo "SUMMARY_JOB=${summary_job}"
echo "OUTPUT_ROOT=${out_root}"
