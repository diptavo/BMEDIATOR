#!/bin/bash
set -euo pipefail

repo="${1:-/data/Dutta_lab/BMEDIATOR/BMEDIATOR}"
config="${2:-$repo/sim/configs/classification_calibration_default.json}"
out_root="${3:-/data/Dutta_lab/BMEDIATOR/sim_benchmark/default_run}"
benchmark="${4:-all}"

manifest_dir="$out_root/manifests"
log_dir="$out_root/logs"
mkdir -p "$manifest_dir" "$log_dir"

manifest="$manifest_dir/benchmark_manifest.tsv"

cd "$repo"
python3 sim/build_benchmark_manifest.py \
  --config "$config" \
  --benchmark "$benchmark" \
  --out "$manifest"

n_tasks=$(( $(wc -l < "$manifest") - 1 ))
if [[ "$n_tasks" -le 0 ]]; then
  echo "Manifest has no tasks: $manifest" >&2
  exit 1
fi

chmod +x \
  "$repo/scripts/run_sim_benchmark_array.sh" \
  "$repo/scripts/run_sim_benchmark_summary.sh"

array_job=$(
  sbatch --parsable \
    -J bmed_sim \
    --array="1-${n_tasks}%20" \
    --cpus-per-task=4 \
    --mem=32g \
    --time=02:45:00 \
    -o "$log_dir/bmed_sim.%A_%a.out" \
    -e "$log_dir/bmed_sim.%A_%a.err" \
    --export=ALL,MANIFEST="$manifest",CONFIG="$config",BMEDIATOR_DIR="$repo",OUT_ROOT="$out_root" \
    "$repo/scripts/run_sim_benchmark_array.sh"
)

summary_job=$(
  sbatch --parsable \
    -J bmed_sim_summary \
    --dependency="afterok:${array_job}" \
    --cpus-per-task=2 \
    --mem=8g \
    --time=00:45:00 \
    -o "$log_dir/bmed_sim_summary.%j.out" \
    -e "$log_dir/bmed_sim_summary.%j.err" \
    --export=ALL,OUT_ROOT="$out_root",BMEDIATOR_DIR="$repo" \
    "$repo/scripts/run_sim_benchmark_summary.sh"
)

echo "ARRAY_JOB=${array_job}"
echo "SUMMARY_JOB=${summary_job}"
