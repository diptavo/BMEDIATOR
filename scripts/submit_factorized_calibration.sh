#!/bin/bash
set -euo pipefail

repo="${1:-/data/Dutta_lab/BMEDIATOR/BMEDIATOR}"
config="${2:-$repo/sim/configs/factorized_calibration_smoke.json}"
out_root="${3:-/data/Dutta_lab/BMEDIATOR/sim_benchmark/factorized_calibration_smoke}"
array_throttle="${4:-20}"
block_size="${5:-1}"
benchmark="${6:-all}"

coordinator_job=$(sbatch --parsable \
  --export=ALL,BMEDIATOR_DIR="$repo",CONFIG="$config",OUT_ROOT="$out_root",ARRAY_THROTTLE="$array_throttle",BLOCK_SIZE="$block_size",BENCHMARK="$benchmark" \
  "$repo/scripts/coordinate_factorized_calibration.sh")

echo "COORDINATOR_JOB=${coordinator_job}"
echo "OUTPUT_ROOT=${out_root}"
echo "Child job IDs will be written by the compute-node coordinator to:"
echo "${out_root}/manifests/job_ids.tsv"
