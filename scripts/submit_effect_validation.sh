#!/bin/bash
set -euo pipefail

repo="${1:-/data/Dutta_lab/BMEDIATOR/BMEDIATOR}"
config="${2:-$repo/sim/configs/factorized_effect_validation.json}"
out_root="${3:-/data/Dutta_lab/BMEDIATOR/sim_benchmark/factorized_effect_validation_20260906_v1}"
array_throttle="${4:-40}"
block_size="${5:-1}"

coordinator_job=$(sbatch --parsable \
  --export=ALL,BMEDIATOR_DIR="$repo",CONFIG="$config",OUT_ROOT="$out_root",ARRAY_THROTTLE="$array_throttle",BLOCK_SIZE="$block_size" \
  "$repo/scripts/coordinate_effect_validation.sh")

echo "COORDINATOR_JOB=${coordinator_job}"
echo "OUTPUT_ROOT=${out_root}"
echo "Child job IDs: ${out_root}/manifests/job_ids.tsv"
