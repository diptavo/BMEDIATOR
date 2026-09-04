#!/bin/bash
set -euo pipefail

repo="${1:-/data/Dutta_lab/BMEDIATOR/BMEDIATOR}"
out_root="${2:-/data/Dutta_lab/BMEDIATOR/sim_benchmark/factorized_ld_main}"
replicates="${3:-100}"
array_throttle="${4:-20}"
block_size="${5:-1}"

coordinator_job=$(sbatch --parsable \
  --export=ALL,BMEDIATOR_DIR="$repo",OUT_ROOT="$out_root",REPLICATES="$replicates",ARRAY_THROTTLE="$array_throttle",BLOCK_SIZE="$block_size" \
  "$repo/scripts/coordinate_factorized_ld_validation.sh")

echo "COORDINATOR_JOB=${coordinator_job}"
echo "OUTPUT_ROOT=${out_root}"
echo "Child job IDs will be written by the compute-node coordinator to:"
echo "${out_root}/manifests/job_ids.tsv"
