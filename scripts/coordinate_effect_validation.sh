#!/bin/bash
#SBATCH --job-name=bmed_effect_coord
#SBATCH --partition=quick
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --time=00:10:00

set -euo pipefail

: "${BMEDIATOR_DIR:?}"
: "${CONFIG:?}"
: "${OUT_ROOT:?}"
: "${ARRAY_THROTTLE:?}"
: "${BLOCK_SIZE:?}"

manifest_dir="$OUT_ROOT/manifests"
log_dir="$OUT_ROOT/logs"
mkdir -p "$manifest_dir" "$log_dir"
manifest="$manifest_dir/benchmark_manifest.tsv"
submitted_config="$manifest_dir/submitted_config.json"

if find "$OUT_ROOT" -name task_metrics.tsv -print -quit | grep -q .; then
  echo "Output root already contains task results: $OUT_ROOT" >&2
  exit 1
fi

cp "$CONFIG" "$submitted_config"
python3 "$BMEDIATOR_DIR/sim/build_benchmark_manifest.py" \
  --config "$submitted_config" \
  --benchmark calibration \
  --out "$manifest"

if [[ "$BLOCK_SIZE" -le 0 ]]; then
  echo "BLOCK_SIZE must be positive" >&2
  exit 1
fi
n_tasks=$(( $(wc -l < "$manifest") - 1 ))
if [[ "$n_tasks" -le 0 ]]; then
  echo "Manifest has no tasks: $manifest" >&2
  exit 1
fi
array_tasks=$(( (n_tasks + BLOCK_SIZE - 1) / BLOCK_SIZE ))
if [[ "$array_tasks" -gt 1000 ]]; then
  echo "Array requires $array_tasks elements; increase BLOCK_SIZE to stay at or below 1000" >&2
  exit 1
fi

build_job=$(sbatch --parsable \
  -J bmed_effect_build \
  -o "$log_dir/bmed_effect_build.%j.out" \
  -e "$log_dir/bmed_effect_build.%j.err" \
  --export=ALL,BMEDIATOR_DIR="$BMEDIATOR_DIR" \
  "$BMEDIATOR_DIR/sim/slurm/build_and_test.sbatch")

array_job=$(sbatch --parsable \
  --dependency="afterok:${build_job}" \
  -J bmed_effect_sim --partition=quick \
  --array="1-${array_tasks}%${ARRAY_THROTTLE}" \
  --cpus-per-task=2 --mem=4g --time=00:30:00 \
  -o "$log_dir/bmed_effect_sim.%A_%a.out" \
  -e "$log_dir/bmed_effect_sim.%A_%a.err" \
  --export=ALL,MANIFEST="$manifest",CONFIG="$submitted_config",BMEDIATOR_DIR="$BMEDIATOR_DIR",OUT_ROOT="$OUT_ROOT",BLOCK_SIZE="$BLOCK_SIZE" \
  "$BMEDIATOR_DIR/scripts/run_factorized_array.sh")

summary_job=$(sbatch --parsable \
  --dependency="afterok:${array_job}" \
  -o "$log_dir/bmed_effect_summary.%j.out" \
  -e "$log_dir/bmed_effect_summary.%j.err" \
  --export=ALL,BMEDIATOR_DIR="$BMEDIATOR_DIR",OUT_ROOT="$OUT_ROOT",CONFIG="$submitted_config" \
  "$BMEDIATOR_DIR/scripts/run_effect_validation_summary.sh")

printf 'BUILD_JOB\t%s\nARRAY_JOB\t%s\nSUMMARY_JOB\t%s\n' \
  "$build_job" "$array_job" "$summary_job" | tee "$manifest_dir/job_ids.tsv"
