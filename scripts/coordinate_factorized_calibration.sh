#!/bin/bash
#SBATCH --job-name=bmed_factor_coord
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
: "${BENCHMARK:?}"

manifest_dir="$OUT_ROOT/manifests"
log_dir="$OUT_ROOT/logs"
mkdir -p "$manifest_dir" "$log_dir"
manifest="$manifest_dir/benchmark_manifest.tsv"
submitted_config="$manifest_dir/submitted_config.json"

if find "$OUT_ROOT" -name task_metrics.tsv -print -quit | grep -q .; then
  echo "Output root already contains task results: $OUT_ROOT" >&2
  echo "Choose a new output root to avoid mixing validation runs." >&2
  exit 1
fi

cp "$CONFIG" "$submitted_config"
python3 "$BMEDIATOR_DIR/sim/build_benchmark_manifest.py" \
  --config "$submitted_config" \
  --benchmark "$BENCHMARK" \
  --out "$manifest"

n_tasks=$(( $(wc -l < "$manifest") - 1 ))
if [[ "$n_tasks" -le 0 ]]; then
  echo "Manifest has no tasks: $manifest" >&2
  exit 1
fi
array_tasks=$(( (n_tasks + BLOCK_SIZE - 1) / BLOCK_SIZE ))

array_job=$(sbatch --parsable \
  -J bmed_factor_cal --partition=quick \
  --array="1-${array_tasks}%${ARRAY_THROTTLE}" \
  --cpus-per-task=2 --mem=4g --time=00:30:00 \
  -o "$log_dir/bmed_factor_cal.%A_%a.out" \
  -e "$log_dir/bmed_factor_cal.%A_%a.err" \
  --export=ALL,MANIFEST="$manifest",CONFIG="$submitted_config",BMEDIATOR_DIR="$BMEDIATOR_DIR",OUT_ROOT="$OUT_ROOT",BLOCK_SIZE="$BLOCK_SIZE" \
  "$BMEDIATOR_DIR/scripts/run_factorized_array.sh")

summary_job=$(sbatch --parsable \
  --dependency="afterok:${array_job}" \
  -o "$log_dir/bmed_factor_cal_summary.%j.out" \
  -e "$log_dir/bmed_factor_cal_summary.%j.err" \
  --export=ALL,BMEDIATOR_DIR="$BMEDIATOR_DIR",OUT_ROOT="$OUT_ROOT",CONFIG="$submitted_config" \
  "$BMEDIATOR_DIR/scripts/run_factorized_calibration_summary.sh")

sensitivity_job=$(sbatch --parsable \
  --dependency="afterok:${summary_job}" \
  -o "$log_dir/bmed_factor_prior_sensitivity.%j.out" \
  -e "$log_dir/bmed_factor_prior_sensitivity.%j.err" \
  --export=ALL,BMEDIATOR_DIR="$BMEDIATOR_DIR",OUT_ROOT="$OUT_ROOT" \
  "$BMEDIATOR_DIR/scripts/run_factorized_prior_sensitivity.sh")

printf 'ARRAY_JOB\t%s\nSUMMARY_JOB\t%s\nSENSITIVITY_JOB\t%s\n' \
  "$array_job" "$summary_job" "$sensitivity_job" | tee "$manifest_dir/job_ids.tsv"
