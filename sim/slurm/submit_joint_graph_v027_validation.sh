#!/bin/bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 RUN_ROOT SOURCE_ARCHIVE OLD_RESULT_ROOT" >&2
  exit 2
fi

run_root=$1
source_archive=$2
old_result_root=$3
repo="$run_root/repo"
results="$run_root/results"
manifest="$run_root/family_manifest.tsv"
seed_base=80400000

mkdir -p "$run_root/logs" "$results"

build_job=$(sbatch --parsable \
  --job-name=jg027_build \
  --cpus-per-task=4 --mem=8g --time=00:20:00 \
  --output="$run_root/logs/build.%j.out" \
  --wrap="set -euo pipefail; source /etc/profile.d/modules.sh || true; module load R/4.5.2; mkdir -p '$repo'; tar -xzf '$source_archive' -C '$repo'; cd '$repo'; make test")

replay_job=$(sbatch --parsable \
  --dependency="afterok:$build_job" \
  --job-name=jg027_replay66 \
  --cpus-per-task=1 --mem=8g --time=02:00:00 \
  --output="$run_root/logs/replay66.%j.out" \
  --export="ALL,BMEDIATOR_DIR=$repo,OLD_RESULT_ROOT=$old_result_root,OUTPUT=$results/replay66.tsv,SEED_BASE=50400000" \
  "$run_root/replay_joint_graph_failures.sbatch")

manifest_job=$(sbatch --parsable \
  --dependency="afterok:$build_job" \
  --job-name=jg027_manifest \
  --cpus-per-task=1 --mem=2g --time=00:05:00 \
  --output="$run_root/logs/manifest.%j.out" \
  --wrap="set -euo pipefail; source /etc/profile.d/modules.sh || true; module load R/4.5.2; cd '$repo'; Rscript research/build_joint_graph_v021_family_manifest.R '$manifest' 10 100")

family_job=$(sbatch --parsable \
  --dependency="afterok:$manifest_job" \
  --array=1-100 --job-name=jg027_family \
  --cpus-per-task=1 --mem=6g --time=00:30:00 \
  --output="$run_root/logs/family.%A_%a.out" \
  --export="ALL,BMEDIATOR_DIR=$repo,MANIFEST=$manifest,OUT_ROOT=$results/families,SEED_BASE=$seed_base" \
  "$run_root/run_joint_graph_v021_family_array.sbatch")

summary_job=$(sbatch --parsable \
  --dependency="afterok:$family_job" \
  --job-name=jg027_summary \
  --cpus-per-task=1 --mem=4g --time=00:10:00 \
  --output="$run_root/logs/summary.%j.out" \
  --wrap="set -euo pipefail; source /etc/profile.d/modules.sh || true; module load R/4.5.2; cd '$repo'; mkdir -p '$results/summary'; Rscript research/summarize_joint_graph_v021_families.R '$results/families' '$results/summary/family_summary.tsv' 10; Rscript research/summarize_joint_graph_numerics.R '$results/families' '$results/summary/numerical_summary.tsv'")

printf 'build_job=%s\nreplay_job=%s\nmanifest_job=%s\nfamily_job=%s\nsummary_job=%s\n' \
  "$build_job" "$replay_job" "$manifest_job" "$family_job" "$summary_job"
