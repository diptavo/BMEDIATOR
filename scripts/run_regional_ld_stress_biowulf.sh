#!/bin/bash
set -euo pipefail

: "${BMEDIATOR_DIR:?Set BMEDIATOR_DIR to the repository directory}"
: "${OUT_DIR:?Set OUT_DIR to the stress-test output directory}"

REPLICATES="${REPLICATES:-500}"
LD_LEVELS="${LD_LEVELS:-0.3 0.6 0.9 0.98}"
WORK_DIR="${LSCRATCH:-$OUT_DIR/work}"

cd "$BMEDIATOR_DIR"
make -j "${SLURM_CPUS_PER_TASK:-2}"

read -r -a ld_args <<< "$LD_LEVELS"
python3 sim/run_regional_ld_stress.py \
  --binary "$BMEDIATOR_DIR/bmediator" \
  --outdir "$OUT_DIR" \
  --workdir "$WORK_DIR" \
  --replicates "$REPLICATES" \
  --ld "${ld_args[@]}"
