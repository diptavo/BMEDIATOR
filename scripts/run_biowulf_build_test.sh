#!/bin/bash
#SBATCH --job-name=bmed_build_test
#SBATCH --partition=quick
#SBATCH --cpus-per-task=2
#SBATCH --mem=4g
#SBATCH --time=00:15:00

set -euo pipefail

: "${BMEDIATOR_DIR:?Set BMEDIATOR_DIR to the repository path}"

cd "$BMEDIATOR_DIR"
make clean
make -j "${SLURM_CPUS_PER_TASK:-2}"
make test
