#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8g
#SBATCH --time=00:20:00

set -euo pipefail
source /etc/profile.d/modules.sh || true

: "${BMEDIATOR_DIR:?}"
cd "$BMEDIATOR_DIR"
make clean
make -j "${SLURM_CPUS_PER_TASK:-1}"
make test
