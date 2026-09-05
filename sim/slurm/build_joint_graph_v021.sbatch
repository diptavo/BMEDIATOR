#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8g
#SBATCH --time=00:10:00

set -euo pipefail
source /etc/profile.d/modules.sh || true

: "${BMEDIATOR_DIR:?}"
cd "$BMEDIATOR_DIR"
make joint-graph-v02-prototype
