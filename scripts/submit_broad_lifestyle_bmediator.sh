#!/usr/bin/env bash
set -euo pipefail

BASE=/data/Dutta_lab/BMEDIATOR
BM_DIR=${BASE}/BMEDIATOR
WORK=${BASE}/analysis/bladder_broad_lifestyle_ukbb_combined_b37
INPUT=${BASE}/generated/bmediator_inputs/bladder_broad_lifestyle
GSMR=/data/Dutta_lab/bladder/gsmr_mr
MANIFEST=${GSMR}/broad_lifestyle_exposures.tsv
PROTEIN_SUMSTAT=${BASE}/generated/UKBB_COMBINED_resource/merged/ukbb_combined_oid_cis_tss_pm6_r2_0.1.tsv
PROTEIN_INFO=${BASE}/generated/UKBB_COMBINED_resource/merged/ukbb_combined_oid_cis_tss_pm6_r2_0.1.protein_info.tsv

mkdir -p "${INPUT}/exposures" "${INPUT}/qc" "${WORK}/manifests" "${WORK}/logs" "${WORK}/results" "${WORK}/combined"

cat > "${WORK}/prepare_broad_lifestyle_inputs.sbatch" <<SBATCH
#!/usr/bin/env bash
#SBATCH --job-name=bmed_life_prep
#SBATCH --partition=norm
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --output=${WORK}/logs/%x_%j.out
#SBATCH --error=${WORK}/logs/%x_%j.err

set -euo pipefail

python3 "${BM_DIR}/scripts/prepare_broad_lifestyle_bmediator_inputs.py" \\
  --manifest "${MANIFEST}" \\
  --gsmr-dir "${GSMR}" \\
  --bim /data/Dutta_lab/REF/EUR.bim \\
  --out-dir "${INPUT}/exposures" \\
  --summary "${INPUT}/qc/preparation_summary.tsv"
SBATCH

python3 "${BM_DIR}/scripts/build_broad_lifestyle_bmediator_manifest.py" \
  --lifestyle-manifest "${MANIFEST}" \
  --exposure-dir "${INPUT}/exposures" \
  --outcome-eur "${BASE}/generated/bmediator_inputs/bladder_bmi_smoking/outcomes/bladder_eur_b37_bmediator.tsv" \
  --outcome-multi "${BASE}/generated/bmediator_inputs/bladder_bmi_smoking/outcomes/bladder_multiancestry_b37_bmediator.tsv" \
  --protein-sumstat "${PROTEIN_SUMSTAT}" \
  --protein-info "${PROTEIN_INFO}" \
  --results-dir "${WORK}/results" \
  --out-manifest "${WORK}/manifests/bladder_broad_lifestyle_manifest.tsv"

N=$(($(wc -l < "${WORK}/manifests/bladder_broad_lifestyle_manifest.tsv") - 1))

cat > "${WORK}/run_broad_lifestyle_array.sbatch" <<'SBATCH'
#!/usr/bin/env bash
#SBATCH --job-name=bmed_blca_life
#SBATCH --partition=norm
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --output=/data/Dutta_lab/BMEDIATOR/analysis/bladder_broad_lifestyle_ukbb_combined_b37/logs/%x_%A_%a.out
#SBATCH --error=/data/Dutta_lab/BMEDIATOR/analysis/bladder_broad_lifestyle_ukbb_combined_b37/logs/%x_%A_%a.err

set -euo pipefail

export BMEDIATOR_DIR=/data/Dutta_lab/BMEDIATOR/BMEDIATOR
export MANIFEST=/data/Dutta_lab/BMEDIATOR/analysis/bladder_broad_lifestyle_ukbb_combined_b37/manifests/bladder_broad_lifestyle_manifest.tsv
export BFILE=/data/Dutta_lab/REF/EUR
export RF_P_THRESH=5e-6
export CIS_P_THRESH=5e-6
export CIS_WINDOW_KB=500
export CLUMP_KB=1000

bash "${BMEDIATOR_DIR}/scripts/run_bmediator_combo_array.sh"
SBATCH

PREP_JOB_ID=$(sbatch --parsable "${WORK}/prepare_broad_lifestyle_inputs.sbatch")
ARRAY_JOB_ID=$(sbatch --parsable --dependency=afterok:${PREP_JOB_ID} --array=1-${N}%4 "${WORK}/run_broad_lifestyle_array.sbatch")
{
  echo "prep_job_id=${PREP_JOB_ID}"
  echo "array_job_id=${ARRAY_JOB_ID}"
  echo "tasks=${N}"
} > "${WORK}/logs/latest_job_ids.txt"
echo "Submitted BMEDIATOR lifestyle prep ${PREP_JOB_ID}; dependent array ${ARRAY_JOB_ID} with ${N} tasks"
