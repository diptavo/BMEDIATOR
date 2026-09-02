#!/bin/bash
set -euo pipefail

ROOT=/data/Dutta_lab/BMEDIATOR
REPO="$ROOT/BMEDIATOR"
INPUT_ROOT="$ROOT/generated/bmediator_inputs/bladder_bmi_smoking"
ANALYSIS_ROOT="$ROOT/analysis/bladder_bmi_smoking_ukbb_combined_b37"
BFILE=/data/Dutta_lab/REF/EUR
BIM="${BFILE}.bim"

BMI_URL='https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/k1t9zK3v2ab8B_UgF8szJ21cRJ_7sEz8q8Pirxacb4cXO7JYATg3Tm689itldnZh/n/ieup4/b/igd/o/ieu-b-40/ieu-b-40.vcf.gz'
SMOKING_4877_URL='https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/1Tjn9cwbSagqducRbJmjR_pZjG6NS2kiSPanj-BJ7Wn_C2nY9Qbx_XMwMa7p83P8/n/ieup4/b/igd/o/ieu-b-4877/ieu-b-4877.vcf.gz'
SMOKING_142_URL='https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/FaEUn6T3INZvlUadVZeGi0kShF-DxtjT9bKoO64kxqbA_tI4LRH_oH1PphFe74mb/n/ieup4/b/igd/o/ieu-b-142/ieu-b-142.vcf.gz'

mkdir -p \
  "$INPUT_ROOT/raw/exposures" \
  "$INPUT_ROOT/exposures" \
  "$INPUT_ROOT/outcomes" \
  "$ANALYSIS_ROOT/logs" \
  "$ANALYSIS_ROOT/jobs" \
  "$ANALYSIS_ROOT/manifests"

prep_job="$ANALYSIS_ROOT/jobs/prepare_bladder_bmi_smoking_inputs.sh"
cat > "$prep_job" <<'EOF'
#!/bin/bash
#SBATCH --job-name=prep_blca_bmed
#SBATCH --cpus-per-task=2
#SBATCH --mem=24g
#SBATCH --time=12:00:00
#SBATCH --output=/data/Dutta_lab/BMEDIATOR/analysis/bladder_bmi_smoking_ukbb_combined_b37/logs/prep.%j.out
#SBATCH --error=/data/Dutta_lab/BMEDIATOR/analysis/bladder_bmi_smoking_ukbb_combined_b37/logs/prep.%j.err
set -euo pipefail

ROOT=/data/Dutta_lab/BMEDIATOR
REPO="$ROOT/BMEDIATOR"
INPUT_ROOT="$ROOT/generated/bmediator_inputs/bladder_bmi_smoking"
ANALYSIS_ROOT="$ROOT/analysis/bladder_bmi_smoking_ukbb_combined_b37"
BIM=/data/Dutta_lab/REF/EUR.bim

download_one() {
  local url="$1"
  local out="$2"
  if [[ ! -s "$out" ]]; then
    wget -nv -c -O "$out" "$url"
  fi
}

download_one "__BMI_URL__" "$INPUT_ROOT/raw/exposures/ieu-b-40.vcf.gz"
download_one "__SMOKING_4877_URL__" "$INPUT_ROOT/raw/exposures/ieu-b-4877.vcf.gz"
download_one "__SMOKING_142_URL__" "$INPUT_ROOT/raw/exposures/ieu-b-142.vcf.gz"

python3 "$REPO/scripts/prepare_bladder_bmediator_sumstats.py" \
  --source ieu_vcf \
  --input "$INPUT_ROOT/raw/exposures/ieu-b-40.vcf.gz" \
  --output "$INPUT_ROOT/exposures/bmi_ieu_b_40_b37_bmediator.tsv" \
  --bim "$BIM" \
  > "$ANALYSIS_ROOT/logs/prep_bmi_ieu_b_40.summary.tsv"

python3 "$REPO/scripts/prepare_bladder_bmediator_sumstats.py" \
  --source ieu_vcf \
  --input "$INPUT_ROOT/raw/exposures/ieu-b-4877.vcf.gz" \
  --output "$INPUT_ROOT/exposures/smoking_ieu_b_4877_b37_bmediator.tsv" \
  --bim "$BIM" \
  > "$ANALYSIS_ROOT/logs/prep_smoking_ieu_b_4877.summary.tsv"

python3 "$REPO/scripts/prepare_bladder_bmediator_sumstats.py" \
  --source ieu_vcf \
  --input "$INPUT_ROOT/raw/exposures/ieu-b-142.vcf.gz" \
  --output "$INPUT_ROOT/exposures/smoking_ieu_b_142_b37_bmediator.tsv" \
  --bim "$BIM" \
  > "$ANALYSIS_ROOT/logs/prep_smoking_ieu_b_142.summary.tsv"

python3 "$REPO/scripts/prepare_bladder_bmediator_sumstats.py" \
  --source bladder_csv \
  --input /data/Dutta_lab/bladder/MVP-NCI-bladderCA-grch37-EUR.txt.gz \
  --output "$INPUT_ROOT/outcomes/bladder_eur_b37_bmediator.tsv" \
  --bim "$BIM" \
  > "$ANALYSIS_ROOT/logs/prep_bladder_eur.summary.tsv"

python3 "$REPO/scripts/prepare_bladder_bmediator_sumstats.py" \
  --source bladder_csv \
  --input /data/Dutta_lab/bladder/MVP-NCI-bladderCA-grch37-Multiancestry.txt.gz \
  --output "$INPUT_ROOT/outcomes/bladder_multiancestry_b37_bmediator.tsv" \
  --bim "$BIM" \
  > "$ANALYSIS_ROOT/logs/prep_bladder_multiancestry.summary.tsv"

python3 "$REPO/scripts/build_bladder_bmediator_manifest.py"
wc -l "$INPUT_ROOT"/exposures/*_bmediator.tsv "$INPUT_ROOT"/outcomes/*_bmediator.tsv
EOF

python3 - "$prep_job" "$BMI_URL" "$SMOKING_4877_URL" "$SMOKING_142_URL" <<'PY'
from pathlib import Path
import sys
path = Path(sys.argv[1])
text = path.read_text()
text = text.replace("__BMI_URL__", sys.argv[2])
text = text.replace("__SMOKING_4877_URL__", sys.argv[3])
text = text.replace("__SMOKING_142_URL__", sys.argv[4])
path.write_text(text)
PY
chmod +x "$prep_job"

run_job="$ANALYSIS_ROOT/jobs/run_bladder_bmediator_array.sh"
cat > "$run_job" <<'EOF'
#!/bin/bash
#SBATCH --job-name=blca_bmed
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=08:00:00
#SBATCH --output=/data/Dutta_lab/BMEDIATOR/analysis/bladder_bmi_smoking_ukbb_combined_b37/logs/run.%A_%a.out
#SBATCH --error=/data/Dutta_lab/BMEDIATOR/analysis/bladder_bmi_smoking_ukbb_combined_b37/logs/run.%A_%a.err
set -euo pipefail

export MANIFEST=/data/Dutta_lab/BMEDIATOR/analysis/bladder_bmi_smoking_ukbb_combined_b37/manifests/bladder_bmi_smoking_manifest.tsv
export BMEDIATOR_DIR=/data/Dutta_lab/BMEDIATOR/BMEDIATOR
export BFILE=/data/Dutta_lab/REF/EUR
export RF_P_THRESH="${RF_P_THRESH:-5e-6}"
export CIS_P_THRESH="${CIS_P_THRESH:-5e-6}"
export CIS_WINDOW_KB="${CIS_WINDOW_KB:-500}"
export CLUMP_KB="${CLUMP_KB:-10000}"

"$BMEDIATOR_DIR/scripts/run_bmediator_combo_array.sh"
EOF
chmod +x "$run_job"

prep_id=$(sbatch --parsable "$prep_job")
run_id=$(sbatch --parsable --dependency=afterok:"$prep_id" --array=1-6%3 "$run_job")

echo "PREP_JOB=${prep_id}"
echo "RUN_JOB=${run_id}"
echo "MANIFEST=$ANALYSIS_ROOT/manifests/bladder_bmi_smoking_manifest.tsv"
echo "RESULTS=$ANALYSIS_ROOT/results"
