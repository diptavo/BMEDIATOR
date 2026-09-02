#!/bin/bash
set -euo pipefail

BASE=/data/Dutta_lab/BMEDIATOR/exposure_sumstats/2026_04_16

mkdir -p \
  "$BASE/raw/anthropometric/mvp" \
  "$BASE/raw/anthropometric/giant" \
  "$BASE/manifests" \
  "$BASE/logs" \
  "$BASE/jobs"

cat > "$BASE/manifests/exposure_sources_initial.tsv" <<'EOF'
trait	source	label	url	outdir	outfile
BMI	MVP	mvp_bmi_gcst90475156	http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475156	/data/Dutta_lab/BMEDIATOR/exposure_sumstats/2026_04_16/raw/anthropometric/mvp	GCST90475156
BMI	GIANT	giant_bmi_speliotes2010	https://giant-consortium.web.broadinstitute.org/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz	/data/Dutta_lab/BMEDIATOR/exposure_sumstats/2026_04_16/raw/anthropometric/giant	SNP_gwas_mc_merge_nogc.tbl.uniq.gz
EOF

cat > "$BASE/jobs/download_initial_exposures.sh" <<'EOF'
#!/bin/bash
#SBATCH --job-name=exp_dl_init
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=04:00:00
#SBATCH --output=/data/Dutta_lab/BMEDIATOR/exposure_sumstats/2026_04_16/logs/exp_dl_init.%j.out
#SBATCH --error=/data/Dutta_lab/BMEDIATOR/exposure_sumstats/2026_04_16/logs/exp_dl_init.%j.err
set -euo pipefail

BASE=/data/Dutta_lab/BMEDIATOR/exposure_sumstats/2026_04_16
MANIFEST="$BASE/manifests/exposure_sources_initial.tsv"

while IFS=$'\t' read -r trait source label url outdir outfile; do
  if [[ "$trait" == "trait" ]]; then
    continue
  fi
  mkdir -p "$outdir"
  echo "[$(date)] downloading $label"
  wget -nv -c -O "$outdir/$outfile" "$url"
  ls -lh "$outdir/$outfile"
done < "$MANIFEST"
EOF

chmod +x "$BASE/jobs/download_initial_exposures.sh"
sbatch "$BASE/jobs/download_initial_exposures.sh"
