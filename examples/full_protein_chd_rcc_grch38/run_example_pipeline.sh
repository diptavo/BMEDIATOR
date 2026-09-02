#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_example_pipeline.sh <ukb_ppp_combined|decode> <chd|rcc>

Required environment variables:
  DATA_ROOT       Root of the full_protein_chd_rcc_grch38 example data

Optional environment variables:
  BMEDIATOR_BIN   BMEDIATOR executable (default: repository/bmediator)
  LD_PREFIX       GRCh38 PLINK prefix (default: DATA_ROOT/ld_reference/G1000plink)
  OUT_ROOT        Output directory (default: DATA_ROOT/results)
  THREADS         Number of BMEDIATOR threads (default: 1)
  DRY_RUN         Set to 1 to print commands without running BMEDIATOR
EOF
}

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 2
fi

platform=$1
outcome_name=$2
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/../.." && pwd)

: "${DATA_ROOT:?Set DATA_ROOT to the full_protein_chd_rcc_grch38 directory}"
BMEDIATOR_BIN=${BMEDIATOR_BIN:-$repo_root/bmediator}
OUT_ROOT=${OUT_ROOT:-$DATA_ROOT/results}
THREADS=${THREADS:-1}
DRY_RUN=${DRY_RUN:-0}
LD_PREFIX=${LD_PREFIX:-$DATA_ROOT/ld_reference/G1000plink}

case "$platform" in
    ukb_ppp_combined)
        manifest="$script_dir/manifests/ukb_ppp_combined.protein_gwas_manifest.txt"
        protein_info="protein_panels/ukb_ppp_combined/protein_info.tsv"
        ;;
    decode)
        manifest="$script_dir/manifests/decode.protein_gwas_manifest.txt"
        protein_info="protein_panels/decode/protein_info.tsv"
        ;;
    *)
        echo "Unknown platform: $platform" >&2
        usage >&2
        exit 2
        ;;
esac

case "$outcome_name" in
    chd)
        outcome="outcome/chd_finngen_r12_grch38.bmediator.tsv"
        ;;
    rcc)
        outcome="outcome/rcc_multiancestry_grch38.bmediator.tsv"
        ;;
    *)
        echo "Unknown outcome: $outcome_name" >&2
        usage >&2
        exit 2
        ;;
esac

rf="rf/bmi_giant_locke_eur_grch38.bmediator.tsv"

for required in "$manifest" "$DATA_ROOT/$rf" "$DATA_ROOT/$outcome" "$DATA_ROOT/$protein_info"; do
    if [[ ! -s "$required" ]]; then
        echo "Missing or empty input: $required" >&2
        exit 1
    fi
done

while IFS=$'\t ' read -r protein sumstat_file _; do
    [[ -z "${protein:-}" || "$protein" == \#* ]] && continue
    if [[ ! -s "$DATA_ROOT/$sumstat_file" ]]; then
        echo "Manifest target is missing or empty: $DATA_ROOT/$sumstat_file" >&2
        exit 1
    fi
    if ! awk -v protein="$protein" 'NR > 1 && $1 == protein { found = 1 } END { exit !found }' \
        "$DATA_ROOT/$protein_info"; then
        echo "Protein is absent from annotation file: $protein" >&2
        exit 1
    fi
done < "$manifest"

if [[ "$DRY_RUN" != 1 ]]; then
    if [[ ! -x "$BMEDIATOR_BIN" ]]; then
        echo "BMEDIATOR executable not found: $BMEDIATOR_BIN" >&2
        exit 1
    fi
fi

for extension in bed bim fam; do
    if [[ ! -s "${LD_PREFIX}.${extension}" ]]; then
        echo "Missing LD reference file: ${LD_PREFIX}.${extension}" >&2
        exit 1
    fi
done

run_root="$OUT_ROOT/$platform/$outcome_name"
single_manifest_dir="$run_root/single_protein_manifests"
mkdir -p "$single_manifest_dir" "$run_root/per_protein" "$run_root/manifest"

run_command() {
    if [[ "$DRY_RUN" == 1 ]]; then
        printf 'DRY RUN:'
        printf ' %q' "$@"
        printf '\n'
    else
        "$@"
    fi
}

cd "$DATA_ROOT"

# Run the same full-mode model once for each protein using a one-row manifest.
while IFS=$'\t ' read -r protein sumstat_file _; do
    [[ -z "${protein:-}" || "$protein" == \#* ]] && continue
    single_manifest="$single_manifest_dir/${protein}.protein_gwas_manifest.txt"
    printf '%s\t%s\n' "$protein" "$sumstat_file" > "$single_manifest"
    mkdir -p "$run_root/per_protein/$protein"

    run_command "$BMEDIATOR_BIN" \
        --rf-sumstat "$rf" \
        --protein-gwas-list "$single_manifest" \
        --cancer-sumstat "$outcome" \
        --protein-info "$protein_info" \
        --bfile "$LD_PREFIX" \
        --p-thresh-rf 5e-6 \
        --p-thresh-cis 5e-6 \
        --cis-window 1000 \
        --clump-kb 10000 \
        --clump-r2 0.1 \
        --threads "$THREADS" \
        --out "$run_root/per_protein/$protein/bmi_${outcome_name}_${protein}"
done < "$manifest"

# Analyze all five proteins together. Fixed priors are the default; this run also
# provides the across-protein ranking and local/cumulative FDR columns.
run_command "$BMEDIATOR_BIN" \
    --rf-sumstat "$rf" \
    --protein-gwas-list "$manifest" \
    --cancer-sumstat "$outcome" \
    --protein-info "$protein_info" \
    --bfile "$LD_PREFIX" \
    --p-thresh-rf 5e-6 \
    --p-thresh-cis 5e-6 \
    --cis-window 1000 \
    --clump-kb 10000 \
    --clump-r2 0.1 \
    --threads "$THREADS" \
    --out "$run_root/manifest/bmi_${outcome_name}_${platform}_all5"
