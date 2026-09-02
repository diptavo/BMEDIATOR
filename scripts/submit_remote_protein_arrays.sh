#!/bin/bash
set -euo pipefail

repo=/data/Dutta_lab/BMEDIATOR/BMEDIATOR
generated=/data/Dutta_lab/BMEDIATOR/generated
ld_prefix=/data/Dutta_lab/REF/EUR
variant_dir=/data/Dutta_lab/decode/variants

chmod +x \
  "$repo/scripts/run_protein_pqtl_chunk_array.sh" \
  "$repo/build/build_aric_pqtl_instruments_multi" \
  "$repo/build/build_ukbb_pqtl_instruments_multi"

aric_manifest="$generated/ARIC_resource/chunks/aric_seqid_tss500kb_chunk_manifest.tsv"
ukbb_eur_manifest="$generated/UKBB_EUR_resource/chunks/ukbb_eur_oid_tss500kb_chunk_manifest.tsv"
ukbb_combined_manifest="$generated/UKBB_COMBINED_resource/chunks/ukbb_combined_oid_tss500kb_chunk_manifest.tsv"

aric_n=$(( $(wc -l < "$aric_manifest") - 1 ))
ukbb_eur_n=$(( $(wc -l < "$ukbb_eur_manifest") - 1 ))
ukbb_combined_n=$(( $(wc -l < "$ukbb_combined_manifest") - 1 ))

aric_job=$(
  sbatch --parsable \
    -J aric_pqtl \
    --array="1-${aric_n}%5" \
    --cpus-per-task=4 \
    --mem=24g \
    --time=05:00:00 \
    -o "$generated/ARIC_resource/logs/aric_pqtl.%A_%a.out" \
    -e "$generated/ARIC_resource/logs/aric_pqtl.%A_%a.err" \
    --export=ALL,MANIFEST="$aric_manifest",BUILDER_BIN="$repo/build/build_aric_pqtl_instruments_multi",LD_PREFIX="$ld_prefix",OUT_ROOT="$generated/ARIC_resource/instruments",RESOURCE_MODE=aric \
    "$repo/scripts/run_protein_pqtl_chunk_array.sh"
)

ukbb_eur_job=$(
  sbatch --parsable \
    -J ukbb_eur_pqtl \
    --array="1-${ukbb_eur_n}%4" \
    --cpus-per-task=4 \
    --mem=32g \
    --time=05:00:00 \
    -o "$generated/UKBB_EUR_resource/logs/ukbb_eur_pqtl.%A_%a.out" \
    -e "$generated/UKBB_EUR_resource/logs/ukbb_eur_pqtl.%A_%a.err" \
    --export=ALL,MANIFEST="$ukbb_eur_manifest",BUILDER_BIN="$repo/build/build_ukbb_pqtl_instruments_multi",LD_PREFIX="$ld_prefix",OUT_ROOT="$generated/UKBB_EUR_resource/instruments",RESOURCE_MODE=ukbb,VARIANT_DIR="$variant_dir" \
    "$repo/scripts/run_protein_pqtl_chunk_array.sh"
)

ukbb_combined_job=$(
  sbatch --parsable \
    -J ukbb_combined_pqtl \
    --array="1-${ukbb_combined_n}%4" \
    --cpus-per-task=4 \
    --mem=32g \
    --time=05:00:00 \
    -o "$generated/UKBB_COMBINED_resource/logs/ukbb_combined_pqtl.%A_%a.out" \
    -e "$generated/UKBB_COMBINED_resource/logs/ukbb_combined_pqtl.%A_%a.err" \
    --export=ALL,MANIFEST="$ukbb_combined_manifest",BUILDER_BIN="$repo/build/build_ukbb_pqtl_instruments_multi",LD_PREFIX="$ld_prefix",OUT_ROOT="$generated/UKBB_COMBINED_resource/instruments",RESOURCE_MODE=ukbb,VARIANT_DIR="$variant_dir" \
    "$repo/scripts/run_protein_pqtl_chunk_array.sh"
)

echo "ARIC_JOB=${aric_job}"
echo "UKBB_EUR_JOB=${ukbb_eur_job}"
echo "UKBB_COMBINED_JOB=${ukbb_combined_job}"
