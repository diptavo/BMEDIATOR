#!/bin/bash
set -euo pipefail

repo=/data/Dutta_lab/BMEDIATOR/BMEDIATOR
generated=/data/Dutta_lab/BMEDIATOR/generated

module load R/4.5

mkdir -p \
  "$generated/ARIC_resource/resources" \
  "$generated/ARIC_resource/chunks" \
  "$generated/ARIC_resource/logs" \
  "$generated/UKBB_EUR_resource/resources" \
  "$generated/UKBB_EUR_resource/chunks" \
  "$generated/UKBB_EUR_resource/logs" \
  "$generated/UKBB_COMBINED_resource/resources" \
  "$generated/UKBB_COMBINED_resource/chunks" \
  "$generated/UKBB_COMBINED_resource/logs"

python3 "$repo/scripts/build_aric_resource.py" \
  --seqid-table /data/Dutta_lab/aric/seqid.txt \
  --sumstat-dir /data/Dutta_lab/aric/EA \
  --out "$generated/ARIC_resource/resources/aric_seqid_resource.tsv"

Rscript "$repo/scripts/build_ukbb_oid_resource.R" \
  /data/Dutta_lab/decode/ukb_ppp/EUR \
  "$generated/UKBB_EUR_resource/resources/ukbb_eur_oid_resource.tsv" \
  UKBB_EUR

Rscript "$repo/scripts/build_ukbb_oid_resource.R" \
  /data/Dutta_lab/decode/ukb_ppp/COMBINED \
  "$generated/UKBB_COMBINED_resource/resources/ukbb_combined_oid_resource.tsv" \
  UKBB_COMBINED

python3 "$repo/scripts/split_resource_chunks.py" \
  --resource "$generated/ARIC_resource/resources/aric_seqid_resource.tsv" \
  --chunk-size 1000 \
  --out-dir "$generated/ARIC_resource/chunks" \
  --prefix aric_seqid_tss500kb

python3 "$repo/scripts/split_resource_chunks.py" \
  --resource "$generated/UKBB_EUR_resource/resources/ukbb_eur_oid_resource.tsv" \
  --chunk-size 500 \
  --out-dir "$generated/UKBB_EUR_resource/chunks" \
  --prefix ukbb_eur_oid_tss500kb

python3 "$repo/scripts/split_resource_chunks.py" \
  --resource "$generated/UKBB_COMBINED_resource/resources/ukbb_combined_oid_resource.tsv" \
  --chunk-size 500 \
  --out-dir "$generated/UKBB_COMBINED_resource/chunks" \
  --prefix ukbb_combined_oid_tss500kb

wc -l \
  "$generated/ARIC_resource/resources/aric_seqid_resource.tsv" \
  "$generated/UKBB_EUR_resource/resources/ukbb_eur_oid_resource.tsv" \
  "$generated/UKBB_COMBINED_resource/resources/ukbb_combined_oid_resource.tsv" \
  "$generated/ARIC_resource/chunks/aric_seqid_tss500kb_chunk_manifest.tsv" \
  "$generated/UKBB_EUR_resource/chunks/ukbb_eur_oid_tss500kb_chunk_manifest.tsv" \
  "$generated/UKBB_COMBINED_resource/chunks/ukbb_combined_oid_tss500kb_chunk_manifest.tsv"
