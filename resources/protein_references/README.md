# Protein Reference Resources

This directory contains the small protein reference resources that were created
with the Biowulf BMEDIATOR workflows. They are included in the GitHub release
because they are compact metadata/instrument resources, not raw GWAS summary
statistics or full pQTL datasets.

## Contents

| Directory | Files | Purpose |
|-----------|-------|---------|
| `aric/` | `aric_seqid_resource.tsv`, `aric_seqid_tss500kb_chunk_manifest.tsv`, `aric_pilot_resource.tsv` | ARIC SeqId-to-gene/TSS resource and small pilot/chunk manifests |
| `ukbb_eur/` | `ukbb_eur_oid_resource.tsv`, `ukbb_eur_oid_tss500kb_chunk_manifest.tsv`, `ukbb_eur_pilot_resource.tsv` | UK Biobank PPP EUR OID protein resource and pilot/chunk manifests |
| `ukbb_combined/` | `ukbb_combined_oid_resource.tsv`, `ukbb_combined_oid_tss500kb_chunk_manifest.tsv`, `ukbb_combined_pilot_resource.tsv` | UK Biobank PPP combined-resource OID protein resource and pilot/chunk manifests |
| `decode/` | `decode_seqid_to_gene_uniprot_tss_hg38.tsv`, missing-mapping reports, build manifest, session info, and build script | deCODE SeqId-to-gene/UniProt/TSS reference and mapping provenance |
| `ea/` | `ea_protein_info_from_seqid.tsv`, `ea_protein_gwas_manifest.txt`, `ea_pqtl_instruments_p5e-6_r2_0.1.tsv` | EA protein-info table, per-protein GWAS manifest, and compact clumped pQTL instrument table |
| `debug/` | tiny one-protein and three-protein fixtures | Small fixtures for checking parser and manifest behavior |

## Usage Notes

- Tables with `sumstat_file` or `sumstat_dir` columns preserve the Biowulf paths
  used when the resources were created. Rewrite those columns to match your
  local or cluster filesystem before running production analyses.
- The `ea_pqtl_instruments_p5e-6_r2_0.1.tsv` file is a compact clumped
  instrument table. It is included as a reference/input resource, but it is not
  a replacement for the full source pQTL summary statistics.
- Raw pQTL, GWAS, PLINK reference panels, and BMEDIATOR output files remain
  intentionally excluded from this GitHub release.
