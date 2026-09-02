# Data Policy

This repository is intended to contain source code, documentation, workflow scripts,
small synthetic test data, and compact curated reference metadata only.

The included `resources/protein_references/` directory is allowed because it
contains small protein mapping/resource tables, compact manifests, and clumped
instrument/reference fixtures created by the BMEDIATOR workflows. These files
are not raw GWAS summary statistics, raw pQTL datasets, PLINK panels, or full
analysis outputs.

Do not commit:

- full GWAS summary statistics
- raw pQTL files
- outcome summary statistics
- PLINK reference panels
- generated BMEDIATOR analysis outputs
- controlled-access, individual-level, or otherwise sensitive data

Use external storage, institutional filesystems, or release artifacts with clear access controls for large data resources.
