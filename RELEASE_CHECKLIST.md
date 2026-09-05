# Release Checklist

Use this checklist before creating a GitHub release.

- Confirm `README.md` has the correct public repository URL.
- Confirm `CITATION.cff` has the correct repository URL and release date.
- Confirm the release status in `docs/PRODUCTION_READINESS.md` is current.
- Confirm `bmediator-joint` is used for joint mediation inference; legacy and
  factorized outputs must be labeled as compatibility/development results.
- Confirm the target ancestry and genome build are stated for each real-data analysis.
- Confirm RF, protein, outcome, and signed-LD alleles are harmonized to the
  same effect allele and genome build.
- Confirm A/B/C roles, LD blocks, role scales, and orientation probabilities
  were obtained independently of the analyzed outcome.
- Confirm sampling-error correlations are declared and zero is justified when
  used.
- Confirm the protein family and priors were fixed before outcome inspection.
- Confirm the manifest has no failures and every row says
  `family_complete=TRUE` and `posterior_fdr_status=AVAILABLE_COMPLETE_MANIFEST`.
- Confirm the exact-aligned-pleiotropy identification boundary is disclosed.
- Review `docs/IDENTIFICATION.md` and multi-signal simulation results before removing the `-dev` suffix.
- Run `make clean && make && make test`.
- Run `make test-regional-stress`.
- Optionally run `make clean && make USE_OPENMP=1 && make test` on a Linux system with OpenMP support.
- Confirm no large data files are present:

```bash
find . -type f -size +20M -print
```

- Confirm no generated result files are present:

```bash
find . \( -name "*.mediation" -o -name "*.regional" -o -name "*.hyp" -o -name "*.instruments" \) -print
```

- Tag the release:

```bash
git tag -a v1.2.0 -m "BMEDIATOR v1.2.0"
git push origin main --tags
```
