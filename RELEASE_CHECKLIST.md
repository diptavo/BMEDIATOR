# Release Checklist

Use this checklist before creating a GitHub release.

- Confirm `README.md` has the correct public repository URL.
- Confirm `CITATION.cff` has the correct repository URL and release date.
- Confirm the target ancestry and genome build are stated for each real-data analysis.
- Confirm an ancestry-matched PLINK LD reference and unpruned regional pQTL/outcome statistics are available.
- Confirm protein/outcome sample overlap is absent, negligible, or explicitly disclosed as a limitation.
- Review `docs/IDENTIFICATION.md` and current genotype-based LD stress results before removing the `-dev` suffix.
- Run `make clean && make && make test`.
- Run `make test-regional-stress`.
- Optionally run `make clean && make USE_OPENMP=1 && make test` on a Linux system with OpenMP support.
- Confirm no large data files are present:

```bash
find . -type f -size +20M -print
```

- Confirm no generated result files are present:

```bash
find . \( -name "*.mediation" -o -name "*.hyp" -o -name "*.instruments" \) -print
```

- Tag the release:

```bash
git tag -a v1.2.0 -m "BMEDIATOR v1.2.0"
git push origin main --tags
```
