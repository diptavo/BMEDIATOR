# Release Checklist

Use this checklist before creating a GitHub release.

- Confirm `README.md` has the correct public repository URL.
- Confirm `CITATION.cff` has the correct repository URL and release date.
- Run `make clean && make && make test`.
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
git tag -a v1.0.0 -m "BMEDIATOR v1.0.0"
git push origin main --tags
```

