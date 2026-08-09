# Contributing

Contributions are welcome through GitHub issues and pull requests.

## Development Setup

```bash
make
make test
```

OpenMP is optional. Use it only with compilers that support `-fopenmp`:

```bash
make clean
make USE_OPENMP=1
make test
```

## Pull Request Checklist

- Run `make test`.
- Keep large GWAS, pQTL, outcome, PLINK, and generated analysis files out of git.
- Add or update smoke-test inputs only when they are small and synthetic.
- Document user-facing behavior changes in `README.md`.
- Avoid committing machine-specific paths unless they are clearly marked as Biowulf examples.

## Code Style

The core executable is C++17. Prefer small, explicit changes and avoid broad refactors unless they are necessary for correctness or maintainability.

