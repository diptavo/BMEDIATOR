# Frozen LD stress validation

This directory contains the genotype/LD stress-test outputs from the frozen
factorized validation completed on Biowulf on 2026-09-03.

- Slurm array: `29021342` (all tasks completed)
- Summary job: `29021347`
- Analyses: 1,200
- Binary SHA-256: `bd39c66663ca3d7b993ad4964f89344cff765be17d351df0153875c91cbbe193`

`factorized_ld_stress_summary.tsv` is the compact release summary;
`factorized_ld_stress.tsv` contains the run-level results.

Matched-reference distinct loci were usually assigned to H3, while shared
mediation and exactly aligned same-variant pleiotropy remained observationally
similar at H4. Deliberate target/reference LD mismatch inflated Set B and
produced false Bayesian support. See `docs/VALIDATION.md` for the complete
interpretation and release limitations.
