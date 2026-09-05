# JG-0.2.1 held-out family calibration plan

## Purpose and freeze

This is the first family-level release gate for `JG-0.2.1`. It is frozen after
the patch-development matrix and before any held-out family is generated. The
run contains 10 cells, 50 family replicates per cell, and 100 proteins per
family: 50,000 protein analyses in total.

Each protein uses independent seeds beginning at `20400000`. These seeds do not
overlap any earlier joint-graph simulation. Computation must run through
Biowulf batch jobs or an allocated interactive node, never the login node.

## Selection rules

For protein `j`, define posterior local false-discovery probability

```text
lfdr_j = 1 - PP_two_path_j.
```

The primary family rule sorts `lfdr`, selects the largest prefix whose mean is
at most `0.05`, and reports its realized false-discovery proportion and power.
A secondary fixed rule selects `PP_two_path >= 0.80`. Both rules are frozen.
Failed or numerically suppressed analyses are never selected.

The exact-alignment exclusion remains part of every claim. Aligned-pleiotropy
proteins are counted as nonmediators in the sensitivity cell to quantify the
unavoidable failure if that exclusion is false.

## Family cells

| Cell | Composition or perturbation | Status |
| --- | --- | --- |
| `baseline` | 10 moderate mediators, 90 nulls | confirmatory |
| `rare` | 2 moderate mediators, 98 nulls | confirmatory |
| `composite_null` | 20 each null, XM-only, MY-only, sparse, directional | confirmatory |
| `mixed` | 10 mediation, 10 mediation+sparse, 20 sparse, 20 directional, 40 null | confirmatory |
| `strong_ld` | baseline composition, within-block AR(1) `rho=0.70` | confirmatory |
| `ld_mismatch` | true LD `rho=0.70`, analyzed LD `rho=0.50` | diagnostic |
| `scale_uncertainty` | baseline with log-normal external-scale error, log-SD `0.25` | diagnostic |
| `undeclared_overlap` | baseline with nonzero unreported sampling correlations | diagnostic |
| `weak_paths` | 10 weak mediators and 90 nulls | diagnostic |
| `aligned_sensitivity` | 20 exactly aligned pleiotropic proteins and 80 nulls | boundary |

## Frozen acceptance criteria

Using the primary 5% posterior-FDR rule:

- `baseline`: mean FDR at most 0.05 and mean power at least 0.70;
- `rare`: mean FDR at most 0.05 and mean power at least 0.50;
- `composite_null`: mean false discoveries at most 0.05 per family and no more
  than 5% of families with any discovery;
- `mixed`: mean FDR at most 0.05 and mean power at least 0.60;
- `strong_ld`: mean FDR at most 0.05 and mean power at least 0.65.

Every confirmatory cell must have at least 49 of 50 fully completed families.
No diagnostic cell can fail or validate the method. Criteria will not be
changed after results are read.

## Remaining interpretation boundary

Passing this grid would support family-level calibration under the simulated
model and selected perturbations. It would not eliminate exact aligned
pleiotropy, validate the external scale estimates in real cohorts, or replace
head-to-head comparisons with established methods. Those remain separate
release requirements.
