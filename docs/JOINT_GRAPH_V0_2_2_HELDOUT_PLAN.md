# JG-0.2.2 final numerical-repair calibration plan

## Freeze and purpose

This plan is frozen before any `JG-0.2.2` family result is generated. The
`JG-0.2.1` family run exposed four marginal failures of its
adaptive-versus-Laplace diagnostic. That run is now development evidence and
will not be reused as final validation.

`JG-0.2.2` replaces that diagnostic with convergence across mode-centered
adaptive Gauss-Hermite rules. Every graph state with preliminary posterior
probability above `1e-6` is evaluated at orders 3 and 5. A difference above
`0.05` triggers order 7 and, if still above `0.05`, order 9. A posterior is
suppressed when the last two evaluated orders differ by more than `0.10` log
units. The structural likelihood, graph states, priors, and selection rule are
otherwise unchanged.

## Independent run

The run uses seed base `30400000`, which is disjoint from all JG-0.1,
JG-0.2, JG-0.2.1 development, and first family-run seeds. It contains the same
10 scenarios, 50 family replicates per scenario, and 100 proteins per family:
50,000 protein analyses. It must run on Biowulf batch or allocated interactive
nodes, never the login node.

## Frozen selection and criteria

The primary rule is unchanged: sort `lfdr = 1 - PP_two_path` and select the
largest prefix with mean `lfdr <= 0.05`. The secondary diagnostic rule remains
`PP_two_path >= 0.80`. Failed analyses cannot be selected.

The confirmatory criteria are unchanged from the first family plan:

- `baseline`: mean FDR at most 0.05 and mean power at least 0.70;
- `rare`: mean FDR at most 0.05 and mean power at least 0.50;
- `composite_null`: at most 0.05 false discoveries per family and at most 5%
  of families with any discovery;
- `mixed`: mean FDR at most 0.05 and mean power at least 0.60;
- `strong_ld`: mean FDR at most 0.05 and mean power at least 0.65;
- at least 49 of 50 fully completed families in every confirmatory scenario.

The `ld_mismatch`, `scale_uncertainty`, `undeclared_overlap`, and `weak_paths`
scenarios remain diagnostic. `aligned_sensitivity` remains an explicit
identification-boundary test and cannot validate or invalidate the model.

## Decision scope

Passing establishes family-level calibration under the stated data-generating
model and perturbations. It does not establish validity under exactly aligned
pleiotropy, outcome-driven instrument selection, incorrect allele
harmonization, arbitrary missingness, or an unrepresentative external scale
panel. Competitor benchmarks and a real-data input construction pipeline
remain separate release gates.
