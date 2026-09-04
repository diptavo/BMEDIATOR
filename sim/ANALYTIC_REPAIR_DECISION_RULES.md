# Prespecified analytical-repair decision rules

These rules were fixed after the 20-replicate pilot and before reading the
200-replicate confirmatory results.

## Confirmatory family

The primary method is the balanced/InSIDE 2-of-2 p-value with AdaFilter-BH at
`alpha=0.05`. The direct balanced BY result is the arbitrary cross-protein
dependence comparator. Directional BY, the two e-BH tracks, and fixed-prior
posterior FDR remain comparators rather than replacement targets.

Before reading results from the separate p-to-e run, the fixed balanced
p-to-e/e-BH method was designated a secondary dependence-robust candidate. It
uses the frozen equal-weight `kappa={0.10,0.25,0.50,0.75}` calibrator mixture.
It must satisfy the same valid-model and robustness-cell FDR criteria below.
Its power is descriptive and neither its grid nor its threshold may be changed
in response to this run.

After the p-to-e run showed low power and before any balanced Student-t e-value
results were generated, the direct balanced Student-t/e-BH method was frozen as
the next dependence-robust candidate. It uses the existing prespecified shift
grid `{2,4,6}`, scale grid `{2,4,8}`, and equal weights. It must satisfy the
same valid-model and robustness-cell FDR criteria. The broad and narrow power
targets remain 0.80 and 0.50, respectively, and no grid or threshold may be
changed after this point.

After those results showed valid but inadequate power, ordinary BH applied to
the balanced partial-conjunction p-value was frozen as a separate candidate
before any BH result was generated. Its theorem requires independence or PRDS
across protein hypotheses, but it does not require independence between the
two leg p-values. It uses no estimated calibration map, filtering threshold,
or tunable effect grid. The new-seed run uses 200 families in each mixed or
robustness cell and 1,000 families in each pure-null cell. The mixed-cell
criteria remain mean FDP no greater than 0.05, `mean FDP + 1.96*MCSE` no
greater than 0.06, broad power at least 0.80, and narrow power at least 0.50.
For each pure-null cell, the point FDR estimate must not exceed 0.05 and the
one-sided 95% Wilson upper bound must not exceed 0.065. Stress cells retain the
mean-FDP criterion of 0.05. No rule or threshold may be changed after this
paragraph is committed to the validation snapshot.

## Observed frozen balanced-BH result

The new-seed run was executed from commit `dc6eb00` without changing the
criteria above. All 3,600 families completed, for 1.8 million protein analyses.
At `alpha=0.05`, broad and narrow mean FDP were 0.0376 and 0.0298, while mean
power was 0.9720 and 0.5945. Least-favorable-null and global-null FDR were
0.009 and 0, with one-sided 95% Wilson upper bounds 0.0154 and 0.0027.
Balanced heterogeneity, sparse outliers, and dense cross-protein dependence
had mean FDP 0.0043, 0.0082, and 0.0232. Every frozen criterion passed.

This result designates balanced partial-conjunction BH as the leading
statistical candidate only where cross-protein independence or PRDS is
defensible. The dense-dependence cell is a stress test and does not establish a
general dependence theorem. Directional pleiotropy and exactly proportional
M5 remain, respectively, an assumption violation and an identification
boundary. Machine-readable decisions are in
`sim/results/analytic_bh_validation_20260905_v1/summary/analytic_bh_frozen_decision.tsv`.

## Valid-model cells

The following cells evaluate the stated balanced/InSIDE working model:

- `broad_balanced`
- `narrow_balanced`
- `rare_mediation`
- `least_favorable_null`
- `global_null`

For each mixed valid-model cell, the primary requirement is mean
replicate-level FDP no greater than 0.05 and a normal-approximation upper bound
`mean FDP + 1.96*MCSE` no greater than 0.06. Pure-null cells must have an FDR
estimate no greater than 0.05; because their FDP is binary, the corresponding
quantity is the familywise probability of at least one rejection.

Power is descriptive except that `broad_balanced` should retain at least 80%
mean power and `narrow_balanced` at least 50% mean power. No minimum power is
set for the 1% mediation cell.

## Robustness cells

`balanced_heterogeneity`, `sparse_outliers`, and
`dense_cross_protein_dependence` are stress tests outside at least one exact
finite-sample theorem. Passing requires mean FDP no greater than 0.05;
otherwise the failure must define an explicit software/manuscript limitation
and cannot be repaired by fitting a threshold to these results.

## Identification boundaries

`directional_pleiotropy` violates the mean-zero/InSIDE assumption.
`aligned_proportional_boundary` makes the M5 direct effect exactly
proportional to pQTL strength and is likelihood-equivalent to mediation for the
tested summaries. High false selection in these cells is expected. They must
remain visibly labeled as assumption-violating or nonidentifiable and may not
be counted as evidence that the balanced method is generally calibrated.

## LD validation

The genotype-based LD run is descriptive for one-protein raw evidence. Under
`no_second_stage`, the balanced conjunction rejection rate at 0.05 should not
exceed 0.075 in either matched-LD cell. The mismatched-reference cell is a
failure-mode audit and is not permitted to justify a positive claim.
