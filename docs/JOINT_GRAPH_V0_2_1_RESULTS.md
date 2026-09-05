# JG-0.2.1 patch-validation results

## Frozen run

The `JG-0.2.1` implementation and new-seed plan were frozen at commit
`14fe583` before this run. The matrix contained 21 cells, 50 replicates per
cell, and 1,050 protein analyses. All calculations used the C++ evaluator; R
generated matched summary statistics and orchestrated local workers.

All 1,050 analyses completed all 16 graph states without Hessian
regularization. The maximum relevant-state adaptive/Laplace discrepancy was
`0.793`, below the frozen one-log-unit failure threshold. Large discrepancies
in unsupported graph states were retained diagnostically but did not suppress
stable posterior-supported states.

## Confirmatory cells

| Cell | Correct | Rate | High two-path | Mean sparse PP | Mean directional PP |
| --- | ---: | ---: | ---: | ---: | ---: |
| Matched null | 50/50 | 1.00 | 0.00 | 0.030 | 0.004 |
| Moderate off-grid mediation | 44/50 | 0.88 | 0.88 | 0.117 | 0.003 |
| Sparse pleiotropy | 41/50 | 0.82 | 0.00 | 0.863 | 0.005 |
| Directional pleiotropy | 49/50 | 0.98 | 0.00 | 0.045 | 1.000 |
| Sparse plus directional pleiotropy | 41/50 | 0.82 | 0.02 | 0.831 | 1.000 |
| Mediation plus sparse, 30 blocks | 45/50 | 0.90 | 0.96 | 0.949 | 0.006 |
| Declared overlap null | 50/50 | 1.00 | 0.00 | 0.032 | 0.011 |
| Signed-LD null | 50/50 | 1.00 | 0.00 | 0.029 | 0.003 |
| Signed-LD mediation | 42/50 | 0.84 | 0.84 | 0.100 | 0.005 |
| Fourfold scale, declared null | 50/50 | 1.00 | 0.00 | 0.026 | 0.005 |
| 70%-accurate orientations | 50/50 | 1.00 | 0.00 | 0.036 | 1.000 |
| Skeptical-prior mediation | 44/50 | 0.88 | 0.88 | 0.100 | 0.002 |
| Diffuse-prior null | 50/50 | 1.00 | 0.00 | 0.046 | 0.006 |

This repairs the specific `JG-0.1` off-grid aliasing, the original `JG-0.2`
orientation failure, and the fail-closed instability caused by irrelevant
states. The paired default and skeptical prior cells both classified the same
moderate-mediation datasets in 44/50 replicates. The paired default and diffuse
null cells both passed 50/50.

## Information and assumption diagnostics

- Weak mediation produced high two-path support in 8% of replicates. This is a
  low-power regime and must normally be reported as unresolved.
- Twenty-block mediation-plus-sparse truth had high two-path support in 78%; 30
  blocks raised joint-pattern recovery to 90%. Coexistence requires more
  independent signals than either mechanism alone.
- Six-block sparse pleiotropy had mean sparse PP `0.617`, demonstrating weak
  mixture information without creating high two-path calls.
- Correctly declared quarter-scale mediation had high two-path support in 84%.
- Undeclared sampling overlap did not create high two-path calls in this one
  null configuration. This does not justify omitting overlap covariance.
- Fourfold variance misspecification produced mean sparse PP `0.571`, while
  correctly supplied fourfold scales passed the null in all replicates.

Near-aligned and exactly aligned pleiotropy produced high two-path support in
44% and 92% of replicates. This is not repaired because exact alignment is
observationally nonidentified; every result continues to carry the explicit
exclusion condition.

## Runtime

The complete local matrix used eight workers. Median standalone C++ time was
approximately 1.2 seconds per 40-variant protein, 1.7 seconds for 30-block
coexistence, and 3.6 seconds when uncertain signs were exactly integrated. No
Biowulf login-node computation was performed.

## Decision

`JG-0.2.1` passes its frozen patch-development matrix. It is a workable joint
model prototype, but this is not yet a production-release result. Remaining
release gates are:

1. a substantially larger held-out family-level calibration with independent
   seeds and frozen selection rules;
2. LD-reference mismatch, winner's curse, missing-trait, and scale-estimation
   uncertainty;
3. head-to-head comparison with the frozen factorized method and external
   colocalization/MR competitors;
4. integration into the main executable only after those gates pass.

The complete summary and replicate results are stored in
`research/joint_graph_v0_2_1_development_summary.tsv` and its
`_replicates.tsv` companion.
