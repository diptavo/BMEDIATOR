# JG-0.1 adversarial validation plan

## Frozen purpose

This plan defines step 5 of joint-core development before its results are
examined. It tests the frozen `JG-0.1` likelihood and priors without changing
its model, quadrature, posterior thresholds, or identification label.

The experiment has three aims:

1. replicate the four matched strong-signal controls;
2. measure loss of information near weak-instrument and identification
   boundaries;
3. expose false mediation-compatible evidence when assumptions omitted from
   `JG-0.1` are violated.

This is adversarial model validation, not proteome-wide calibration. A failed
cell will be retained and documented. It will not be repaired by changing a
threshold after inspecting the result.

## Execution

The R driver generates summary statistics with fixed seeds and invokes the
standalone C++ posterior evaluator for every dataset. There are 50 replicates
per cell. Cell `j`, replicate `r` uses seed

```text
20262000 + 1000 j + r.
```

The four matched controls use 60 variants in each of roles A, B, and C. Unless
specified below, standard errors are `0.02`, the residual RF-outcome path is
the first nonzero seven-node quadrature point under its prior, and strong path
coefficients are the first nonzero seven-node point under the `a`, `b`, and
`lambda` priors. These choices reproduce the locked implementation fixtures.

## Posterior patterns

The prespecified patterns are:

| Pattern | Required posterior behavior |
| --- | --- |
| `null` | `PP_XM <= 0.20` and `PP_global_MY <= 0.20` |
| `mediation` | `PP_XM >= 0.80`, `PP_global_MY >= 0.80`, and `PP_nonaligned_P <= 0.20` |
| `pleiotropy` | `PP_XM >= 0.80`, `PP_global_MY <= 0.20`, and `PP_nonaligned_P >= 0.80` |
| `coexistence` | all three marginal mechanism probabilities are at least `0.70` |

For cells with a declared pattern, the minimum acceptable correct-pattern rate
is `0.80`. In every cell with `b=0`, a high mediation-compatible claim is
defined independently as `PP_two_path >= 0.80`; its rate is always reported.
It is not assigned a universal 0.05 guarantee because posterior thresholds are
not frequentist tests.

## Frozen cells

| Cell | Truth or perturbation | Status | Pattern tested |
| --- | --- | --- | --- |
| `matched_null` | no `a`, `b`, or pleiotropy | matched | null |
| `matched_mediation` | strong `a` and `b` | matched | mediation |
| `matched_pleiotropy` | strong `a`, `b=0`, `q=0.35` | matched | pleiotropy |
| `matched_coexistence` | strong `a`, `b`, and `q=0.35` pleiotropy | matched | coexistence |
| `moderate_mediation` | `a=b=0.40`, between quadrature nodes | numerical stress | mediation |
| `weak_mediation` | `a=b=0.20`, SE `0.05` | weak information | diagnostic |
| `few_B_mediation` | 60 A, 8 B, no C variants | weak information | mediation |
| `no_M_anchor_mediation` | 60 A, no B or C variants | identification boundary | diagnostic |
| `sparse_pleiotropy` | `b=0`, `q=0.05`, below fitted q grid | misspecified mixture | diagnostic |
| `near_aligned_pleiotropy` | `b=0`, `q=0.95` | near boundary | diagnostic |
| `exact_aligned_pleiotropy` | `b=0`, `q=1` | nonidentified boundary | diagnostic |
| `directional_d_pleiotropy` | outcome shift aligned with `abs(d)` | misspecified mean | diagnostic |
| `sparse_outcome_outliers` | 5% large directional outcome shifts | misspecified tails | diagnostic |
| `overlap_null` | null paths, error correlations 0.50 | undeclared sample overlap | diagnostic |
| `ld_null` | null paths, AR(1) signed LD `rho=0.80` | ignored LD | diagnostic |
| `ld_mediation` | strong mediation, AR(1) LD `rho=0.80` | ignored LD | diagnostic |
| `low_variance_null` | latent role variances multiplied by 0.25 | scale misspecification | diagnostic |
| `high_variance_null` | latent role variances multiplied by 4 | scale misspecification | diagnostic |
| `low_variance_mediation` | mediation with variances multiplied by 0.25 | scale misspecification | diagnostic |
| `role_contamination_pleiotropy` | 20% of A and B labels exchanged | role misspecification | diagnostic |

The aligned cells are not expected to validate a mediation/pleiotropy label.
Their required result is preservation of the explicit output scope
`CONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY`; their posterior behavior measures
how the restricted model aliases the observationally equivalent mechanism.

## Interpretation rules

- A matched-cell failure challenges the frozen model or implementation.
- A numerical-stress failure challenges the seven-node evidence calculation.
- A weak-information failure indicates inadequate information, not necessarily
  model misspecification.
- A boundary cell cannot establish the nonidentified distinction.
- A misspecified cell quantifies why that omitted component must be added
  before production use.
- Results from these seeds cannot be used as held-out confirmation of a later
  repaired model.
