# JG-0.1 adversarial validation results

## Scope

This is step 5 of joint-core development. The plan, seeds, posterior patterns,
and minimum 80% pattern rate were committed before the full run in commit
`b75989a`. The frozen `JG-0.1` priors, quadrature, and posterior thresholds were
not changed after inspecting results.

The run used 20 cells, 50 replicates per cell, and 1,000 datasets total. R
generated the data and orchestrated execution; the C++ evaluator computed all
posteriors. The complete summary and replicate-level results are in
`research/joint_graph_v0_1_adversarial_summary.tsv` and
`research/joint_graph_v0_1_adversarial_summary_replicates.tsv`.

## Prespecified results

| Cell | Correct | High two-path rate | Mean PP XM | Mean PP global MY | Mean PP pleiotropy |
| --- | ---: | ---: | ---: | ---: | ---: |
| Matched null | 50/50 | 0.00 | 0.132 | 0.048 | 0.048 |
| Matched mediation | 50/50 | 1.00 | 1.000 | 1.000 | 0.048 |
| Matched pleiotropy | 50/50 | 0.00 | 1.000 | 0.050 | 1.000 |
| Matched coexistence | 50/50 | 1.00 | 1.000 | 1.000 | 1.000 |
| Moderate mediation | 0/50 | 0.16 | 0.530 | 0.298 | 1.000 |
| Few-B mediation | 49/50 | 1.00 | 1.000 | 1.000 | 0.071 |

The four strong matched controls passed. Eight strong molecular anchors were
sufficient in this particular few-B cell. These results reproduce behavior
under the model's easiest assumptions; they are not general calibration.

The moderate mediation cell failed completely. Its true `a=b=0.40` values lie
between the frozen seven-node quadrature points. The evaluator assigned mean
`PP_nonaligned_P` of 1.00, while only 16% of replicates had high two-path
support. This is numerical model aliasing, not an identification theorem. A
general implementation cannot use this fixed coarse quadrature.

## Information boundaries

With `a=b=0.20` and larger standard errors, no replicate had high two-path
support. Mean `PP_two_path` was 0.006 and the pleiotropy posterior was unstable
(mean 0.436; 24% at least 0.80). `JG-0.1` therefore has poor weak-path behavior
under its present quadrature and priors.

The cell named `no_M_anchor_mediation` removed B and C variants, but it did not
realize the exact mathematical boundary. The generator and fitted model retain
the frozen positive residual molecular variance `vM_A=0.0025`; 60 A variants
therefore still contain modeled `d` variation. The resulting 100% high
two-path rate must not be interpreted as evidence that A-only data solve the
identification problem. It instead shows that inference can depend strongly on
the assumed role variances. The exact zero-`d` boundary remains established by
the algebraic identification result, not this cell.

## Pleiotropy boundaries and violations

| Cell with `b=0` | High two-path rate | Mean PP global MY | Mean PP pleiotropy |
| --- | ---: | ---: | ---: |
| Sparse pleiotropy, `q=0.05` | 0.00 | 0.048 | 0.658 |
| Near-aligned pleiotropy, `q=0.95` | 1.00 | 0.999 | 0.567 |
| Exactly aligned pleiotropy, `q=1` | 1.00 | 1.000 | 0.048 |
| Directional `d` pleiotropy | 1.00 | 1.000 | 1.000 |
| Sparse outcome outliers | 0.00 | 0.048 | 1.000 |
| Role-contaminated pleiotropy | 0.00 | 0.051 | 1.000 |

Exactly aligned pleiotropy was assigned to the global molecular-outcome path
in every replicate. This is the expected observational equivalence and is why
the output scope explicitly excludes exact aligned pleiotropy. Near alignment
also produced high two-path support in every replicate.

The directional `abs(d)` perturbation produced high two-path and high
pleiotropy support in every replicate even though `b=0`. Factorial graph states
prevent M1/M5 mutual exclusion, but they do not by themselves distinguish a
global causal slope from every directional pleiotropic distribution. This is a
serious assumption violation requiring a directional mean component and
sensitivity analysis before production use.

Sparse outcome outliers and 20% A/B role-label contamination activated the
pleiotropy state without creating high two-path evidence in these cells.

## LD, overlap, and scale diagnostics

The simple null sample-overlap cell and AR(1) LD null cell produced no high
two-path calls. One of 50 LD-null replicates had high pleiotropy support. The LD
mediation cell retained 100% high two-path support. These cells do not validate
an independence likelihood under LD: they examine only one correlation
structure and do not test posterior calibration or uncertainty coverage.

Multiplying all latent variances by four caused high pleiotropy support in all
null replicates, although it did not produce high two-path calls. Low-variance
null and mediation cells preserved their coarse classifications. This confirms
that fixed role variances can be mistaken for a biological mechanism and must
be estimated or integrated out.

## Runtime

The 1,000-dataset matrix completed locally with eight workers in approximately
20 seconds of wall time. Median standalone process time was about 0.06 to 0.13
seconds per protein, depending on the number of variants and concurrent load.
No Biowulf login-node computation was used.

## Decision

Step 5 rejects `JG-0.1` as a production inference engine. It establishes a
useful structural prototype and an efficient C++ evaluator, but identifies
three blocking issues:

1. fixed seven-node quadrature aliases ordinary off-grid effect sizes with
   pleiotropy;
2. exact and near-aligned pleiotropy remain nonidentified or weakly separated;
3. directional pleiotropy and misspecified role variances can produce confident
   but incorrect mechanism probabilities.

The next model version should first replace fixed quadrature with adaptive
integration or a validated Laplace/variational calculation. It must then add
estimated role scales, signed-LD and declared sampling covariance, and an
explicit directional-pleiotropy component. The same seeds are development data
and cannot serve as held-out validation of those repairs.
