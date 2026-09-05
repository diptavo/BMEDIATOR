# JG-0.2.4 final held-out family results

## Provenance and decision

The implementation and plan were frozen at `2031cb7`. Biowulf build job
`29157704`, array job `29157710`, and summary job `29157723` ran under
`/data/Dutta_lab/BMEDIATOR_runs/JG_0_2_4_2031cb7`. The source archive SHA-256
was `b576c200677d2b42e61543bcfadc041d2ee13c31cb91361d5ad926a2c4d5c2b2`.
All computation ran on Slurm batch nodes.

All 500 families and 50,000 protein analyses ran. The five frozen statistical
criteria passed among reportable analyses, but the confirmatory completion
criterion failed. `JG-0.2.4` is therefore rejected as a production engine.

| Scenario | Complete families | Failed proteins | Mean FDR | Mean power |
| --- | ---: | ---: | ---: | ---: |
| baseline | 40/50 | 10/5000 | 0.0000 | 0.864 |
| rare | 48/50 | 2/5000 | 0.0000 | 0.590 |
| composite null | 41/50 | 11/5000 | 0.0000 | NA |
| mixed | 38/50 | 15/5000 | 0.0081 | 0.953 |
| strong LD | 46/50 | 4/5000 | 0.0000 | 0.830 |

The frozen criterion required at least 49 of 50 complete families in every
confirmatory scenario. The observed range was 38 to 48. Across all ten
scenarios, 66 of 50,000 protein fits were suppressed and 59 of 500 families
were incomplete.

Diagnostic scenarios had 4 failed fits under LD mismatch, 7 under external
scale error, 4 under undeclared overlap, 6 with weak paths, and 3 in the
exact-alignment sensitivity arm. The exact-alignment arm retained 100% false
discovery as expected from the stated nonidentifiability boundary.

Batch diagnostic job `29159048` replayed all 66 suppressed fits with the
correct frozen seed base. Every suppression was caused by the posterior-error
limit: estimated posterior error ranged from `0.01006` to `0.02879`, while the
maximum relevant state discrepancy ranged from `0.0400` to `0.2030`, below the
one-unit state cap in all 66 cases. The failures included 30 nulls, 16
non-mediation pleiotropic/pathway cases, and 20 true mediation or weak-mediation
cases. They are therefore neither confined to nulls nor ignorable for power.

## Interpretation

The posterior-level stability rule fixed the conceptual problem in the
per-state `JG-0.2.3` gate, but it did not make tensor adaptive Gauss-Hermite
integration sufficiently reliable. Family-wide posterior-FDR selection must
be unavailable when any protein fails.

The next numerical implementation should replace tensor quadrature rather
than relax the frozen reportability threshold. Suitable candidates include
deterministic sparse-grid/cubature integration or a validated bridge-sampling
or importance-sampling estimator with independent error estimates. Any
replacement must be frozen and tested on new seeds; these results cannot be
reused as its final validation.
