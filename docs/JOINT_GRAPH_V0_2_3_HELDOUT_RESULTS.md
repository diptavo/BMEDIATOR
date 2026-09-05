# JG-0.2.3 held-out family results

## Provenance and decision

The implementation and plan were frozen at `4771692`. Biowulf build job
`29156019`, array job `29156034`, and summary job `29156111` ran under
`/data/Dutta_lab/BMEDIATOR_runs/JG_0_2_3_4771692`. The source archive SHA-256
was `88533c023165f85ff2d5d08076a63dee4d2dbb7cc4390bf0b8c6aaa701593c94`.
All computation ran on batch nodes.

Exact integration over `q` reduced numerical suppressions substantially, and
all FDR/power criteria passed among reportable fits. The frozen completion
criterion still failed:

| Scenario | Failed proteins | Mean FDR | Mean power |
| --- | ---: | ---: | ---: |
| baseline | 51/5000 | 0.0000 | 0.878 |
| rare | 50/5000 | 0.0000 | 0.590 |
| composite null | 21/5000 | 0.0000 | NA |
| mixed | 30/5000 | 0.0100 | 0.958 |
| strong LD | 55/5000 | 0.0000 | 0.818 |

The remaining diagnostic cells had 31 to 64 suppressed proteins. Only 13 to
33 of 50 families per scenario were fully complete, so `JG-0.2.3` is rejected
as a production engine under its frozen plan.

## Failure analysis

A representative failed null had only one problematic state: posterior mass
`0.00144` and an order-11-to-13 log-evidence change of `0.115`. Its exact
maximum one-state change on the normalized posterior scale was approximately
`0.00018`; the summed bound across all states was `0.00212`. Suppressing the
entire protein because of the relative evidence error in that low-mass state
was unnecessarily conservative.

`JG-0.2.4` therefore retains high-order quadrature but gates the summed
normalized-posterior perturbation at `0.01`, with a separate hard one-log-unit
state cap. No result from this failed run is reused as final validation.
