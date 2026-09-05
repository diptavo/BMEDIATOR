# JG-0.2.1 held-out family results

## Provenance

The frozen family plan was committed at `9b669cf`; explicit Biowulf R-module
loading was added at `3beb163`. The analysis ran as Slurm build job `29152927`,
array job `29152929`, and summary job `29152941` under
`/data/Dutta_lab/BMEDIATOR_runs/JG_0_2_1_3beb163`. The transferred source
archive SHA-256 was
`32abe80a6434452962f74dfb514c935322d04535763940ad0cf7f620be1f7121`.
No analysis was run on the Biowulf login node.

All 500 families and 50,000 protein analyses ran. The frozen confirmatory
criteria all passed:

| Scenario | Complete families | Mean FDR | Mean power |
| --- | ---: | ---: | ---: |
| baseline | 50/50 | 0.0000 | 0.874 |
| rare | 50/50 | 0.0000 | 0.710 |
| composite null | 50/50 | 0.0000 | NA |
| mixed | 50/50 | 0.0052 | 0.967 |
| strong LD | 49/50 | 0.0000 | 0.794 |

The diagnostic scenarios gave power `0.680` under LD mismatch, `0.806` under
external-scale error, `0.862` under undeclared overlap, and `0.000` for weak
paths. Exactly aligned pleiotropy was selected in every family, as expected
from the stated nonidentifiability boundary.

## Numerical suppressions

Four of 50,000 fits were suppressed by the adaptive-versus-Laplace
one-log-unit diagnostic: one strong-LD mediator and three nulls. Their
differences were `1.022`, `1.089`, `1.004`, and `1.039`. Relaxing the gate only
for diagnosis produced the expected classifications, but the rule was not
changed post hoc.

This result passed its frozen gate but exposed that disagreement with Laplace
was not a direct test of adaptive-quadrature convergence. It therefore became
development evidence for the numerical-only `JG-0.2.2` repair. It is not the
final release calibration.
