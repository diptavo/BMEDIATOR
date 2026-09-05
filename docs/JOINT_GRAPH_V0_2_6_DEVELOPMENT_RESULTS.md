# JG-0.2.6 small-scale development results

## Decision

`JG-0.2.6` improves numerical completion and passed one fresh 10,000-protein
development experiment, but it is **not production-ready**. One of the 66
historical `JG-0.2.4` numerical failures remains nonreportable, tensor
quadrature has a severe runtime tail, and the successive-order discrepancy is
a stability diagnostic rather than a proven integration-error bound.

No thresholds were relaxed. Reporting still requires posterior perturbation
at most `0.01`, a relevant-state evidence discrepancy at most one log unit, no
regularized Hessian, and convergence of all 16 states.

## Diagnosed defect and repair

`JG-0.2.4` decided whether to refine using each state's last evidence
difference but decided whether to report using the aggregate perturbation of
the normalized 16-state posterior. Several moderate state errors could
therefore fail the posterior gate without requesting additional integration.

`JG-0.2.5` added posterior-aware refinement: while aggregate error exceeds
`0.01`, the posterior-relevant state with the largest
probability-times-discrepancy contribution advances one quadrature order,
through order 17. `JG-0.2.6` extends this bounded fallback through orders 19
and 21. Every result records the maximum order and number of posterior-aware
refinements. The family runner now retains those diagnostics, and
`research/summarize_joint_graph_numerics.R` summarizes them.

## Provenance

All cluster computation used Slurm batch nodes. No analysis or simulation was
run on the Biowulf login node.

- `JG-0.2.5`: commit `5d8a509`, archive SHA-256
  `bee4e32008b893baa4734762e0cd9a89e097e2021d72e01fa1aa4dc37c96a199`;
  build `29167795`, historical replay `29168036`, family array `29168048`,
  summary `29168274`.
- `JG-0.2.6`: commit `f1d08cd`, archive SHA-256
  `20a407d19dd51196dfd7e758f7534e39be7bd5f077bf75e64c0126ffce0ea12c`;
  build `29171402`, targeted replay `29171456`, fresh family array `29171466`,
  summary `29171467`, historical replay `29172388`.
- Results are under
  `/data/Dutta_lab/BMEDIATOR_runs/JG_0_2_5_DEV_5d8a509` and
  `/data/Dutta_lab/BMEDIATOR_runs/JG_0_2_6_DEV_f1d08cd`.

## Numerical results

`JG-0.2.5` made 62/66 historical failures reportable. Its seed-604 development
run reported 9,997/10,000 proteins; one rare-family null, one composite-null
sparse case, and one scale-uncertainty null remained suppressed.

`JG-0.2.6` repaired all three of those exact cases, with final posterior errors
between `0.00868` and `0.00975`. In the full historical replay it made 65/66
cases reportable. The remaining LD-mismatch null failed at posterior error
`0.01489` after `434.9 s`; extending tensor order further is not justified.

The independent seed-704 experiment reported all 10,000 proteins. Eleven
proteins required posterior-aware refinement. None reached order 19 or 21 in
this seed set. Across all scenarios, median protein runtime was `1.63 s`, the
95th percentile was `9.19 s`, and the maximum was `233.6 s`. Four
composite-null family tasks took 31.5 to 38.6 minutes for 100 proteins, versus
roughly 3 to 5 minutes for most families.

## Calibration and power

The fresh experiment used 10 families of 100 proteins per scenario and seed
base `70400000`. Bayesian-FDR selection used the complete 100-protein family.

| Scenario | Complete | Mean family FDR | Power | Interpretation |
| --- | ---: | ---: | ---: | --- |
| Baseline | 10/10 | 0.000 | 0.78 | Moderate mediation |
| Rare mediators | 10/10 | 0.000 | 0.65 | Two mediators per family |
| Composite null | 10/10 | 0.000 | NA | No discoveries |
| Mixed mechanisms | 10/10 | 0.000 | 0.95 | Mediation plus pleiotropy |
| Strong LD | 10/10 | 0.000 | 0.81 | Within-block LD 0.70 |
| LD mismatch | 10/10 | 0.000 | 0.56 | Diagnostic assumption violation |
| Scale uncertainty | 10/10 | 0.000 | 0.72 | Diagnostic assumption violation |
| Undeclared overlap | 10/10 | 0.000 | 0.88 | Diagnostic assumption violation |
| Weak paths | 10/10 | 0.000 | 0.00 | Important power failure |
| Exact alignment | 10/10 | 1.000 | NA | Expected nonidentifiability |

For context, exact 95% binomial intervals for aggregate moderate-path power
are `0.686-0.857` (baseline), `0.408-0.846` (rare), `0.910-0.976` (mixed), and
`0.719-0.882` (strong LD). Zero composite-null selections among 1,000 proteins
has an exact upper 95% bound of about `0.004` on the per-protein selection
rate; this is not an FDR guarantee for broader data-generating mechanisms.

## Remaining critical deficiencies

1. **Numerical integration remains blocking.** Tensor cost grows exponentially
   with state dimension, one historical case still fails, and the error radius
   is not a theorem-backed upper bound. The next integrator should be a
   validated sparse-grid/cubature or importance/bridge-sampling method with an
   independent accuracy assessment.
2. **Weak-path power is inadequate.** No weak mediator was selected in this
   experiment. This must be characterized over instrument strength, block
   count, sample size, path size, and priors rather than hidden by pooling with
   moderate effects.
3. **Exact aligned pleiotropy is unidentifiable.** Its 100% false-call result is
   expected under the stated model boundary. Software or more LD data cannot
   remove it; interpretation requires external assumptions or designs.
4. **Simulation scope is narrow.** Most scenarios are generated from or close
   to the fitted model. Broader non-aligned pleiotropy, winner's curse,
   role-selection error, ancestry mismatch, and independent simulators remain
   necessary.
5. **No effect-size posterior is reported.** `PP_two_path` is model evidence,
   not an estimate or interval for the mediated effect `a*b`.

The statistical point results justify continued development. They do not
justify routine real-data use, a calibrated-production claim, or a methods
paper performance claim.
