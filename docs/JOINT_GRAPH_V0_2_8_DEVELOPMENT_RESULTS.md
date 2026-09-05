# JG-0.2.8 small-scale development results

## Decision

`JG-0.2.8` passed its fresh small-scale numerical-completion and
family-calibration experiment on Biowulf, but it is not production-ready. It
repaired the two JG-0.2.7 sparse-integration failures without changing the
likelihood, priors, or reportability thresholds. A broader historical replay
still left one of 66 difficult cases unreportable, so the complete validation
chain correctly failed closed.

This is development evidence, not a production-release decision. The run used
10 families per scenario and mechanisms generated from the fitted model. A
larger frozen confirmatory experiment, independent input-preparation
validation, real-data replication, and effect-estimation output remain open.

## Frozen run

- Source commit: `01851fc`
- Source archive SHA-256:
  `78c3ba76b49fb6f7e26a8dcfb7247adfba1d4afe0b8ef8073a094141cdaa29ed`
- Biowulf root:
  `/data/Dutta_lab/BMEDIATOR_runs/JG_0_2_8_DEV_01851fc`
- Fresh seed base: `90400000`
- Design: 10 scenarios, 10 families per scenario, 100 proteins per family
- Build/test job: `29185298`
- Historical replay jobs: `29185299` and `29185300`
- Family array: `29185302`
- Summary job: `29185303`
- Sanitizer job: `29185709`

All compilation, tests, simulations, and summaries ran as Slurm jobs on
compute nodes. The Biowulf login node was used only for file transfer and
scheduler/status operations.

## Repair under test

Successive Smolyak levels now cache state-specific tensor components instead
of recomputing shared components. The bounded fallback can continue through
sparse level 15, corresponding to a maximum one-dimensional Gauss-Hermite
order of 31. The fixed normalized-posterior error limit remains `0.01`, the
relevant-state log-evidence difference limit remains one, and unreliable fits
still fail closed.

## Numerical completion

All 10,000 fresh protein fits were reportable. There were no incomplete
families. Across all fits, median runtime was 1.44 seconds, the 95th percentile
was 4.76 seconds, and the maximum was 57.8 seconds. The maximum accepted
posterior-error estimate was `0.00995`.

Nineteen fits invoked sparse-grid integration and eight required
posterior-aware refinement. None of the fresh fits needed sparse level 13 or
higher. High levels were exercised by permanent historical regression cases:

| Historical case | Truth | `PP_two_path` | Posterior error | Sparse level | Seconds |
|---|---|---:|---:|---:|---:|
| JG-0.2.7 mixed replicate 2, P097 | mediation plus sparse pleiotropy | 0.9586 | 0.00869 | 14 | 149.5 |
| JG-0.2.7 strong-LD replicate 9, P080 | null | 0.00065 | 0.00946 | 15 | 8.2 |

Both were previously suppressed at level 12. They now pass the original
thresholds with the expected qualitative inference.

The broader replay made 65 of 66 historical JG-0.2.4 failures reportable. The
remaining LD-mismatch null (replicate 36, P031) stopped at sparse level 15
after 1,754 seconds. Its posterior-error estimate was `0.03744`, so no
posterior was emitted. For comparison, its level-12 JG-0.2.7 error was
`0.04306`; the small improvement, long runtime, and cancellation increase from
`0.999866` to `0.999942` do not justify blindly adding higher sparse levels.
The chained completion job was cancelled by dependency as designed.

The GCC AddressSanitizer and UndefinedBehaviorSanitizer test job completed
successfully, including the level-14 and level-15 regression cases.

## Statistical behavior

Bayesian-FDR results below are averages over only 10 complete families per
scenario, so they are small development estimates rather than precise
performance claims.

| Scenario | Mean Bayesian FDR | Mean power |
|---|---:|---:|
| baseline | 0.000 | 0.840 |
| rare mediators | 0.000 | 0.550 |
| composite null | no discoveries | not applicable |
| mixed | 0.020 | 0.955 |
| strong LD | 0.000 | 0.740 |
| LD mismatch | 0.000 | 0.650 |
| scale uncertainty | 0.000 | 0.830 |
| undeclared overlap | 0.000 | 0.920 |
| weak paths | 0.000 | 0.000 |

The exact-aligned-pleiotropy sensitivity arm produced FDR 1.0, as required by
the stated nonidentifiability boundary. It is not counted as a calibrated
identifiable scenario and cannot be repaired numerically.

## Consequence

The two fresh JG-0.2.7 blockers are repaired, and fresh numerical completion
improved from 9,998/10,000 to 10,000/10,000. However, a severe-cancellation
historical case still blocks numerical release. The next integrator should use
a positive-weight cubature or importance/bridge-sampling fallback with an
independent accuracy assessment, rather than extending this signed sparse
sequence indefinitely. JG-0.2.8 must remain labeled developmental until that
blocker and the larger validation items in `PRODUCTION_READINESS.md` are
complete.
