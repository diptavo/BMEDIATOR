# JG-0.2.7 small-scale development results

## Decision

`JG-0.2.7` is not production-ready. It replaced the high-order tensor fallback
with deterministic Smolyak sparse-grid integration and substantially improved
the earlier pathological regressions, but two of 10,000 fresh analyses still
failed the unchanged numerical-reportability gate.

## Frozen run

- Source commit: `fda164d`
- Source archive SHA-256:
  `ca6658126d02a16dbb34a97819e29559bb1cb048a6202935bc5af7c3443fc94b`
- Biowulf root:
  `/data/Dutta_lab/BMEDIATOR_runs/JG_0_2_7_DEV_fda164d`
- Fresh seed base: `80400000`
- Design: 10 scenarios, 10 families per scenario, 100 proteins per family
- Build/test job: `29183862`
- Family array: `29183868`
- Summary job: `29183870`

The build and complete test suite passed on a compute node. All 100 family
tasks exited normally, but the runner retained two per-protein numerical
failures and therefore withheld family-wide selection for those families.

## Numerical completion

Across 10,000 proteins, 9,998 were reportable and two were suppressed:

| Scenario | Protein | Truth | Posterior error | Final sparse level | Elapsed seconds |
|---|---|---|---:|---:|---:|
| mixed, replicate 2 | P097 | mediation plus sparse pleiotropy | 0.01352 | 12 | 1023.7 |
| strong LD, replicate 9 | P080 | null | 0.01260 | 12 | 649.2 |

Both exceeded the fixed posterior-error maximum of `0.01`; no posterior or
selection was reported for either. Mixed and strong-LD completion was 9/10
families; every other scenario completed 10/10. Among successful proteins,
median runtime was 1.44 seconds, the 95th percentile was 4.85 seconds, and the
maximum was 47.8 seconds. The two suppressed fits created a much larger tail.

## Statistical behavior

The following are development estimates and are conditional on complete
families; they do not override the failed completion decision.

| Scenario | Bayesian FDR | Power |
|---|---:|---:|
| baseline | 0.000 | 0.920 |
| rare mediators | 0.000 | 0.800 |
| composite null | no discoveries | not applicable |
| mixed | 0.017 | 0.939 |
| strong LD | 0.000 | 0.789 |
| LD mismatch | 0.000 | 0.690 |
| scale uncertainty | 0.000 | 0.890 |
| undeclared overlap | 0.000 | 0.960 |
| weak paths | 0.000 | 0.000 |

The exact-aligned-pleiotropy sensitivity arm remained nonidentified, with FDR
1.0 as expected from the stated identification boundary.

## Consequence

Inspection showed that successive nonnested sparse levels recomputed lower
level tensor components. JG-0.2.8 caches those components and extends the
bounded sparse sequence through level 15. The two failures above are permanent
regression fixtures for that repair.
