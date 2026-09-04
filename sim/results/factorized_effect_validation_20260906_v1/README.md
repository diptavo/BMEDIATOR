# Frozen post-selection effect validation

This directory contains the compact summaries from the held-out effect-
estimation validation frozen in `sim/configs/factorized_effect_validation.json`
and `sim/EFFECT_VALIDATION_DECISION_RULES.md`.

## Provenance

- Code commit: `bc5d8f1683bda4f243d604addc27270c117f8889`
- Simulation seed: `20260906`
- Biowulf array: `29101196`
- Biowulf summary/evaluation job: `29101198`
- Design: 2,700 simulated families of 500 proteins, or 1.35 million protein
  analyses
- Frozen decision: `YES` (all 87 component criteria passed)

## Confirmatory results

At a target FCR of 0.05, mean replicate false coverage proportions were:

| Cell | BH | BY |
| --- | ---: | ---: |
| Strong independent instruments | 0.01157 | 0.00027 |
| Moderate independent instruments | 0.00165 | descriptive |
| Four independent instruments | 0.00000 | descriptive |
| Balanced heterogeneity | 0.02801 | 0.00114 |
| Dense cross-protein dependence | descriptive | 0.00033 |

Every selected family in every cell had a complete effect set, and no selected
set was unbounded. All confirmatory cells made discoveries in every replicate.
For true-M1 proteins in the four independent confirmatory cells, all three
effect estimates were available in 100% of analyses. The lowest nominal
coverage for either causal leg was 0.94487, and the largest absolute
standardized bias among BH-selected M1 effects was 0.18997. When sample overlap
was declared, all 39,923 tested proteins failed closed and no finite selected
confidence set was emitted.

## Limits

These results validate finite-grid operating characteristics under the frozen
assumptions; they do not prove those assumptions in real data. In particular:

- Under directional pleiotropy, which violates balanced/InSIDE, BH mean FCP
  was 0.29491. This is a documented failure mode, not part of the confirmatory
  claim.
- The near-threshold weak-instrument cell made no selections, so it cannot
  validate weak-instrument effect coverage.
- BY made calls in only 80 of 300 moderate-instrument families and 5 of 300
  four-instrument families; those BY cells were prespecified as descriptive.
- Same-sample instrument selection produced only 25 BH-selected families, and
  unreported sample overlap remains an assumption-violation diagnostic.
- The FCR result is conditional on the marginal finite-instrument Student-t
  profile intervals having their stated coverage. The intervals are not
  weak-instrument-exact.

`summary/effect_validation_frozen_decision.tsv` is the authoritative frozen
pass/fail record. The other files preserve effect, FCR, scenario, and
replicate-level summaries needed to audit that decision.
