# Factorized pilot, 2026-09-03

This frozen pilot used `sim/configs/factorized_identification_smoke.json` with
10 replicates and 20 independent datasets per truth state per cell. There are
200 observations per cell/state. Each protein was run separately with
`--structural-method factorized`; no simulated truth was used by the
estimator or its thresholds.

The identified cell used 6 Set A and 4 Set B instruments. The balanced
pleiotropy cell used 8 Set A and 8 Set B instruments with Set B direct-effect
SD 0.04 relative to outcome SE 0.045. The weak cell used weaker effects and
larger standard errors.

At nominal conjunction p <= 0.05, false two-stage rejection was 0%-0.5% in
the identified cell and 1%-2.5% in the balanced-pleiotropy cell across M2, M4,
and M5. Power was 38.5% and 68%, respectively. The small-sample t reference is
intentionally conservative.

Requiring both fixed-prior leg Bayes factors to exceed 10 gave 81.5% and 85.5%
power in the identified and balanced-pleiotropy cells. Corresponding non-M1
rates were 0%-0.5% and 0.5%-4%. In the weak cell, both procedures had almost no
power, correctly exposing rather than concealing weak identification.

This is implementation evidence, not release-grade validation. It does not
include an LD reference, winner's curse, ancestry mismatch, sample overlap,
or the regional H3/H4 gate. The exact aligned-pleiotropy boundary remains
nonidentifiable and is documented separately.
