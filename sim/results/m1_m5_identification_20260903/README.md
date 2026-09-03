# M1 versus M5 identification study: 2026-09-03

This directory contains the aggregate results of 36,000 independent protein
simulations generated with `m1_m5_identification_main.json`. The design used 50
replicates, 40 M1 and 40 M5 proteins per replicate, and nine design cells. Each
protein was generated and analyzed in a separate dataset to prevent
protein-specific RF instruments from entering another protein's model.

The executable used the production Set A/B/C structural model with fixed
default priors. `discrimination.tsv` reports both conditional M1-versus-M5
metrics based on `P_M1/(P_M1+P_M5)` and the more stringent six-state argmax
metrics. `pair_calibration.tsv` compares the pairwise score with the observed
M1 fraction in ten bins.

## Main findings

- Across the seven identifiable or assumption-stress cells, pairwise AUC was
  0.911 and pairwise accuracy was 0.834, but six-state accuracy was only 0.492.
- Six-state M1 recall was 0.662 and six-state M5 recall was 0.323. M5 was often
  assigned to M2 or M0 even when `P_M5 > P_M1`.
- Increasing Set B instruments improved pairwise AUC from 0.908 with two to
  0.999 with eight, but M5 six-state recall remained approximately 0.39.
- In the weaker second-stage cell, pairwise accuracy was 0.721 and six-state
  accuracy was 0.263.
- In the Set C-heavy cell, pairwise AUC was 0.661, pair scores were defined for
  only 62.1% of observations, and six-state M5 recall was 0.036.
- At the exactly aligned-pleiotropy boundary, 96.7% of M5 datasets were called
  M1 in the pairwise comparison. This is expected because this data-generating
  model is observationally equivalent to M1.
- Pairwise probabilities were not calibrated. Aggregate ten-bin expected
  calibration error was 0.089, with substantial overconfidence in several
  middle and upper bins.

These results fail a production release criterion for reliable M1/M5
classification. They show that multiple Set B instruments provide ranking
information under the model restrictions, but they do not validate the current
`P_M1` and `P_M5` values as calibrated probabilities.
