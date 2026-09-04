# Frozen factorized validation

This directory contains the compact summaries from the prespecified
factorized validation run completed on Biowulf on 2026-09-03.

- Configuration: `sim/configs/factorized_final_validation.json`
- Base seed: `20261603`
- Slurm array: `29021148` (3,200 tasks; all completed)
- Summary job: `29021152`
- Prior-sensitivity job: `29021157`
- Protein analyses: 1,360,000
- Binary SHA-256: `bd39c66663ca3d7b993ad4964f89344cff765be17d351df0153875c91cbbe193`
- Configuration SHA-256: `108075810bb87239e3024d571985e12b944f7445fc74ada06332ad5e89d5723d`

`factorized_decision_summary.tsv`, `factorized_posterior_calibration.tsv`,
`factorized_scenario_summary.tsv`, and
`factorized_bayes_prior_sensitivity.tsv` are the release summaries. The
replicate-level file is retained for auditability, but is not required to run
BMEDIATOR.

The run supports the implementation checks described in
`docs/VALIDATION.md`. It does not establish globally calibrated posterior
probabilities or production readiness. The analytically valid Student-t
density-ratio e-BH procedure made no discoveries in the tested families.
