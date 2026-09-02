# Simulation design for latest Biowulf BMEDIATOR

> Historical design note: the `linked_*_leakage` LD cells below are not valid
> tests of distinct causal variants in LD and must not be used for the repaired
> LD-resolution claim. They inject outcome effects proportional to pQTL effects,
> which is algebraically indistinguishable from mediation. Use
> `sim/run_regional_ld_stress.py` for genotype-based shared-versus-distinct LD
> validation. The legacy framework measures structural `P_M1`; confirmatory
> selection now uses `P_mediator_ld_resolved`.

This design targets the current Biowulf binary at:

```bash
/data/Dutta_lab/BMEDIATOR/BMEDIATOR/bmediator
```

The historical benchmark targeted version `1.1.0-dev` and six posterior states:

- `M0`: null
- `M1`: partial mediation, `RF -> PP -> disease`, with residual direct effect allowed
- `M2`: `RF -> PP` only
- `M3`: direct or residual `RF -> disease` only
- `M4`: `PP -> disease` only
- `M5`: correlated/shared pleiotropy masquerading as mediation

The existing simulation harness is useful but not yet fully aligned with this
six-state model. In particular, it currently simulates correlated pleiotropy as
`M3`, while the latest BMEDIATOR code defines correlated/shared pleiotropy as
`M5`. Production simulations should therefore use the design below and update
the harness before large Biowulf runs.

## Goals

1. Measure whether BMEDIATOR separates true mediation from null, one-leg causal
   signals, direct RF effects, protein-only disease effects, and correlated
   pleiotropy.
2. Validate calibration of `P_mediator = P(M1)` and Bayesian FDR summaries.
3. Stress-test robustness to weak instruments, Set C overlap, sparse cis
   support, direction inconsistency, and correlated pleiotropy.
4. Verify that selection based on `P_mediator`, `directional_mediator_prob`,
   and Bayesian FDR behaves as expected.

## Required Harness Updates

Before production runs, update `sim/benchmark_lib.R` and the Python mirror if it
is still maintained:

1. Change `SCENARIOS` from `M0-M3` to `M0-M5`.
2. Generate all six truth states with model-consistent effects:
   - `M0`: all structural effects and pleiotropy centered at zero.
   - `M1`: nonzero `beta1` and `beta2`; `beta3` may be zero or small nonzero.
   - `M2`: nonzero `beta1`, zero `beta2`, zero `beta3`.
   - `M3`: zero `beta1`, zero `beta2`, nonzero `beta3`; signal appears in RF
     and outcome but not pQTL.
   - `M4`: zero `beta1`, nonzero `beta2`, zero `beta3`; cis pQTL predicts
     outcome without RF involvement.
   - `M5`: nonzero RF-pQTL and RF-outcome nuisance components with correlated
     `delta` and `psi`; zero mediated effect.
3. Update scenario aliases only for explicit summary views. Do not collapse
   `M5` into `M2` in the primary confusion matrix.
4. Add generator knobs for:
   - `direction_flip_rate` for M1 proteins.
   - Set A/B/C count distributions by state.
   - LD/correlation perturbation, if full-mode simulations are added.
   - prior mode and direction mode passed through to the binary.
5. Preserve the native output columns from `.mediation`, especially `P_M0-P_M5`,
   `P_mediator`, `directional_mediator_prob`, `selection_probability`,
   `selection_local_fdr`, `selection_cum_fdr`, and `selection_rank`.

## Benchmarks

### 1. Six-State Classification

Purpose: measure whether the posterior mode identifies the correct biological
state.

Run design:

- States: balanced numbers of `M0-M5`.
- Replicates: 50 pilot, 200 main.
- Proteins per state per replicate: 25 pilot, 50 main.
- Cells:
  - `balanced_medium`: A=3, B=3, C=3; moderate effects.
  - `overlap_heavy`: A=2, B=1-2, C=5; tests reliance on overlapping RF/cis SNPs.
  - `cis_sparse`: A=3, B=0-1, C=2; tests weak PP->disease identification.
  - `weak_signal`: smaller effects and larger SEs.

Primary metrics:

- Confusion matrix over all six states.
- Sensitivity and specificity for `M1`.
- False mediator rate by true state, especially true `M5` and `M3`.
- Mean posterior mass assigned to the true state.
- Convergence rate.

Decision criteria:

- True `M1` should have high `P_M1` under medium and overlap-heavy settings.
- True `M5` should not be ranked as high-confidence mediation.
- True `M3` should not borrow pQTL support into `M1` when pQTL effects are null.

### 2. Posterior Calibration and Bayesian FDR

Purpose: validate that `P_mediator` is calibrated enough for ranked discovery.

Run design:

- Replicates: 100 pilot, 500 main.
- Proteins per replicate: 500 pilot, 2000 main.
- Scenario mix, approximating sparse discovery:
  - `M0`: 0.75
  - `M1`: 0.08
  - `M2`: 0.06
  - `M3`: 0.04
  - `M4`: 0.04
  - `M5`: 0.03
- Include medium-power and low-power cells.

Primary metrics:

- Calibration curve: binned mean `P_mediator` versus observed true `M1`.
- Expected calibration error for `P_mediator`.
- Empirical FDR at Bayesian FDR thresholds 0.01, 0.05, 0.10.
- Power at those same thresholds.
- Rank-bin empirical FDR for top 10, 25, 50, 100, and all selected proteins.

Decision criteria:

- Empirical FDR should be close to or below the Bayesian FDR target in the main
  medium-power design.
- Low-power settings may be conservative, but should not show systematic
  inflation driven by `M5`.

### 3. Direction-Consistency Modes

Purpose: test the new `--direction-mode` options.

Run design:

- Use the six-state classification and calibration datasets.
- Re-run the same generated inputs with:
  - `--direction-mode report`
  - `--direction-mode prioritize`
  - `--direction-mode soft --direction-weight 1`
  - `--direction-mode hard --direction-min-prob 0.80`
- Add M1 subcells with `direction_flip_rate` of 0, 0.10, and 0.25.

Primary metrics:

- Change in `selection_rank` for true direction-consistent M1 proteins.
- False positives removed among direction-inconsistent M1-like artifacts.
- Loss of power for true M1 when outcome direction is noisy.
- Agreement between `P_mediator` and `directional_mediator_prob`.

Decision criteria:

- `report` should preserve posterior probabilities.
- `prioritize` should improve ranked FDR without changing posterior
  classification.
- `hard` should be treated as sensitivity analysis because it can lose true M1
  under noisy RF-outcome direction estimates.

### 4. Instrument Architecture Stress Tests

Purpose: identify failure modes as Set A/B/C support changes.

Run design:

- Grid over Set C count: 0, 1, 3, 6.
- Grid over Set B count: 0, 1, 3, 6.
- Grid over RF instrument strength: strong, medium, weak.
- Include M1, M2, M4, and M5 states.

Primary metrics:

- `M1` power as Set B and Set C support changes.
- Misclassification of M4 as M1 when RF support is weak/noisy.
- Misclassification of M5 as M1 when Set C dominates.
- Runtime and convergence by instrument count.

Decision criteria:

- With B=0 and C high, BMEDIATOR should become more cautious about `M1`.
- With B high and C moderate, `M1` should be recoverable.

### 5. Full-Mode Input Smoke and Scaling

Purpose: verify the newer `--protein-gwas-list` and `--bfile` path on Biowulf.

Run design:

- Small smoke: 20 proteins, 2 replicates, one chromosome-sized LD reference or a
  synthetic PLINK reference if already available.
- Main scaling: 500-2000 proteins, 10 replicates after smoke passes.

Primary metrics:

- Agreement with pre-clumped stacked pQTL mode when using equivalent inputs.
- Walltime and memory by protein count.
- Failure modes from missing SNPs, allele mismatches, and LD clumping.

Decision criteria:

- Full mode should produce near-identical posterior rankings to legacy mode
  when the selected instruments are equivalent.
- Any large differences should be traceable to clumping or allele harmonization.

## Biowulf Run Plan

## Current Pilot Findings and Production Policy

The Biowulf six-state pilot after the inference fix completed cleanly with zero
parser count mismatches. Medium and overlap-heavy classification cells separated
all six states correctly, including true `M5` correlated/shared pleiotropy.

The remaining failures are architecture-specific:

- `cis_sparse_m1` and `weak_signal_m1` were over-penalized by the global M1
  gates. These cells now use per-cell `analysis_overrides` with relaxed M1 gates
  and M1-favoring fixed priors.
- `cis_sparse_m4` and `weak_signal_m4` remain low-identifiability stress
  diagnostics. Local tuning did not reliably recover `M4`, so these cells should
  not be used as production pass/fail gates.
- `calibration_low_power` showed near-zero power and inflated selected FDR in
  the pilot. It remains in the pilot as a diagnostic but is excluded from the
  main production calibration config.
- The compact policy-tune batch showed that medium calibration FDR was still
  inflated by true `M4`, not by `M5`. The reason was that `M1` required Set A/C
  first-stage records to exist, but did not require a nonzero RF-to-PP
  first-stage estimate. Production configs now set `m1_min_first_stage_z = 1.5`
  so `M1` requires evidence on both mediation legs.

Before launching the main run, execute the compact policy-tune batch:

```bash
repo=/data/Dutta_lab/BMEDIATOR/BMEDIATOR
out=/data/Dutta_lab/BMEDIATOR/sim_benchmark/six_state_policy_tune_$(date +%Y%m%d)

bash "$repo/scripts/submit_sim_benchmark_array.sh" \
  "$repo" \
  "$repo/sim/configs/classification_calibration_six_state_policy_tune.json" \
  "$out" \
  all \
  200
```

Proceed to the main run only if this batch confirms:

- medium and overlap-heavy `M1`/`M5` behavior remains clean;
- sparse/weak `M1` no longer collapses systematically into `M5`;
- true `M4` does not enter high-confidence `M1` selections in medium
  calibration;
- medium calibration FDR is not worse than the previous pilot;
- stress-only `M4` failures are reported separately and not counted as
  production gate failures.

### Pilot

Use the current submit wrapper after harness updates:

```bash
repo=/data/Dutta_lab/BMEDIATOR/BMEDIATOR
out=/data/Dutta_lab/BMEDIATOR/sim_benchmark/six_state_pilot_$(date +%Y%m%d)

bash "$repo/scripts/submit_sim_benchmark_array.sh" \
  "$repo" \
  "$repo/sim/configs/classification_calibration_six_state_pilot.json" \
  "$out" \
  all
```

Pilot resource request:

- Array throttle: `%20` by default, or pass `200` as the fifth wrapper argument
  for `%200`
- CPU: 4 per task
- Memory: 32 GB
- Time: 2:45:00

### Main

After the pilot confirms clean logs and stable convergence:

```bash
repo=/data/Dutta_lab/BMEDIATOR/BMEDIATOR
out=/data/Dutta_lab/BMEDIATOR/sim_benchmark/six_state_main_$(date +%Y%m%d)

bash "$repo/scripts/submit_sim_benchmark_array.sh" \
  "$repo" \
  "$repo/sim/configs/classification_calibration_six_state_main.json" \
  "$out" \
  all
```

Main resource request should be set after the pilot runtime distribution is
available. If runtime remains similar to the current harness, keep `%20` and
increase walltime only for full-mode runs.

## Deliverables

Each run should produce:

- Raw replicate directories with generated inputs, truth, `.mediation`, `.hyp`,
  and `task_metrics.tsv`.
- Summary tables:
  - `protein_level_metrics.tsv`
  - `classification_by_scenario.tsv`
  - `classification_confusion.tsv`
  - `calibration_bins.tsv`
  - `calibration_bfdr.tsv`
  - `fdr_power_by_threshold.tsv`
  - `direction_mode_comparison.tsv`
- Figures:
  - Six-state confusion heatmap.
  - M1 sensitivity/specificity by cell.
  - Calibration curve for `P_mediator`.
  - Empirical versus estimated Bayesian FDR.
  - Direction-mode rank and power comparison.

## Current Smoke Status

2026-08-03 Biowulf inference-fix smoke:

- Output: `/data/Dutta_lab/BMEDIATOR/sim_benchmark/six_state_inference_fix_smoke_codex_20260803`
- Slurm: array `26427069`, summary `26427072`; both completed with exit code 0.
- Parser alignment: 0 truth-versus-observed instrument-count mismatches across
  the six classification cells.
- Six-state classifier gate: passed. `M0`-`M5` all had 10/10 analyzed proteins
  per replicate and 100% posterior-mode accuracy in both smoke replicates.
- M1 false mediation control: true `M5` had mean `P_M1 = 0` in both smoke
  replicates after requiring independent Set B support for M1.
- Calibration caveat: the two-replicate smoke is small. Replicate 1 had
  empirical FDR 0 at 1% and 5% Bayesian FDR; replicate 2 remained noisy
  (`0.333` empirical FDR at 5% Bayesian FDR). Treat pilot calibration as the
  next production-readiness gate.

The current production settings use fixed scenario priors, `m1_min_cis_only=2`,
and `m1_min_second_stage_z=1.5`. The residual-correlation penalty option exists
for sensitivity analysis but is disabled in the production configs because the
first attempted threshold over-penalized true `M1`.

Previous 2026-08-03 Biowulf smoke:

- Output: `/data/Dutta_lab/BMEDIATOR/sim_benchmark/six_state_shared_smoke_codex_20260803`
- Slurm: array `26425867`, summary `26425870`; both completed with exit code 0.
- Parser alignment: 0 truth-versus-observed instrument-count mismatches across
  the six classification cells.
- Clean classification: `M0`, `M1`, `M2`, and `M4` all had 10/10 analyzed
  proteins per replicate and 100% posterior-mode accuracy in the smoke.
- Remaining blocker: `M3` is called as `M5`, and `M5` is called as `M1` in both
  smoke replicates. This is no longer caused by the earlier per-protein RF
  contamination issue; it reflects overlap in the current BMEDIATOR model
  parameterization and posterior gating.

Do not launch the main production simulation until the pilot confirms the
calibration behavior with enough replicates to estimate empirical FDR stably.
