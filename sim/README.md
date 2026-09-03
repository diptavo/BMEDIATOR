# BMEDIATOR simulation benchmarks

## Current LD benchmark

Use the genotype-based regional stress runner for mediation-versus-LD claims:

```bash
python3 sim/run_regional_ld_stress.py \
  --binary ./bmediator \
  --outdir build/regional_ld_stress \
  --replicates 100 \
  --ld 0.3 0.6 0.9
```

The runner generates an explicit PLINK reference, correlated regional
genotypes, marginal protein/outcome associations from shared or distinct causal
variants, and two independent RF instruments. It tests three cases:

- `shared_mediation`: a shared regional causal variant under the mediation data-generating model
- `distinct_ld`: different protein and outcome causal variants in LD
- `same_variant_pleiotropy`: the documented non-identifiable boundary, observationally equivalent to mediation without the exclusion assumption

The old JSON benchmark framework remains useful for six-state structural-model
tests, sample overlap, weak instruments, and contamination. Its historical
`linked_*_leakage` cells are deprecated because adding an outcome effect
proportional to the pQTL effect simulates mediation algebraically rather than a
distinct causal variant in LD. The task runner now refuses those cells instead
of producing a misleading LD benchmark. The old framework also runs legacy
pre-clumped mode and therefore validates `P_M1`, not
`P_mediator_ld_resolved`.

## Dedicated M1-versus-M5 study

The M1/M5 validation uses one independent generated dataset and one BMEDIATOR
invocation per protein. This is required because the outcome effects attached
to RF instruments differ between the M1 and M5 data-generating mechanisms;
placing several such proteins in one synthetic outcome GWAS would mix their RF
instruments and invalidate the comparison.

The study varies the number of protein-specific Set B instruments, second-stage
strength, M5 residual correlation, sample overlap, and Set C burden. Two cells
are labeled `nonidentifiable` and are not scored for classification accuracy:
one omits Set B, and one gives M5 outcome-direct effects exactly proportional
to protein effects. The latter has the same observable distribution as M1 and
is included to verify the theoretical identification boundary.

Local smoke run:

```bash
python3 sim/run_m1_m5_task.py \
  --config sim/configs/m1_m5_identification_smoke.json \
  --cell identified_setb4 \
  --replicate 1 \
  --outdir build/m1_m5_smoke \
  --binary ./bmediator
python3 sim/summarize_m1_m5.py --outdir build/m1_m5_smoke
```

Biowulf study:

```bash
bash scripts/submit_m1_m5_validation.sh
```

The primary outputs are `summary/m1_m5/discrimination.tsv` and
`summary/m1_m5/pair_calibration.tsv`. Pairwise scores normalize M1 and M5
support as `P_M1/(P_M1+P_M5)`; all-state accuracy additionally requires M1 or
M5 to beat the other four structural states.

This directory contains a config-driven simulation harness for BMEDIATOR method
benchmarks:

- classification across null (`M0`), mediated (`M1`), and non-mediated states
- calibration of posterior `P(M1)` and Bayesian FDR summaries
- Nature Genetics-style stress tests against summary-statistic competitors

The scripts generate summary-statistics inputs in the native BMEDIATOR file
format, run the compiled `bmediator` binary, and aggregate benchmark metrics.

## Files

- `benchmark_lib.py`: shared simulation and evaluation logic
- `run_benchmark.py`: generate datasets, run BMEDIATOR, and write summaries
- `run_competitor_benchmark.py`: score completed simulation outputs with
  BMEDIATOR, two-step MR, IVW, Egger, SMR/top-cis, coloc-proxy, and MVMR-style
  baselines
- `summarize_benchmark.py`: rebuild summary tables and figures from raw outputs
- `configs/classification_calibration_default.json`: starter benchmark config
- `configs/classification_calibration_main3.json`: main paper-facing three-state config
- `configs/ng_technical_report_stress_smoke.json`: small local check for the
  technical-report stress suite
- `configs/ng_technical_report_stress_main.json`: Biowulf-scale stress suite
  covering baseline, LD/colocalization confounding, sample overlap, pleiotropy,
  weak instruments, contamination, and Set C ablations

## Quick start

From the `BMEDIATOR` directory:

```bash
python3 sim/run_benchmark.py \
  --config sim/configs/ng_technical_report_stress_smoke.json \
  --benchmark all \
  --outdir sim/output/ng_stress_smoke
python3 sim/run_competitor_benchmark.py \
  --outdir sim/output/ng_stress_smoke \
  --benchmark all
```

Recompute summaries only:

```bash
python3 sim/summarize_benchmark.py \
  --outdir sim/output/ng_stress_smoke \
  --rebuild
python3 sim/run_competitor_benchmark.py \
  --outdir sim/output/ng_stress_smoke \
  --benchmark all
```

## Output structure

Each run creates:

- `classification/<cell>/rep_XXXX/`: generated inputs and BMEDIATOR outputs
- `calibration/<cell>/rep_XXXX/`: generated inputs and BMEDIATOR outputs
- `summary/classification/`: BMEDIATOR six-state classification metrics
- `summary/calibration/`: posterior calibration and Bayesian FDR metrics
- `summary/competitors/`: protein-level competitor scores and method summaries

## Notes

- The generator simulates summary statistics directly rather than individual
  genotypes because BMEDIATOR operates on summary data.
- Optional stress-test knobs include correlated sampling error
  (`sample_overlap`), LD-like leakage (`ld_r`, `linked_*_leakage`), and
  contamination rates for pQTL/RF/outcome effects. These preserve the simulated
  truth state while creating realistic false-positive pressure.
- Set membership is explicit:
  - Set A: RF-significant outside cis
  - Set B: cis-only pQTL
  - Set C: overlapping RF-significant cis instruments
- The primary classifier target for competitor summaries is `M1` versus all
  non-mediation states.
- The smoke config is intentionally small so it can run locally. Use the main
  config for paper-grade Biowulf runs.

## Biowulf

The repository also includes a manifest-driven Biowulf workflow:

```bash
bash scripts/submit_sim_benchmark_array.sh \
  /data/Dutta_lab/BMEDIATOR/BMEDIATOR \
  /data/Dutta_lab/BMEDIATOR/BMEDIATOR/sim/configs/ng_technical_report_stress_main.json \
  /data/Dutta_lab/BMEDIATOR/sim_benchmark/ng_stress_main
```

This submits:

- one array task per `benchmark x cell x replicate`
- `32g` memory and `02:45:00` walltime per simulation task
- a dependent summary job that rebuilds aggregate tables, figures, and
  competitor-method summaries
