# BMEDIATOR: Bayesian Mediation MR

**Bayesian framework for evaluating candidate plasma-protein mediators between risk factors and disease outcomes using GWAS summary statistics.**

Version 1.2.0-dev 

> **Development status:** The dedicated 36,000-run M1/M5 study found useful
> pairwise ranking but inadequate six-state classification (M5 recall 32.3%)
> and probability miscalibration. This version is for research evaluation, not
> production-calibrated mediation claims. See
> [Validation Status](docs/VALIDATION.md).

> **Current factorized validation:** Balanced partial-conjunction BH passed all
> frozen criteria in a new-seed 1.8-million-analysis run: broad and narrow
> power were 0.972 and 0.595, with mean FDP 0.0376 and 0.0298. The
> least-favorable and global-null family error rates were 0.009 and 0.000 over
> 1,000 replicates each. This rule requires valid balanced/InSIDE leg tests and
> independence or PRDS across protein-level p-values. Directional pleiotropy,
> exact proportional pleiotropy, and LD-reference mismatch remain explicit
> limitations. Version 1.2.0-dev is not yet a production or publication-ready
> release.

> **Experimental repair:** `--structural-method factorized` separates the two
> causal legs, residual RF effect, and pleiotropy instead of forcing them into
> mutually exclusive M1-M5 states. It reports nuisance-aware fixed-prior
> Bayes factors, model-conditional conjunction p-values, BY q-values,
> safe e-values with e-BH, fixed-prior posterior expected-FDR values, and a
> separate regional H3/H4 identification gate.
> It is not yet the default. See
> [Factorized two-stage method](docs/FACTORIZED_METHOD.md).

> **Joint-core research prototype:** `JG-0.1` is the first frozen reference
> implementation of a replacement joint Bayesian core. It fits one trivariate
> RF/protein/outcome likelihood over Sets A, B, and C and compares eight
> factorial pathway/pleiotropy states. Independent base-R and C++
> implementations agree on locked development fixtures. This prototype is not
> connected to the production command-line analysis and is not calibrated for
> real data. See [JG-0.1 specification](docs/JOINT_GRAPH_MODEL_V0_1.md) and
> [development results](docs/JOINT_GRAPH_V0_1_RESULTS.md).

---

## Overview

BMEDIATOR evaluates whether intermediate phenotypes (for example, plasma
proteins) are compatible with mediating a risk-factor effect on a disease
outcome, using summary statistics from three sources:

1. **RF GWAS** — genome-wide significant instruments for the risk factor
2. **pQTL study** — protein quantitative trait loci (cis and trans)
3. **Cancer GWAS** — disease outcome summary statistics

### The Model

The structural causal model is:

```
RF ---β₁---> PP ---β₂---> Cancer
 \                         ^
  \----------β₃-----------/
         (direct effect)
```

In legacy mode, BMEDIATOR computes ELBO-based variational approximate model
probabilities for six scenarios. These states are retained for compatibility;
the recommended research workflow in this development release is the
factorized mode described below.

| Scenario | Interpretation | Constraints |
|----------|---------------|-------------|
| M=0 | Null | β₁=0, β₂=0, β₃=0 |
| M=1 | **Partial mediation** | β₁≠0, β₂≠0, β₃ allowed |
| M=2 | RF→PP only | β₁≠0, β₂=0, β₃=0 |
| M=3 | RF→disease direct/residual only | β₁=0, β₂=0, β₃≠0 |
| M=4 | PP→disease only | β₁=0, β₂≠0, β₃=0 |
| M=5 | Correlated/shared pleiotropy masquerading as mediation | β₁≠0, β₂=0, β₃≠0, Cov(δ,ψ)≠0 |

### Key Innovation: Three-Set Instrument Partition

For each protein, instruments are partitioned into:

- **Set A** (RF-only): Genome-wide significant for RF, NOT in the cis-region of the protein gene
- **Set B** (cis-only): Significant cis-pQTL for the protein, NOT genome-wide significant for RF
- **Set C** (overlap): Both RF-significant AND in the cis-region

In factorized mode, Set A and Set B are also cross-clumped at `--clump-r2`.
Any cis-pQTL above that LD threshold with an RF instrument is moved to Set C,
and every RF instrument physically inside the cis window is excluded from Set
A regardless of its protein-association p-value.

Set C instruments are excluded from the factorized confirmatory leg estimates
and retained for regional diagnostics; they cannot by themselves distinguish
mediation from correlated pleiotropy. The identifying evidence for β₂ comes
from multiple independent Set B instruments under the assumption that their
outcome-direct effects are not systematically aligned with their protein
effects. The legacy M5 state permits
correlated RF-to-protein and RF-to-outcome residuals at RF-associated Set A/C
instruments; it does not permit unrestricted protein-effect/outcome-direct
covariance at Set B instruments. Without this restriction, M1 and M5 are
observationally equivalent. See
[Identification and LD Resolution](docs/IDENTIFICATION.md).

### Inference

- **Inner loop**: Coordinate Ascent Variational Inference (CAVI) with spike-and-slab priors on pleiotropy terms
- **Default priors**: Prespecified scenario and effect priors remain fixed, avoiding reuse of the analyzed outcomes to construct their own model priors
- **Optional exploratory mode**: `--empirical-bayes` estimates hyperparameters across analyzed proteins
- **Regional check**: Full mode uses the reference-panel LD matrix to identify conditionally independent protein and outcome signals, then compares H0-H4 for every signal pair

`P_M1` is approximate structural support for mediation within the six-state model.
`P_mediator_ld_resolved` is nonzero only when full regional data support a
shared signal, at least two independent RF-to-protein instruments exist, and
at least two independently matched cis protein-outcome signals are available.
This remains conditional on exclusion and
valid-instrument assumptions; see
[Identification and LD Resolution](docs/IDENTIFICATION.md) and
[Validation Status](docs/VALIDATION.md).

The experimental factorized analysis is enabled with:

```bash
./bmediator [the usual input options] --structural-method factorized
```

It requires at least three exact observed Set A instruments and three independent Set B
instruments for a confirmatory call by default. The key new output columns are
`factor_log_BF_XM`, `factor_log_BF_MY`, `factor_conjunction_p`,
`factor_conjunction_q_BY`, `factor_log_e_mediation`, `factor_e_q_EBH`,
`factor_PP_two_stage`, `factor_posterior_cum_fdr`,
`factor_log_BF_directional_XM`, `factor_log_BF_directional_MY`,
`factor_directional_collinearity_XM`, `factor_directional_collinearity_MY`,
`factor_balanced_conjunction_p`, `factor_balanced_conjunction_q_BH`,
`factor_balanced_conjunction_q_AdaFilter`,
`factor_log_e_XM_balanced`, `factor_log_e_MY_balanced`,
`factor_log_e_mediation_balanced`, `factor_e_q_balanced_EBH`,
`factor_log_e_p2e_balanced_mediation`, `factor_e_q_p2e_balanced_EBH`,
`factor_posterior_status`, `factor_mediation_status`,
`factor_frequentist_status`, and `factor_ebh_status`. The factorized posterior
is a fixed-prior working-model probability for two nonzero slopes, not an
assumption-free M1/M5 probability. The e-BH guarantee concerns the two-leg evidence family; the
regional H3/H4 result is an additional causal-identification assessment.
When sampling-error overlap is declared for either causal leg, factorized BY
and e-BH q-values are unavailable rather than being presented as calibrated.
When slope and oriented-intercept curvature is unstable, statistical evidence
is retained but effect estimates are unavailable and the causal status is
`UNRESOLVED_EFFECT_ESTIMATION`; BMEDIATOR does not silently substitute a
zero-intercept estimate.

Factorized mode skips six-state CAVI for efficiency. The legacy `P_M0`-`P_M5`,
ELBO, and legacy selection columns are therefore `nan` in a factorized run;
they are not silently reused as factorized evidence.

A manuscript-style description is available in
[Introduction and Methods](manuscript/BMEDIATOR_METHODS_DRAFT.md).

---

## Installation

```bash
git clone git@github.com:diptavo/BMEDIATOR.git
cd BMEDIATOR
make
```

Requirements: a C++17 compiler and `make`.

OpenMP is optional. It is disabled by default for portability:

```bash
make clean
make USE_OPENMP=1
```

On Linux with GCC, `USE_OPENMP=1` uses the compiler's native `-fopenmp`
support. Apple Clang requires the separate OpenMP runtime:

```bash
brew install libomp
make clean
make USE_OPENMP=1
```

The Makefile detects Homebrew `libomp` automatically. For another installation
location, set `OPENMP_ROOT=/path/to/libomp`. A normal serial build does not
require OpenMP.

---

## Test

Run the bundled smoke test:

```bash
make test
```

The test uses the small synthetic files in `testdata/` and writes temporary
outputs under `build/test/`.

---

## Bundled Reference Resources

Small protein reference resources are included under
`resources/protein_references/`:

- ARIC SeqId resource, pilot resource, and chunk manifest
- UK Biobank PPP EUR and combined OID resources, pilot resources, and chunk manifests
- deCODE SeqId-to-gene/UniProt/TSS mapping files and mapping provenance
- EA protein-info, per-protein manifest, and compact clumped pQTL instrument reference
- tiny debug fixtures for parser and manifest checks

Some manifests preserve Biowulf filesystem paths in `sumstat_file` or
`sumstat_dir` columns. Rewrite those paths for your local or cluster filesystem
before running full analyses. Raw GWAS, raw pQTL, PLINK reference panels, and
generated BMEDIATOR outputs are intentionally excluded; see
`docs/DATA_POLICY.md`.

---

## Running BMEDIATOR

### One protein

A full-mode analysis of one protein uses a manifest containing one row:

```text
PROTEIN_ID    /path/to/protein_sumstats.tsv
```

Run:

```bash
./bmediator \
    --rf-sumstat RF_sumstats.tsv \
    --protein-gwas-list one_protein_manifest.txt \
    --cancer-sumstat outcome_sumstats.tsv \
    --protein-info protein_info.tsv \
    --bfile LD_reference_prefix \
    --structural-method factorized \
    --out one_protein_result
```

### Protein manifest

To analyze multiple proteins, put one protein and its summary-statistics path on
each manifest row and run the same command:

```bash
./bmediator \
    --rf-sumstat RF_sumstats.tsv \
    --protein-gwas-list protein_manifest.txt \
    --cancer-sumstat outcome_sumstats.tsv \
    --protein-info protein_info.tsv \
    --bfile LD_reference_prefix \
    --structural-method factorized \
    --out manifest_result
```

Full mode requires the `.bed`, `.bim`, and `.fam` files represented by
`LD_reference_prefix`. All inputs must use the same genome build.

For three-sample instrument selection, add `--factor-independent-selection`
and provide `P_SELECT` (or `SELECT_P`) in the RF and protein files. `P_SELECT`
must come from a discovery sample independent of the effect estimates used in
the causal-leg likelihoods. Without that flag, selection uses `P` and the
output records `factor_selection_design=same-sample`.

---

## Input File Formats

All input files are tab/space-delimited with headers.

### RF GWAS Summary Statistics (`--rf-sumstat`)

```
SNP     A1  A2  FREQ    BETA      SE        P         CHR  BP
rs12345 A   G   0.30    0.052     0.008     1.2e-10   1    12345678
```

### Cancer GWAS Summary Statistics (`--cancer-sumstat`)

Same format as RF GWAS.

### pQTL Summary Statistics (`--pqtl-sumstat`)

Multi-protein format with an additional PROTEIN column:

```
PROTEIN  SNP      A1  A2  FREQ  BETA      SE        P         CHR  BP
IL6      rs12345  A   G   0.30  0.105     0.012     2.1e-18   1    12345678
IL6      rs67890  T   C   0.15  0.083     0.015     3.4e-8    1    12400000
TNF      rs11111  G   A   0.25  0.072     0.010     5.0e-13   6    31543000
```

**Important**: This file should contain pre-clumped instruments (r² < 0.1). All SNPs for all proteins are included in a single file.

BMEDIATOR also supports a full mode with a per-protein GWAS manifest. In that
mode, use `--protein-gwas-list <manifest>` instead of `--pqtl-sumstat <file>`.
Each non-comment manifest line contains a protein identifier and a path:

```
IL6   /path/to/IL6.PHENO1.glm.linear
TNF   /path/to/TNF.PHENO1.glm.linear
```

The per-protein files use PLINK 2 `glm` tabular columns: `#CHROM`, `POS`, `ID`,
`REF`, `ALT`, `A1`, `A1_FREQ`, `TEST`, `BETA`, `SE`, and `P`; `ERRCODE` is
optional. Only `TEST=ADD` rows are read. Full mode requires unpruned regional
variants, not only lead pQTLs, plus an ancestry-matched `--bfile` reference.

### Protein Annotation (`--protein-info`)

```
PROTEIN  GENE   CHR  START      END
IL6      IL6    7    22725889   22732002
TNF      TNF    6    31543344   31546112
```

---

## Output Files

### `.mediation` — Main Results

Tab-delimited, one row per protein. Legacy mode sorts by
`selection_probability`; factorized mode sorts by finite
`factor_posterior_rank` ascending.

| Column | Description |
|--------|-------------|
| Protein | Protein identifier |
| Gene | Gene name |
| nA, nB, nC | Number of instruments in each set |
| P_M0 ... P_M5 | ELBO-based variational approximate scenario probabilities |
| P_mediator | P(M1), approximate structural support for partial mediation |
| P_mediator_ld_resolved | P(M1) when regional and instrument requirements pass; otherwise 0 |
| P_mediator_identified | Deprecated compatibility alias for P_mediator_ld_resolved |
| P_protein_disease | P(M1)+P(M4), posterior support that the protein affects disease |
| P_rf_responsive | P(M1)+P(M2)+P(M5), posterior support for RF→PP involvement |
| P_rf_direct | P(M1)+P(M3)+P(M5), posterior support for residual/direct RF→disease involvement |
| beta1, se_beta1 | RF→PP causal effect (estimated under M=1) |
| beta2, se_beta2 | PP→Cancer causal effect (estimated under M=1) |
| beta3, se_beta3 | residual/direct RF→Cancer effect (estimated under M=1) |
| ivw_rf_to_outcome_* | IVW RF→disease summary using RF instruments |
| indirect_direction | Sign of β₁×β₂ |
| rf_to_outcome_direction | Sign of the RF→disease IVW estimate |
| direction_consistent | Whether the point estimate signs agree |
| direction_consistency_prob | Posterior-normal approximation to Pr(sign(β₁×β₂)=sign(RF→disease)) |
| proportion_mediated | (β₁×β₂)/(RF→disease IVW estimate); use as a sanity-check, not a strict bound |
| directional_mediator_prob | P(M1) × direction_consistency_prob |
| selection_probability | Probability/score used for mode-specific selection |
| selection_local_fdr, selection_cum_fdr, selection_rank | Mode-specific selection diagnostics |
| regional_PP_shared, regional_PP_distinct | Posterior probabilities of shared and distinct regional causal configurations |
| regional_shared_given_both | P(shared configuration given that both traits are regionally associated) |
| regional_method | Regional inference method (`ld-multisignal` by default) |
| regional_protein_signals, regional_outcome_signals | Number of conditionally independent signals detected for each trait |
| regional_signal_pairs | Number of protein-outcome signal pairs tested |
| regional_max_credible_set_pair_r2 | Maximum cross-trait credible-set LD for the selected signal pair |
| mediation_identifiability | Explicit LD-resolution and conditional-identification state |
| factor_beta1, factor_beta1_se, factor_p_XM, factor_log_BF_XM | Factorized RF-to-protein adjusted-profile estimate and oriented-intercept-adjusted evidence |
| factor_beta2, factor_beta2_se, factor_p_MY, factor_log_BF_MY | Factorized protein-to-outcome adjusted-profile estimate and oriented-intercept-adjusted evidence from Set B only |
| factor_beta1_ci_lower/upper, factor_beta2_ci_lower/upper | Approximate finite-instrument t intervals with `n-2` degrees of freedom |
| factor_fcr_alpha_BH, factor_beta1/beta2/indirect_fcr_ci_*_BH | Selected-protein confidence sets targeting 5% FCR under independent protein vectors |
| factor_fcr_alpha_BY, factor_beta1/beta2/indirect_fcr_ci_*_BY | Harmonic-adjusted selected-protein confidence sets for arbitrary cross-protein dependence |
| factor_beta3, factor_beta3_se, factor_p_XY, factor_log_BF_XY | Factorized residual RF-to-outcome estimate and evidence |
| factor_log_e_XM, factor_log_e_MY, factor_log_e_XY | Safe log e-values from a fixed Student t density-ratio mixture after eliminating common scale and oriented-intercept nuisance terms |
| factor_tau_XM, factor_tau_MY, factor_tau_XY | Profile estimates of residual heterogeneity SD for each leg |
| factor_indirect, factor_indirect_se | Factorized mediated effect and second-order product SE |
| factor_indirect_ci_lower/upper | Conservative simultaneous-confidence-rectangle interval for the product |
| factor_conjunction_p, factor_conjunction_q_BY | Intersection-union p-value and proteome-wide BY q-value |
| factor_min_log_BF | Smaller of the two causal-leg log Bayes factors; not a joint BF |
| factor_log_e_mediation, factor_e_q_EBH | Minimum causal-leg log e-value and proteome-wide e-BH adjusted value |
| factor_p_XM_balanced, factor_p_MY_balanced | Experimental residual-scaled leg p-values under mean-zero pleiotropy/InSIDE |
| factor_balanced_conjunction_p, factor_balanced_conjunction_q_BH, factor_balanced_conjunction_q_BY | Experimental 2-of-2 balanced-pleiotropy p-value, independence/PRDS BH adjustment, and arbitrary-dependence BY adjustment |
| factor_balanced_conjunction_q_AdaFilter | Experimental AdaFilter-BH adjustment; conditional on independent legs and weak cross-protein dependence |
| factor_log_e_XM_balanced, factor_log_e_MY_balanced, factor_log_e_mediation_balanced, factor_e_q_balanced_EBH | Balanced/InSIDE Student-t density-ratio e-values and arbitrary-dependence e-BH adjustment |
| factor_log_e_p2e_balanced_mediation, factor_e_q_p2e_balanced_EBH | Fixed p-to-e calibration of the balanced 2-of-2 p-value and arbitrary-dependence e-BH adjustment |
| factor_log_e_XM_adaptive, factor_log_e_MY_adaptive, factor_log_e_mediation_adaptive, factor_e_q_adaptive_EBH | Strict known-covariance Gaussian-mixture e-values and e-BH adjustment |
| factor_p_XM_strict, factor_p_MY_strict, factor_strict_conjunction_p, factor_strict_conjunction_q_BY | Gaussian score p-values after projecting out an allele-oriented pleiotropic intercept, with BY adjustment |
| factor_log_BF_directional_XM, factor_log_BF_directional_MY | Evidence for an allele-oriented pleiotropic intercept on each leg, averaged over slope presence |
| factor_log_BF_slope_only_*, factor_log_BF_directional_only_*, factor_log_BF_slope_directional_* | Component-model log BFs versus the neither-component model, retained for prior sensitivity |
| factor_PP_directional_XM, factor_PP_directional_MY | Fixed-prior posterior probabilities of the directional component |
| factor_directional_intercept_XM, factor_directional_intercept_MY (and `_se`) | Joint adjusted-profile estimates of the allele-oriented intercepts |
| factor_directional_collinearity_XM, factor_directional_collinearity_MY | Weighted design collinearity between the causal slope and oriented intercept; values near 1 indicate weak separation |
| factor_log_e_p2e_mediation, factor_e_q_p2e_EBH | Fixed p-to-e strict-model sensitivity and e-BH adjusted value |
| factor_PP_XM, factor_PP_MY, factor_PP_two_stage | Fixed-prior working-model posterior probabilities for each leg and both legs |
| factor_posterior_local_fdr, factor_posterior_cum_fdr, factor_posterior_rank | Bayesian expected-FDR diagnostics across the analyzed protein family |
| factor_nA, factor_nB, factor_ld_source, factor_cross_set_max_r2 | Instruments, LD source, and largest retained Set A/Set B cross-LD `r²` used by factor inference |
| factor_pattern | Descriptive nominal evidence pattern; not a biological state posterior |
| factor_mediation_status | Bayesian two-leg evidence plus instrument and regional identification gates |
| factor_frequentist_status | Conjunction/BY evidence plus the same identification gates |
| factor_ebh_status | Safe-e/e-BH evidence plus the same identification gates |
| factor_balanced_status, factor_balanced_bh_status, factor_fcr_bh_status, factor_fcr_by_status, factor_adafilter_status, factor_balanced_ebh_status, factor_balanced_p2e_status, factor_adaptive_ebh_status | Assumption-labeled decisions for the experimental analytical calibration and selected-confidence-set tracks |
| factor_posterior_status | Fixed-prior posterior expected-FDR selection plus identification gates |
| mediated_effect | β₁×β₂ |
| se_mediated | Delta-method SE for mediated effect |
| ELBO_M0 ... ELBO_M5 | Evidence lower bounds per scenario |
| converged | Whether CAVI converged |

### `.regional` — Signal-pair Results

One row per conditionally independent protein-outcome signal pair. It reports
the lead variants, H0-H4 posterior probabilities, H4/(H3+H4), lead-variant
`r2`, maximum cross-trait credible-set `r2`, and the shared/distinct/ambiguous
classification. The file contains only its header when no pair is testable.

### `.hyp` — Priors and Inference Settings

Scenario/effect priors and inference settings. These remain fixed by default;
with `--empirical-bayes`, the file reports the estimated values.

### Analytical sensitivity calibration

For label-free sensitivity analysis, `scripts/calibrate_bmediator_analytic.py`
computes a composite-null mediation p-value using a prespecified global or
protein-specific bound on LD/pleiotropic bias, followed by BY or BH adjustment.
See [docs/ANALYTIC_CALIBRATION.md](docs/ANALYTIC_CALIBRATION.md) for assumptions,
usage, and limitations. This is distinct from simulation-trained empirical
calibration.

---

## Options Reference

### Required

| Flag | Description |
|------|-------------|
| `--rf-sumstat <file>` | RF GWAS summary statistics |
| `--pqtl-sumstat <file>` | pQTL summary statistics (multi-protein legacy mode) |
| `--protein-gwas-list <file>` | Per-protein GWAS manifest; use instead of `--pqtl-sumstat` |
| `--cancer-sumstat <file>` | Cancer GWAS summary statistics |
| `--protein-info <file>` | Protein gene annotation |
| `--out <prefix>` | Output file prefix |

### Instrument Selection

| Flag | Default | Description |
|------|---------|-------------|
| `--p-thresh-rf <val>` | 5e-6 | p-value threshold for RF instruments |
| `--p-thresh-cis <val>` | 5e-6 | p-value threshold for cis-pQTL instruments |
| `--cis-window <kb>` | 1000 | Cis window in kb (±1Mb default) |
| `--clump-kb <kb>` | 10000 | LD clumping window in kb |
| `--clump-r2 <val>` | 0.1 | LD r-squared threshold for RF and cis instrument clumping |

### Prior Specification

| Flag | Default | Description |
|------|---------|-------------|
| `--prior-p0 <val>` | 0.85 | Prior prob of null |
| `--prior-p1 <val>` | 0.03 | Prior prob of mediation |
| `--prior-p2 <val>` | 0.05 | Prior prob of RF→PP only |
| `--prior-p3 <val>` | 0.03 | Prior prob of RF→disease direct only |
| `--prior-p4 <val>` | 0.02 | Prior prob of PP→disease only |
| `--prior-p5 <val>` | 0.02 | Prior prob of correlated/shared pleiotropy |
| `--sigma2-beta1 <val>` | 0.1 | Prior variance for β₁ |
| `--sigma2-beta2 <val>` | 0.1 | Prior variance for β₂ |
| `--sigma2-beta3 <val>` | 0.1 | Prior variance for β₃ |

### Inference Control

| Flag | Default | Description |
|------|---------|-------------|
| `--structural-method <mode>` | legacy-six-state | Use `factorized` for the experimental analytical engine |
| `--sampling-corr-rf-pqtl <val>` | 0 | RF/pQTL estimation-error correlation |
| `--sampling-corr-rf-outcome <val>` | 0 | RF/outcome estimation-error correlation |
| `--sampling-corr-pqtl-outcome <val>` | 0 | pQTL/outcome estimation-error correlation |
| `--factor-min-set-a <int>` | 3 | Minimum exact Set A associations for confirmatory factor inference |
| `--factor-min-set-b <int>` | 3 | Minimum independent Set B instruments |
| `--factor-alpha <val>` | 0.05 | Leg-test and BY decision threshold |
| `--factor-bf-threshold <val>` | 10 | Minimum BF required for each causal leg |
| `--factor-prior-xm <val>` | 0.50 | Fixed prior probability for an RF-to-protein effect |
| `--factor-prior-my <val>` | 0.25 | Fixed prior probability for a protein-to-outcome effect |
| `--factor-prior-directional <val>` | 0.10 | Fixed prior probability for allele-oriented pleiotropy on either leg |
| `--factor-directional-variance <val>` | 0.01 | Fixed prior variance for the oriented pleiotropic intercept |
| `--factor-pleio-sd-xm <val>` | 0.1 | Half-normal prior SD for RF-to-protein residual heterogeneity |
| `--factor-pleio-sd-my <val>` | 0.1 | Half-normal prior SD for protein-to-outcome residual heterogeneity |
| `--factor-pleio-sd-xy <val>` | 0.1 | Half-normal prior SD for residual RF-to-outcome heterogeneity |
| `--max-cavi-iter <int>` | 200 | Max CAVI iterations per scenario |
| `--elbo-tol <val>` | 1e-6 | ELBO convergence tolerance |
| `--max-eb-iter <int>` | 20 | Max empirical Bayes iterations when enabled |
| `--eb-tol <val>` | 1e-4 | EB convergence tolerance |
| `--empirical-bayes` | off | Estimate priors/scales from analyzed proteins; exploratory |
| `--legacy-adaptive-priors` | off | Reuse IVW evidence in local priors; compatibility only |
| `--regional-prior-pp <val>` | 1e-4 | Per-variant protein association prior |
| `--regional-prior-outcome <val>` | 1e-4 | Per-variant outcome association prior |
| `--regional-prior-shared <val>` | 1e-8 | Shared prior; default is the independence product of trait priors |
| `--regional-min-both <val>` | 0.80 | Minimum P(shared or distinct) before regional classification |
| `--regional-min-shared <val>` | 0.80 | Minimum P(shared given both associated) for shared classification |
| `--regional-method <mode>` | ld-multisignal | LD-aware conditional multi-signal model; use `single` for the previous one-causal calculation |
| `--regional-max-signals <int>` | 10 | Maximum conditionally independent signals per trait |
| `--regional-signal-p <val>` | 5e-6 | Conditional p-value required to retain a signal |
| `--regional-coverage <val>` | 0.95 | Credible-set posterior coverage |
| `--regional-high-ld-r2 <val>` | 0.80 | Threshold used to label distinct signals as high LD |
| `--allow-unresolved-selection` | off | Permit selection without regional resolution |
| `--threads <int>` | 1 | Number of threads |
| `--verbose` | off | Detailed logging |

### Direction Consistency

These options use the RF→disease IVW direction as a benchmark for the mediated
effect direction, sign(β₁×β₂). The default is conservative: report direction
diagnostics but do not change posterior scenario probabilities.

| Flag | Default | Description |
|------|---------|-------------|
| `--direction-mode report` | report | Report direction columns only; rank/select by P(M1) |
| `--direction-mode prioritize` | report | Keep raw posteriors, but rank/select by P(M1)×Pr(direction consistent) |
| `--direction-mode soft` | report | Add a soft M1 penalty: `direction_weight × log(Pr(direction consistent))` before posterior normalization |
| `--direction-mode hard` | report | Remove M1 as an eligible scenario when Pr(direction consistent) is below `direction_min_prob` |
| `--direction-weight <val>` | 1.0 | Strength of the soft-mode direction penalty |
| `--direction-min-prob <val>` | 0.80 | Minimum consistency probability for hard mode |

---

## Example Data

Download the [BMEDIATOR example data from NIH Box](https://nih.app.box.com/folder/414412700653?s=kfsiv6jl3y166fwfbwhoznkodt504jec).
The example uses GRCh38 BMI and CHD summary statistics, UKB-PPP COMBINED protein
summary statistics, and the European GRCh38 `G1000plink` LD reference. Run the
following commands from the top-level example-data directory.

### One-protein example

This command analyzes IL6R as the only candidate mediator:

```bash
/path/to/BMEDIATOR/bmediator \
    --rf-sumstat rf/bmi_giant_locke_eur_grch38.bmediator.tsv \
    --protein-gwas-list manifests/ukb_ppp_combined_il6r.protein_gwas_manifest.txt \
    --cancer-sumstat outcome/chd_finngen_r12_grch38.bmediator.tsv \
    --protein-info protein_panels/ukb_ppp_combined/protein_info.tsv \
    --bfile ld_reference/G1000plink \
    --structural-method factorized \
    --out bmi_chd_il6r
```

The one-protein manifest contains:

```text
OID20385    protein_panels/ukb_ppp_combined/proteins/IL6R__OID20385.tsv
```

### Protein-manifest example

This command analyzes all five UKB-PPP proteins in the supplied manifest:

```bash
/path/to/BMEDIATOR/bmediator \
    --rf-sumstat rf/bmi_giant_locke_eur_grch38.bmediator.tsv \
    --protein-gwas-list manifests/ukb_ppp_combined.protein_gwas_manifest.txt \
    --cancer-sumstat outcome/chd_finngen_r12_grch38.bmediator.tsv \
    --protein-info protein_panels/ukb_ppp_combined/protein_info.tsv \
    --bfile ld_reference/G1000plink \
    --structural-method factorized \
    --out bmi_chd_ukb_ppp_all5
```

To run the supplied five-protein deCODE panel, change the manifest and protein
annotation paths:

```bash
/path/to/BMEDIATOR/bmediator \
    --rf-sumstat rf/bmi_giant_locke_eur_grch38.bmediator.tsv \
    --protein-gwas-list manifests/decode.protein_gwas_manifest.txt \
    --cancer-sumstat outcome/chd_finngen_r12_grch38.bmediator.tsv \
    --protein-info protein_panels/decode/protein_info.tsv \
    --bfile ld_reference/G1000plink \
    --structural-method factorized \
    --out bmi_chd_decode_all5
```

### Expected results

Each command creates four files:

```text
<output prefix>.mediation
<output prefix>.regional
<output prefix>.hyp
<output prefix>.instruments
```

For the one-protein command, `bmi_chd_il6r.mediation` contains one protein row:

| Protein | Gene | factor_nA | factor_nB | log BF XM | log BF MY | PP two-stage | Regional shared | Regional distinct | Mediation status |
|---------|------|-----------|-----------|-----------|-----------|--------------|-----------------|-------------------|------------------|
| OID20385 | IL6R | 75 | 87 | -2.508583 | 11.081333 | 0.075255 | 0.064727 | 0.735029 | NO_TWO_STAGE_BAYES_EVIDENCE |

IL6R has strong protein-to-CHD evidence in this analysis, but no corresponding
BMI-to-IL6R leg (`factor_p_XM=0.535829`). It is therefore not supported as a BMI
mediator. Its regional configuration is also ambiguous. Factorized effect
estimates use the balanced/InSIDE adjusted-profile model; the oriented
directional-intercept fit is reported separately as a sensitivity diagnostic.

For the five-protein command, `bmi_chd_ukb_ppp_all5.mediation` contains:

| Protein | Gene | log BF XM | log BF MY | PP two-stage | Shared signals | Regional interpretation | Final posterior status |
|---------|------|-----------|-----------|--------------|----------------|-------------------------|------------------------|
| OID20235 | PCSK9 | 4.508972 | 22.203667 | 0.989110 | 1 | UNRESOLVED_SINGLE_SHARED_SIGNAL | UNRESOLVED_SINGLE_SHARED_SIGNAL |
| OID20049 | NPPB | 2.841584 | -0.174359 | 0.206692 | 0 | LD_CONFIGURATION_AMBIGUOUS | NOT_SELECTED_BY_POSTERIOR_FDR |
| OID20385 | IL6R | -2.508583 | 11.081333 | 0.075255 | 0 | LD_CONFIGURATION_AMBIGUOUS | NOT_SELECTED_BY_POSTERIOR_FDR |
| OID20407 | ANGPTL3 | -0.471483 | -1.527416 | 0.025931 | 0 | UNRESOLVED_NO_OUTCOME_SIGNAL | NOT_SELECTED_BY_POSTERIOR_FDR |
| OID20650 | VEGFA | 2.295403 | -3.593386 | 0.008254 | 0 | UNRESOLVED_NO_OUTCOME_SIGNAL | NOT_SELECTED_BY_POSTERIOR_FDR |

PCSK9 passes the two-leg working-model posterior threshold, but the regional
model recovers only one shared signal. A single shared variant cannot separate
protein mediation from an exactly aligned direct effect, so the result is
`UNRESOLVED_SINGLE_SHARED_SIGNAL`, not a confirmed mediator. None of the five
proteins is selected by the analytically calibrated e-BH procedure. These
values were generated on Biowulf with the GRCh38 European reference supplied
in the example-data folder and the same-sample instrument-selection setting.

The deCODE manifest run gives the same qualitative conclusion:

| Protein | Gene | log BF XM | log BF MY | PP two-stage | Shared signals | Final status |
|---------|------|-----------|-----------|--------------|----------------|--------------|
| SeqId_5231_79 | PCSK9 | 2.391288 | 21.509874 | 0.916161 | 1 | UNRESOLVED_SINGLE_SHARED_SIGNAL |
| SeqId_15602_43 | IL6R | -1.230570 | 8.791513 | 0.225979 | 0 | NO_TWO_STAGE_BAYES_EVIDENCE |
| SeqId_16751_15 | NPPB | -0.393690 | -0.774358 | 0.053656 | 0 | NO_TWO_STAGE_BAYES_EVIDENCE |
| SeqId_10391_1 | ANGPTL3 | 2.888711 | -1.786530 | 0.050105 | 0 | NO_TWO_STAGE_BAYES_EVIDENCE |
| SeqId_2597_8 | VEGFA | 5.858941 | -2.322421 | 0.031554 | 0 | NO_TWO_STAGE_BAYES_EVIDENCE |

The legacy engine assigned `P(M1)` near 1 to PCSK9, IL6R, and ANGPTL3 in one
or both panels. Those probabilities are not estimates of the factorized target
and should not be compared numerically to `factor_PP_two_stage`. Under the new
analysis, IL6R and ANGPTL3 lack evidence for both causal legs, while PCSK9 has
two-leg evidence but remains nonidentifiable at a single shared signal. No
example protein is a confirmed mediator, and none is selected by e-BH.

### Frozen post-selection effect validation

The post-selection confidence-set procedure was evaluated without post hoc
tuning at commit `bc5d8f1`. The held-out design contained 2,700 simulated
families of 500 proteins (1.35 million analyses). All 87 prespecified component
criteria passed. BH mean false coverage proportion was 0.01157 with strong
instruments, 0.00165 with moderate instruments, 0 with four instruments, and
0.02801 under balanced heterogeneity. BY mean false coverage proportion was
0.00027 with strong instruments, 0.00114 under balanced heterogeneity, and
0.00033 under dense cross-protein dependence. Every selected effect family was
complete and none was unbounded. Declared sample overlap failed closed for all
39,923 tested proteins without emitting a finite selected confidence set.

This is an assumption-conditional, finite-grid validation. BH mean false
coverage proportion rose to 0.29491 under directional pleiotropy, which
violates the balanced/InSIDE model. The near-threshold weak-instrument cell
made no selections and therefore does not validate weak-instrument coverage.
The underlying Student-t profile intervals are not weak-instrument-exact.
Full provenance and compact summaries are in
`sim/results/factorized_effect_validation_20260906_v1/`.

---

## Interpreting Results

For factorized runs, start with `factor_two_stage_status` and the
multiple-testing status matching the inferential framework you prespecified:
`factor_ebh_status` for analytical e-BH, `factor_frequentist_status` for the
residual-scaled BY analysis, or `factor_posterior_status` for the fixed-prior
working Bayesian analysis. A statistical selection is not an identified
mediation result unless the regional and instrument checks also return a
supported conditional status. `UNRESOLVED_*` values are results, not missing
labels.

The following columns apply to the legacy six-state mode and are `nan` in a
factorized run:

- **P(M1) > 0.8**: Strong structural evidence for M1, not necessarily identified mediation
- **P_mediator_ld_resolved > 0.8**: Strong conditional evidence after regional and instrument checks
- **P(M2) dominant**: RF affects protein, but protein does not causally affect cancer
- **P(M3) dominant**: RF appears to affect disease directly/residually, without this protein mediating
- **P(M4) dominant**: Protein appears disease-relevant, but not an RF mediator
- **P(M5) dominant**: Apparent mediation signal driven by correlated/shared pleiotropy
- **P(M0) dominant**: No evidence for protein involvement

In factorized mode, use `factor_beta1`, `factor_beta2`, and `factor_indirect`
only when `factor_effect_estimator` is not an unresolved value. In legacy mode,
`beta2`, `se_beta2`, and `mediated_effect` retain their original meanings.
