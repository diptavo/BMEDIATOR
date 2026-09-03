# BMEDIATOR: Bayesian Mediation MR

**Bayesian framework for identifying mediating plasma proteins between risk factors and disease outcomes using GWAS summary statistics.**

Version 1.2.0-dev 

---

## Overview

BMEDIATOR identifies which intermediate phenotypes (e.g. plasma proteins; PP) are true causal mediators on the pathway from a risk factor (RF) to a disease outcome (e.g., cancer), using only summary statistics from three sources:

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

For each protein, BMEDIATOR computes ELBO-based variational approximate model
probabilities for six scenarios.
Mediation is modeled as **partial mediation**: a protein may carry an identifiable
component of the RF→disease association while a residual RF→disease effect remains.

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

In the default `joint-ld` analysis, Set A enters the independent genome-wide
CAVI block. The complete cis region replaces the pruned Set B and Set C
observations in a joint RF-protein-outcome likelihood. This prevents duplicate
use of cis statistics and allows beta2 to be informed by multiple regional
components while accounting for their signed LD. The earlier Set A/B/C
likelihood remains available through compatibility modes.

### Inference

- **Inner loop**: Coordinate Ascent Variational Inference (CAVI) with spike-and-slab priors on pleiotropy terms
- **Default priors**: Prespecified scenario and effect priors remain fixed, avoiding reuse of the analyzed outcomes to construct their own model priors
- **Optional exploratory mode**: `--empirical-bayes` estimates hyperparameters across analyzed proteins
- **Joint regional likelihood**: Full mode jointly models RF, protein, and outcome effects at conditionally independent cis components using their signed LD and sampling covariance
- **Regional diagnostic**: The same run compares H0-H4 for every conditional protein-outcome signal pair and reports shared, distinct, or ambiguous configurations

`P_M1` is approximate structural support for mediation within the six-state model.
`P_mediator_ld_resolved` is nonzero only when full regional data support a
shared signal, at least two independent RF-to-protein instruments exist, and a
cis-only signal is available. This remains conditional on exclusion and
valid-instrument assumptions; see
[Identification and LD Resolution](docs/IDENTIFICATION.md) and
[Validation Status](docs/VALIDATION.md).

A manuscript-style technical draft with Introduction, Overview of Methods, and
Detailed Methods is available in
[manuscript/BMEDIATOR_METHODS_DRAFT.md](manuscript/BMEDIATOR_METHODS_DRAFT.md).

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
    --out manifest_result
```

Full mode requires the `.bed`, `.bim`, and `.fam` files represented by
`LD_reference_prefix`. All inputs must use the same genome build.
The default joint analysis also requires that RF, protein, and outcome summary
statistics overlap within each protein cis region; variants absent from any one
of the three inputs cannot enter the joint regional likelihood.

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

Tab-delimited, one row per protein, sorted by the configured
`selection_probability` descending (LD-resolved probability by default):

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
| regional_method | Regional inference method (`joint-ld` by default) |
| regional_rf_signals, regional_protein_signals, regional_outcome_signals | Number of conditionally independent signals detected for each trait |
| regional_signal_pairs | Number of protein-outcome signal pairs tested |
| regional_max_credible_set_pair_r2 | Maximum cross-trait credible-set LD for the selected signal pair |
| regional_joint_n_variants | Number of cis variants observed in all three GWAS and the LD reference |
| regional_joint_components | Number of unique conditional lead variants in the joint RF/protein/outcome model |
| regional_joint_condition_number | Condition number of the regularized component LD matrix |
| regional_joint_status | Whether the joint likelihood was evaluated, unresolved, or failed diagnostics |
| regional_joint_log_bf_M0 ... M5 | Regional log Bayes factor for each structural scenario relative to M0 |
| mediation_identifiability | Explicit LD-resolution and conditional-identification state |
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
| `--regional-method <mode>` | joint-ld | Joint LD-aware structural likelihood; `ld-multisignal` and `single` retain earlier compatibility behavior |
| `--regional-max-signals <int>` | 10 | Maximum conditionally independent signals per trait |
| `--regional-signal-p <val>` | 5e-6 | Conditional p-value required to retain a signal |
| `--regional-coverage <val>` | 0.95 | Credible-set posterior coverage |
| `--regional-high-ld-r2 <val>` | 0.80 | Threshold used to label distinct signals as high LD |
| `--regional-prior-var-rf <val>` | 0.04 | Prior variance for latent direct RF component effects |
| `--regional-prior-var-pp <val>` | 0.04 | Prior variance for latent direct protein component effects |
| `--regional-prior-var-outcome <val>` | 0.04 | Prior variance for latent direct outcome component effects |
| `--regional-ld-shrinkage <val>` | 0.05 | Shrinkage of component LD toward the identity matrix |
| `--overlap-rf-protein <val>` | 0 | RF-protein summary-error correlation |
| `--overlap-rf-outcome <val>` | 0 | RF-outcome summary-error correlation |
| `--overlap-protein-outcome <val>` | 0 | Protein-outcome summary-error correlation |
| `--regional-pleiotropy-rho <val>` | 0.5 | Prespecified direct-effect correlation under M5 |
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
    --out bmi_chd_ukb_ppp_all5
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

| Protein | Gene | nA | nB | nC | P_M1 | P_M5 | regional_PP_distinct | regional_PP_shared | Joint components | P_mediator_ld_resolved | mediation_identifiability | mediated_effect |
|---------|------|----|----|----|------|------|----------------------|--------------------|------------------|------------------------|---------------------------|-----------------|
| OID20385 | IL6R | 112 | 87 | 0 | 0.155580 | 0.333591 | 0.735029 | 0.064727 | 11 | 0.000000 | LD_CONFIGURATION_AMBIGUOUS | -0.000195 |

IL6R is not identified as an LD-resolved mediator. The regional diagnostic
finds 10 conditional protein signals and one outcome signal but does not resolve
H3 versus H4. Within the joint structural model, M3 has the largest weight
(0.510829), followed by M5 and M1.

For the five-protein command, `bmi_chd_ukb_ppp_all5.mediation` contains:

| Protein | Gene | P_M1 | P_M5 | regional_PP_distinct | regional_PP_shared | Protein signals | Outcome signals | Joint components | P_mediator_ld_resolved | mediation_identifiability | selected_fdr5 |
|---------|------|------|------|----------------------|--------------------|-----------------|-----------------|------------------|------------------------|---------------------------|---------------|
| OID20235 | PCSK9 | 0.016397 | 0.983476 | 0.000002 | 0.999998 | 10 | 10 | 16 | 0.016397 | LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL | NO |
| OID20385 | IL6R | 0.155580 | 0.333591 | 0.735029 | 0.064727 | 10 | 1 | 11 | 0.000000 | LD_CONFIGURATION_AMBIGUOUS | NO |
| OID20650 | VEGFA | 0.137487 | 0.856451 | 0.000000 | 0.000000 | 10 | 0 | 11 | 0.000000 | UNRESOLVED_NO_OUTCOME_SIGNAL | NO |
| OID20407 | ANGPTL3 | 0.004021 | 0.990973 | 0.000000 | 0.000000 | 10 | 0 | 19 | 0.000000 | UNRESOLVED_NO_OUTCOME_SIGNAL | NO |
| OID20049 | NPPB | 0.002667 | 0.996798 | 0.715503 | 0.000001 | 10 | 1 | 11 | 0.000000 | LD_CONFIGURATION_AMBIGUOUS | NO |

The H3/H4 diagnostic supports a shared PCSK9-CHD signal, but the joint
structural likelihood favors correlated pleiotropy (M5) over mediation (M1).
Therefore, no protein is selected at 5% FDR under the default analysis. These
values were generated on Biowulf with the production code and the GRCh38
European reference supplied in the example-data folder.

---

## Interpreting Results

- **P(M1) > 0.8**: Strong structural evidence for M1, not necessarily identified mediation
- **P_mediator_ld_resolved > 0.8**: Strong conditional evidence after regional and instrument checks
- **P(M2) dominant**: RF affects protein, but protein does not causally affect cancer
- **P(M3) dominant**: RF appears to affect disease directly/residually, without this protein mediating
- **P(M4) dominant**: Protein appears disease-relevant, but not an RF mediator
- **P(M5) dominant**: Apparent mediation signal driven by correlated/shared pleiotropy
- **P(M0) dominant**: No evidence for protein involvement

Check `beta2` and `se_beta2` under M=1 for the estimated causal effect of the protein on cancer. The `mediated_effect` (β₁×β₂) quantifies how much of the RF→Cancer effect flows through this protein.
