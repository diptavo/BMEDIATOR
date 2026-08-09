# BMEDIATOR: Bayesian Mediation MR

**Bayesian framework for identifying mediating plasma proteins between risk factors and disease outcomes using GWAS summary statistics.**

Version 1.1.0-dev

---

## Release Status

BMEDIATOR is distributed as a C++17 command-line research tool with helper R,
Python, and shell scripts. It is not an R package. The repository is intended to
contain source code, documentation, workflow scripts, and small synthetic test
data only; full GWAS, pQTL, outcome, PLINK reference, and generated analysis
files should be stored outside git.

---

## Overview

BMEDIATOR identifies which plasma proteins (PP) are genuine causal mediators on the pathway from a risk factor (RF) to a disease outcome (e.g., cancer), using only summary statistics from three sources:

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

For each protein, BMEDIATOR computes the posterior probability of six scenarios.
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

Set C instruments are the most informative for mediation and are handled with a dedicated likelihood that avoids double-counting. β₂ is identified by triangulating across Set B (clean cis instruments) and Set C (powerful but contaminated overlap instruments).

### Inference

- **Inner loop**: Coordinate Ascent Variational Inference (CAVI) with spike-and-slab priors on pleiotropy terms
- **Outer loop**: Empirical Bayes updates scenario priors (p₀,p₁,p₂,p₃,p₄,p₅) across all proteins

---

## Installation

```bash
git clone <repository-url>
cd BMEDIATOR
make
```

Requirements: a C++17 compiler and `make`.

OpenMP is optional. It is disabled by default for portability:

```bash
make clean
make USE_OPENMP=1
```

Use `USE_OPENMP=1` only with compilers that support `-fopenmp`.

---

## Test

Run the bundled smoke test:

```bash
make test
```

The test uses the small synthetic files in `testdata/` and writes temporary
outputs under `build/test/`.

---

## Quick Start

```bash
./bmediator \
    --rf-sumstat      bmi_gwas.txt \
    --pqtl-sumstat    deCODE_pqtl.txt \
    --cancer-sumstat  breast_cancer_gwas.txt \
    --protein-info    protein_annotations.txt \
    --out             bmi_brca_mediation
```

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

### Protein Annotation (`--protein-info`)

```
PROTEIN  GENE   CHR  START      END
IL6      IL6    7    22725889   22732002
TNF      TNF    6    31543344   31546112
```

---

## Output Files

### `.mediation` — Main Results

Tab-delimited, one row per protein, sorted by P(M1) descending:

| Column | Description |
|--------|-------------|
| Protein | Protein identifier |
| Gene | Gene name |
| nA, nB, nC | Number of instruments in each set |
| P_M0 ... P_M5 | Posterior scenario probabilities |
| P_mediator | P(M1), posterior support for partial mediation |
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
| mediated_effect | β₁×β₂ |
| se_mediated | Delta-method SE for mediated effect |
| ELBO_M0 ... ELBO_M5 | Evidence lower bounds per scenario |
| converged | Whether CAVI converged |

### `.hyp` — Empirical Bayes Hyperparameters

Estimated scenario priors and pleiotropy parameters.

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
| `--max-eb-iter <int>` | 20 | Max empirical Bayes iterations |
| `--eb-tol <val>` | 1e-4 | EB convergence tolerance |
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

## Example Workflow

```bash
# 1. Prepare input files (ensure instruments are pre-clumped at r²<0.1)

# 2. Run BMEDIATOR
bmediator \
    --rf-sumstat      bmi_instruments.txt \
    --pqtl-sumstat    ukb_ppp_pqtl.txt \
    --cancer-sumstat  breast_cancer_overall.txt \
    --protein-info    protein_gene_map.txt \
    --p-thresh-rf     5e-6 \
    --p-thresh-cis    5e-6 \
    --cis-window      1000 \
    --threads         8 \
    --out             bmi_brca

# 3. Identify mediators (P(M1) > 0.5 or Bayesian FDR < 0.05)
awk -F'\t' 'NR==1 {for (i=1; i<=NF; ++i) h[$i]=i; print; next} $(h["P_mediator"]) > 0.5' bmi_brca.mediation

# 4. Bayesian FDR control
# Select proteins with P(M1) > threshold t such that:
#   sum(1-P(M1)) / count(selected) < 0.05
```

---

## Interpreting Results

- **P(M1) > 0.8**: Strong evidence for partial mediation
- **P(M2) dominant**: RF affects protein, but protein does not causally affect cancer
- **P(M3) dominant**: RF appears to affect disease directly/residually, without this protein mediating
- **P(M4) dominant**: Protein appears disease-relevant, but not an RF mediator
- **P(M5) dominant**: Apparent mediation signal driven by correlated/shared pleiotropy
- **P(M0) dominant**: No evidence for protein involvement

Check `beta2` and `se_beta2` under M=1 for the estimated causal effect of the protein on cancer. The `mediated_effect` (β₁×β₂) quantifies how much of the RF→Cancer effect flows through this protein.

---

## Citation

If you use BMEDIATOR, please cite:

> Dey D. (2026) BMEDIATOR: A Bayesian framework for identifying mediating plasma proteins using summary-statistics-based Mendelian randomization. [manuscript in preparation]

Machine-readable citation metadata is available in `CITATION.cff`.

---

## License

GNU General Public License v3.0. See `LICENSE`.
