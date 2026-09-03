# BMEDIATOR: Bayesian mediation analysis with instrument partitioning and LD-aware regional resolution

## Introduction

Molecular traits can lie on the pathway from an epidemiologic risk factor to
disease, but genetic associations alone do not establish mediation. A variant
associated with a risk factor, protein, and disease may support a mediated
pathway, direct risk-factor effects, horizontal pleiotropy, or distinct causal
variants correlated through linkage disequilibrium (LD). These explanations
can produce similar marginal summary statistics.

BMEDIATOR is a summary-statistic framework that separates two inferential
questions. First, a six-state structural model compares mediation with null,
single-path, direct-effect, and correlated-pleiotropy alternatives. Second, an
LD-aware regional model asks whether protein and outcome associations are
consistent with shared or distinct causal variants. The structural and
regional results are reported separately because colocalization is necessary
for a simple cis-mediated interpretation but does not prove mediation.

The central design feature is instrument partitioning. Genome-wide
risk-factor instruments, protein-specific cis-pQTLs, and variants belonging to
both groups have different identifying roles. In particular, the distinction
between mediation and correlated pleiotropy depends on independent
protein-specific instruments and an exclusion restriction. No summary-statistic
method can distinguish a mediated effect from an unrestricted direct effect
that is exactly proportional to the protein effect at every instrument.

## Overview of methods

Let `X` denote the risk factor, `M` a protein, and `Y` the outcome. The pathway
coefficients are `beta1` for `X -> M`, `beta2` for `M -> Y`, and `beta3` for a
residual `X -> Y` path. BMEDIATOR compares six structural states:

| State | Active pathway | Interpretation |
|---|---|---|
| M0 | none | null |
| M1 | `beta1`, `beta2`, optionally `beta3` | partial mediation |
| M2 | `beta1` | risk factor affects protein only |
| M3 | `beta3` | residual/direct risk-factor effect only |
| M4 | `beta2` | protein affects outcome but is not risk-factor responsive |
| M5 | `beta1`, `beta3`, correlated RF-instrument residuals; `beta2=0` | correlated pleiotropy without mediation |

Variants are partitioned into Set A, genome-wide risk-factor instruments
outside the protein cis region; Set B, cis-pQTL instruments that are not
genome-wide risk-factor instruments; and Set C, variants satisfying both
definitions. Set A identifies the risk-factor pathways. Set B provides the
primary identifying evidence for a common protein-to-outcome coefficient.
Set C is modeled once in a joint equation and adds information, but cannot by
itself separate mediation from correlated pleiotropy.

For each protein, coordinate-ascent variational inference approximates the
evidence and coefficient posterior under each state. Normalized evidence
weights are reported as `P_M0` through `P_M5`. These are variational approximate
model weights, not exact posterior probabilities.

In full mode, unpruned cis-region summary statistics and a PLINK LD reference
are used to identify conditionally associated protein and outcome signals.
Each protein-outcome signal pair is evaluated under H0-H4: neither trait,
protein only, outcome only, both traits with distinct causal variants, or both
traits with one shared causal variant. This regional result gates confirmatory
mediation reporting but is not folded into the M0-M5 evidence.

## Detailed methods

### Summary-statistic inputs and harmonization

The required inputs are risk-factor GWAS summary statistics, outcome GWAS
summary statistics, one protein GWAS file per candidate protein, protein gene
coordinates, and an ancestry- and genome-build-matched PLINK reference panel.
Variants are aligned to a common effect allele. Ambiguous or incompatible
alleles, invalid standard errors, duplicates, and variants absent from the LD
reference are removed according to the documented quality-control rules.

Risk-factor instruments are selected at the configured significance threshold
and LD clumped. Cis-pQTL candidates are selected within the configured window
around each protein-coding gene and are subjected to instrument strength,
clumping, optional HEIDI-outlier, and optional Steiger filters. Set membership
is then assigned without duplicating an overlap variant across likelihood
terms.

### Structural equations

For risk-factor instrument `k` in Set A, let `gamma_k` be its effect on `X`,
`alpha_k` its effect on `M`, and `Gamma_k` its effect on `Y`. The model is

```text
alpha_k = beta1 * gamma_k + delta_k + e_Mk
Gamma_k = (beta3 + beta1 * beta2) * gamma_k
          + beta2 * delta_k + psi_k + e_Yk.
```

For protein-specific cis instrument `l` in Set B,

```text
Gamma_l = beta2 * alpha_l + phi_l + e_Yl.
```

For overlap instrument `c` in Set C,

```text
alpha_c = beta1 * gamma_c + delta_c + e_Mc
Gamma_c = beta3 * gamma_c + beta2 * alpha_c + psi_c + e_Yc.
```

The sampling errors are centered normal variables with variances given by the
reported standard errors. LD weights reduce duplicated information among
remaining correlated instruments.

The nuisance effects `delta`, `phi`, and `psi` receive spike-and-slab priors.
They allow sparse departures from the structural pathways in every state, so
model comparisons do not give one state a nuisance flexibility unavailable to
another. Pathway coefficients receive centered normal priors with fixed default
scales. Empirical-Bayes and data-adaptive scales are available only as
exploratory options.

### Definition and identification of M5

Under M5, `beta2=0`, while `beta1` and `beta3` are active. The residual
RF-to-protein effect `delta` and residual RF-to-outcome effect `psi` may have a
correlated bivariate spike-and-slab prior at RF-associated Set A/C instruments.
This represents a shared pleiotropic mechanism linked to risk-factor
instruments. It does not introduce correlation between `alpha_l` and `phi_l`
at protein-specific Set B instruments.

That restriction is the source of M1/M5 identifiability. Under M1, multiple
independent Set B instruments should show an outcome association proportional
to their protein effects through one common `beta2`, apart from independent
sparse direct effects. Under M5, there is no common Set B slope; correlated
residual evidence is confined to the RF-associated instrument pattern.

An unrestricted alternative is nonidentified. If `M=d` at a cis component,
M1 gives `Y=beta2*d+h` and therefore
`Cov(M,Y)=beta2*Var(d)`. A model with `beta2=0` but freely chosen
`Cov(d,h)` can reproduce the same covariance and variances. A comparison of
those two formulations is determined by priors and complexity penalties rather
than by distinct observable implications. Such a model is not used by
BMEDIATOR.

Consequently, a strong M1 interpretation requires multiple independent
RF-to-protein observations, at least one protein-specific cis instrument, and
the assumption that protein-specific direct effects are sparse and independent
of protein instrument strength. More than one independent Set B instrument is
preferable because proportionality can then be assessed across instruments.
Single-signal loci remain vulnerable to same-variant horizontal pleiotropy.

### Variational inference and structural evidence

For each state, BMEDIATOR alternates closed-form or conjugate updates for the
active pathway coefficients, nuisance inclusion probabilities, nuisance effect
moments, and residual quantities until the evidence lower bound (ELBO)
converges or the iteration limit is reached. State priors are prespecified.
The normalized structural weight for state `m` is

```text
P_Mm = exp(log_prior_m + ELBO_m) /
       sum_q exp(log_prior_q + ELBO_q).
```

Because an ELBO is an approximation to log evidence and approximation error
may differ by state, the resulting values require simulation and external
validation before being interpreted as calibrated probabilities. Directional
consistency and frequentist IVW summaries are reported as diagnostics.

### LD-aware regional H0-H4 analysis

The default `ld-multisignal` procedure extracts signed genotype correlations
from the PLINK panel. Protein and outcome signals are selected iteratively from
conditional summary statistics. For each retained lead, the method conditions
on other leads, calculates Wakefield approximate Bayes factors, and constructs
a posterior credible set.

Every protein-outcome signal pair is compared under:

- H0: neither trait has a regional association;
- H1: protein association only;
- H2: outcome association only;
- H3: both traits are associated through distinct causal variants; and
- H4: both traits are associated through the same causal variant.

The default shared prior equals the product of the two marginal association
priors. This represents prior independence of the two causal-status indicators
and avoids an automatic enrichment for H4. `regional_PP_shared` is the H4
posterior; `regional_PP_distinct` is the H3 posterior; and
`regional_shared_given_both` is H4 divided by H3+H4 when the denominator is
positive.

This analysis distinguishes one shared association from distinct variants in
LD. It does not distinguish mediation from horizontal pleiotropy through the
same variant. It is therefore used as a resolution gate rather than as a
structural M1/M5 likelihood.

### Decision outputs

`P_M1` reports structural support for M1. `P_mediator_ld_resolved` equals
`P_M1` only when the regional model supports a shared signal and the minimum RF
and cis instrument conditions are met; otherwise it is zero. The associated
`mediation_identifiability` field states whether the locus was shared,
distinct, ambiguous, missing an association, missing instruments, or missing
LD information.

The mediated effect is estimated as `beta1*beta2`, with uncertainty calculated
by a delta-method approximation. Selection and false-discovery summaries are
based on the configured structural or LD-resolved target. A positive call means
that the data favor mediation within the stated model and assumptions. It is
not assumption-free proof of a biological mechanism.

### Required sensitivity analyses

Primary analyses should report results across reasonable LD references,
conditional-signal thresholds, prior scales, and cis windows. Loci with nearly
collinear credible sets, one effective cis instrument, ancestry mismatch, or
substantial sample overlap should be labeled unresolved or sensitivity-limited.
Validation should include separate M1 and M5 data-generating mechanisms and
must demonstrate discrimination, not merely numerical convergence. Real-data
replication should use independent protein and outcome resources when
available.

