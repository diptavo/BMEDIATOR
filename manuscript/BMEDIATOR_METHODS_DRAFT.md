# BMEDIATOR: joint LD-aware Bayesian mediation analysis using genome-wide summary statistics

## Introduction

Genome-wide association studies identify risk-factor and disease-associated
loci, while large proteomic studies identify genetic determinants of circulating
protein abundance. Integrating these resources can identify proteins that lie
on a causal pathway from a modifiable or biological risk factor to disease.
Such analyses could prioritize mechanisms and intervention targets without
requiring all phenotypes to have been measured in the same individuals.

Mendelian randomization (MR) provides a framework for causal estimation using
genetic instruments, but mediation analysis with summary statistics poses a
more demanding identification problem than a conventional two-sample MR
analysis. Evidence is needed for both the risk-factor-to-protein and
protein-to-outcome paths while allowing a residual risk-factor-to-outcome path.
The same data must also distinguish mediation from horizontal pleiotropy,
correlated pleiotropy, and the appearance of a shared association created by
distinct causal variants in linkage disequilibrium (LD). An analysis based on a
single regional variant cannot, without an exclusion assumption, distinguish
`variant -> protein -> outcome` from `variant -> protein` plus an independent
direct effect of the same variant on outcome.

Existing summary-data approaches address parts of this problem. Two-step and
multivariable MR estimate components of a mediated pathway; colocalization and
conditional fine-mapping distinguish shared from distinct regional signals;
and pleiotropy-robust MR methods test or model violations of the exclusion
restriction. Applying these procedures sequentially, however, does not produce
one coherent comparison of the principal causal structures, and a
colocalization result alone is not evidence that the shared variant acts on the
outcome through the measured protein.

We developed BMEDIATOR, a Bayesian summary-statistic framework that compares
six structural scenarios for each candidate protein. The method combines an
independent genome-wide instrument likelihood with a joint regional likelihood
for risk-factor, protein, and outcome effects. The regional model includes
signed LD, correlated summary-error noise from sample overlap, multiple
conditional signals, and an explicit correlated-pleiotropy scenario. Its core
identification strategy is overidentification: when multiple independent
genetic perturbations are available, one set of pathway coefficients must
explain their cross-trait effect pattern. BMEDIATOR separately reports
structural model support, shared-versus-distinct regional evidence, and whether
the available instruments satisfy the method's conditional identification
gate. This separation is essential because no summary-statistic method can turn
an observationally equivalent one-signal locus into assumption-free evidence
of biological mediation.

## Overview of Methods

BMEDIATOR requires GWAS summary statistics for a risk factor (`X`) and outcome
(`Y`), full cis-region association statistics for each candidate protein
(`M`), genomic annotation for each protein, and an ancestry- and build-matched
PLINK LD reference panel. Alleles are harmonized to the reference panel before
instrument selection or regional analysis.

For each protein, variants are initially classified as RF-only instruments
(Set A), cis protein instruments that are not RF instruments (Set B), or
instruments satisfying both definitions (Set C). Set A supplies independent
genome-wide information about the `X -> M` and `X -> Y` relations. Full
cis-region statistics supply a separate LD-aware block. Under the default
joint analysis, Set B and Set C observations are excluded from the independent
block because they are represented in the regional block.

Six scenarios are compared. M0 contains no pathway coefficient; M1 contains
the `X -> M` and `M -> Y` paths and permits a residual `X -> Y` path; M2 contains
only `X -> M`; M3 contains only `X -> Y`; M4 contains only `M -> Y`; and M5
contains `X -> M` and `X -> Y` together with correlated protein and outcome
direct genetic effects but no mediated `M -> Y` path. M5 is therefore a
specific alternative in which correlated pleiotropy can mimic mediation.

The independent block is fitted by coordinate-ascent variational inference
(CAVI) with spike-and-slab nuisance effects. In the cis region, conditional
association analysis identifies candidate RF, protein, and outcome signals.
Their union defines a set of regional components. Marginal effects are
transformed to approximate joint effects using the regularized signed-LD
matrix. A random-effects structural equation model then evaluates all three
traits across all components simultaneously. Latent component effects are
integrated analytically, pathway coefficients are optimized, and scenario
evidence is approximated by Laplace integration.

The two likelihood blocks share the pathway coefficients. BMEDIATOR combines
their approximate evidence with a prior-overlap correction, rather than adding
two analyses that each count the same coefficient prior. Scenario priors are
prespecified and fixed by default. Approximate posterior scenario weights are
obtained by normalizing the combined evidence across M0-M5.

As an orthogonal regional diagnostic, BMEDIATOR evaluates conditional
protein-outcome signal pairs under H0-H4: neither trait associated, protein
only, outcome only, both traits with distinct causal variants, or both traits
with a shared causal variant. A protein is eligible for an LD-resolved mediation
report only when the shared-signal criterion, RF-instrument criterion, cis
instrument criterion, and multi-component criterion are satisfied. The
reported mediation probability remains the M1 model weight and remains
conditional on the exclusion and model assumptions.

## Detailed Methods

### Notation and input data

Let `X`, `M`, and `Y` denote the risk factor, candidate protein, and outcome.
For variant `j`, let `b_Xj`, `b_Mj`, and `b_Yj` be harmonized marginal genetic
association estimates with standard errors `s_Xj`, `s_Mj`, and `s_Yj`.
All association estimates must refer to the same effect allele and genome
build. The LD reference panel provides signed genotype correlations, not only
`r^2`, because the direction of correlation is required to reconstruct joint
effects.

The protein input must contain unpruned association statistics spanning the
configured cis window. A file containing only genome-wide significant lead
pQTLs is insufficient for the default analysis because conditional signal
selection, LD transformation, and H3/H4 evaluation require the surrounding
regional variants. Sample sizes enter preprocessing and quality-control
calculations where applicable; uncertainty in the likelihood is represented by
the supplied standard errors.

### Allele harmonization and quality control

RF, protein, and outcome alleles are aligned to the PLINK reference alleles.
Variants with incompatible alleles, invalid effect estimates or standard
errors, low minor-allele frequency, or excessive allele-frequency discrepancy
are excluded. Palindromic variants may be removed. RF and cis instruments are
LD clumped using the reference panel and configured physical and `r^2`
thresholds. Optional HEIDI and Steiger filters are available as preprocessing
diagnostics; they are not substitutes for the joint regional likelihood.

### Instrument partition

For a protein with cis region `C`, BMEDIATOR constructs:

- Set A: RF-associated variants outside `C`;
- Set B: protein-associated variants inside `C` that are not selected RF
  instruments; and
- Set C: selected RF instruments inside `C` that are also protein associated.

In compatibility analyses, all three sets enter the CAVI likelihood. In the
default `joint-ld` analysis, Set A remains in the independent block, whereas
the complete cis data replace Set B and Set C in the regional block. This
partition prevents duplicate use of the same cis association estimates.

### Structural scenarios

The pathway coefficients are `beta1` for `X -> M`, `beta2` for `M -> Y`, and
`beta3` for the residual or direct `X -> Y` path. The six scenarios impose:

| Scenario | Free pathway coefficients | Interpretation |
|---|---|---|
| M0 | none | no modeled pathway |
| M1 | `beta1`, `beta2`, `beta3` | partial mediation |
| M2 | `beta1` | RF affects protein only |
| M3 | `beta3` | RF affects outcome directly only |
| M4 | `beta2` | protein affects outcome but is not RF responsive |
| M5 | `beta1`, `beta3` | correlated pleiotropy; `beta2=0` |

Zero coefficients are exact scenario constraints. Free coefficients have
independent zero-mean normal priors with prespecified variances. M1 permits
`beta3` and therefore represents partial, rather than necessarily complete,
mediation.

### Independent genome-wide likelihood

For an RF instrument `k` in Set A, let `gamma_k` be its association with `X`,
`alpha_k` its association with `M`, and `Gamma_k` its association with `Y`.
The working equations are

```text
alpha_k = beta1 gamma_k + delta_k + e_Mk

Gamma_k = (beta3 + beta1 beta2) gamma_k
          + beta2 delta_k + psi_k + e_Yk.
```

`delta_k` and `psi_k` represent direct nuisance effects on protein and outcome.
They receive spike-and-slab priors. Under M5 they have a bivariate inclusion
model with nonzero prior correlation; otherwise nuisance terms are available
under every scenario so that the models are compared with a common nuisance
structure. Measurement errors are normal with variances given by the squared
summary-statistic standard errors. CAVI alternates updates for pathway
coefficients, inclusion probabilities, and nuisance-effect distributions until
the ELBO change is below tolerance or the maximum iteration count is reached.

In compatibility mode, Set B contributes

```text
Gamma_l = beta2 alpha_l + phi_l + e_Yl,
```

and Set C contributes both the first-stage and joint outcome equations. These
cis equations are omitted from CAVI under `joint-ld` because the regional model
uses the same observations with their full LD covariance.

### Conditional regional components

For each trait, marginal z scores are conditioned iteratively on previously
selected lead variants. If `S` indexes selected leads, the conditional score
for candidate `j` is computed from the signed correlations `R_jS`, the inverse
`R_SS^-1`, and the selected marginal scores. Selection stops when no remaining
variant passes `--regional-signal-p` or when `--regional-max-signals` is
reached.

For each retained signal, BMEDIATOR conditions on the other leads and computes
Wakefield approximate Bayes factors over regional variants. Variants are ranked
by their normalized Bayes factors to form a credible set at the configured
coverage. The union of conditional RF, protein, and outcome leads forms `K`
joint regional components. Duplicate lead variants across traits occur once.
The joint components are selected within the three-trait intersection. The
protein-outcome H3/H4 diagnostic is selected separately within the larger
protein-outcome intersection and therefore does not depend on whether an RF
summary-statistic row is available for a diagnostic variant.

### LD regularization and joint-effect transformation

Let `R` be the `K x K` signed correlation matrix among component variants. To
stabilize inversion, BMEDIATOR uses

```text
R_lambda = (1 - lambda) R + lambda I,
```

where `lambda=0.05` by default. The condition number of `R_lambda` is reported.
For trait `t`, the vector of approximate joint component effects is

```text
theta_hat_t = R_lambda^-1 b_t.
```

If the covariance of marginal summary errors is approximated by
`D_t R_lambda D_t`, where `D_t=diag(s_t)`, the covariance after transformation
is

```text
V_tt = R_lambda^-1 D_t R_lambda D_t R_lambda^-1.
```

This preserves uncertainty induced by LD instead of treating transformed
component estimates as independent observations.

### Correlated sampling error from sample overlap

When two GWAS share participants, their summary errors may be correlated. Let
`eta_st` denote the user-supplied error correlation between traits `s` and `t`.
BMEDIATOR uses the cross-trait covariance block

```text
V_st = eta_st R_lambda^-1 D_s R_lambda D_t R_lambda^-1.
```

The RF-protein, RF-outcome, and protein-outcome correlations default to zero.
They are analysis-design parameters, not estimated from the analyzed locus.
BMEDIATOR verifies that the resulting `3 x 3` trait error-correlation matrix is
positive definite.

### Joint regional structural model

For component `j`, define latent direct effects `g_j`, `d_j`, and `h_j` on RF,
protein, and outcome. The true joint component effects obey

```text
theta_Xj = g_j,
theta_Mj = beta1 g_j + d_j,
theta_Yj = (beta3 + beta1 beta2) g_j + beta2 d_j + h_j.
```

Equivalently, with `u_j=(g_j,d_j,h_j)'`,

```text
theta_j = A(beta) u_j,

A(beta) = [ 1                         0      0
            beta1                     1      0
            beta3 + beta1 beta2       beta2  1 ].
```

The latent effects are independent across components and normally distributed
with diagonal variances `sigma_g^2`, `sigma_d^2`, and `sigma_h^2`. These are
set by `--regional-prior-var-rf`, `--regional-prior-var-pp`, and
`--regional-prior-var-outcome`. Under M5 only,

```text
Cov(d_j, h_j) = rho_dh sigma_d sigma_h,
```

where `rho_dh` is fixed by `--regional-pleiotropy-rho`. This alternative allows
direct protein and outcome component effects to track each other without a
nonzero `beta2`.

Conditional on pathway coefficients, integrating `u_j` is analytic. The
genetic covariance for component `j` is

```text
G_j(beta) = A(beta) Sigma_u A(beta)',
```

and the full observed vector
`(theta_hat_X', theta_hat_M', theta_hat_Y')'` has a zero-mean multivariate
normal distribution with covariance equal to the transformed sampling
covariance plus the blockwise genetic covariance.

### Regional evidence and coefficient uncertainty

For each scenario, inactive coefficients are set to zero and active
coefficients are assigned their normal priors. BMEDIATOR maximizes the joint
log posterior from multiple starting points using coordinate pattern search.
A finite-difference Hessian at the optimum estimates local posterior precision.
If the Hessian is positive definite, the regional log evidence is approximated
by

```text
log Z_Rm approximately log p(data_R, beta_hat_m | M_m)
                       + q_m/2 log(2 pi)
                       - 1/2 log |H_m|,
```

where `q_m` is the number of active coefficients and `H_m` is the negative
Hessian. Reported regional log Bayes factors subtract the M0 evidence. Failed
matrix factorization or nonfinite evidence produces an explicit joint-model
failure state.

### Combining regional and genome-wide evidence

Let `Z_Gm` denote the independent-block evidence and `Z_Rm` the regional
evidence. Both integrations use the same prior `p(beta|M_m)`. Assuming the data
blocks are conditionally independent given `beta`, their joint evidence can be
written as

```text
Z_m = Z_Gm Z_Rm integral q_Gm(beta) q_Rm(beta) / p(beta|M_m) d beta,
```

where `q_Gm` and `q_Rm` are the block-specific posterior approximations. The
normal-posterior overlap integral corrects the otherwise duplicated prior and
couples coefficient estimates across blocks. BMEDIATOR adds this term to the
independent-block ELBO and the regional log Bayes factor. Constants common to
all scenarios cancel during normalization.

For a coefficient with block-specific posterior means `mu_G`, `mu_R`,
variances `v_G`, `v_R`, and prior variance `v_0`, the combined normal precision
and mean are

```text
p_C  = 1/v_G + 1/v_R - 1/v_0,
v_C  = 1/p_C,
mu_C = v_C (mu_G/v_G + mu_R/v_R).
```

The M1 mediated-effect estimate is `beta1 beta2`. Its standard error uses a
first-order delta approximation and the marginal normal variances of `beta1`
and `beta2`.

### Posterior scenario weights

Let `E_m` be the combined approximate log evidence and `pi_m` the prespecified
prior probability of scenario `m`. The reported model weight is

```text
P_Mm = pi_m exp(E_m) / sum_r pi_r exp(E_r).
```

Scenario and effect priors remain fixed by default. Empirical-Bayes estimation
is available as an exploratory option but changes the interpretation because
the analyzed protein collection then contributes to its own priors.

### H0-H4 regional diagnostic

Every conditional protein signal is paired with every conditional outcome
signal. Signal-specific Wakefield Bayes factors are combined under H0-H4:

- H0: neither trait is regionally associated;
- H1: protein only;
- H2: outcome only;
- H3: both associated through distinct causal variants; and
- H4: both associated through the same causal variant.

The default shared prior is the product of the protein and outcome association
priors. Thus the prior does not enrich H4 before observing the data. Posterior
H3 and H4 probabilities, H4/(H3+H4), lead-pair `r^2`, and maximum
credible-set-pair `r^2` are reported for each pair.

### Conditional identification and selection

The joint likelihood requires at least two distinct regional components before
its evidence contributes to structural scenario ranking. This is a necessary
overidentification condition, not proof that every instrument is valid. A
single component is labeled `UNRESOLVED_SINGLE_COMPONENT` and cannot produce a
nonzero `P_mediator_ld_resolved`.

The confirmatory mediation score is nonzero only when: (i) a protein-outcome
signal pair passes the shared H4 criterion; (ii) at least two observed
RF-to-protein instruments are available; (iii) at least one cis-only protein
instrument is available; (iv) the LD reference and full regional data are
available; and (v) the regional configuration is not numerically unresolved.
When these conditions pass, `P_mediator_ld_resolved=P_M1`; otherwise it is
zero. This gate prevents a high structural M1 weight from being described as
LD-resolved mediation when the required design information is absent.

Proteins are ranked by the configured selection probability. By default this
is `P_mediator_ld_resolved`. Local false discovery rates are `1-p`, and ordered
cumulative means provide Bayesian FDR summaries. These quantities inherit the
calibration of the approximate model weights and must be validated for the
target analysis design.

### Direction and effect reporting

BMEDIATOR reports the sign of `beta1 beta2`, the sign of the RF-to-outcome IVW
estimate, and a normal-approximation probability that the signs agree. The
default `report` mode does not alter model evidence. Optional modes can use
direction consistency for ranking or penalization, but direction agreement is
not an identification test and cannot distinguish mediation from aligned
pleiotropy.

### Assumptions and limitations

Interpretation requires ancestry-matched LD, harmonized alleles and genome
builds, adequately strong instruments, recovery of the material regional
signals, appropriate effect and nuisance priors, and a plausible summary-error
correlation specification. The exclusion restriction must hold apart from the
forms of pleiotropy represented in the model.

Multiple components improve identification only when they supply genuinely
independent perturbations whose cross-trait pattern is informative. A
same-variant direct effect, or a pleiotropic process proportional across every
component, can be observationally equivalent to mediation. M5 represents one
correlated-pleiotropy structure but cannot exhaust all alternatives. Therefore,
an LD-resolved BMEDIATOR result is conditional statistical evidence and should
be triangulated with external cohorts, alternative proteomic platforms,
functional perturbation, tissue-specific evidence, and sensitivity analyses.

The current implementation uses conditional lead selection rather than a full
multi-variant sparse regression posterior, normal random-effects priors for the
joint regional model, fixed component-effect variances, and Laplace and
variational approximations. Its posterior weights require broad simulation and
empirical validation before they can be interpreted as calibrated frequencies.
