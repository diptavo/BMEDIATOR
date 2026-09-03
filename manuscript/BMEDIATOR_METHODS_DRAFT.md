# BMEDIATOR: factorized Bayesian inference for two-stage molecular mediation using genetic instruments

## Introduction

Molecular traits may transmit part of the effect of an epidemiologic risk
factor to disease, but coincident genetic associations do not by themselves
establish mediation. A variant associated with a risk factor, protein, and
disease can arise from a mediated pathway, a residual risk-factor effect,
horizontal pleiotropy, or distinct causal variants correlated through linkage
disequilibrium (LD). These mechanisms can generate similar marginal summary
statistics, particularly at loci represented by one association signal.

Two-stage Mendelian randomization (MR) provides a natural mediation design:
genetic instruments for the risk factor estimate its effect on a molecular
trait, while separate molecular-trait instruments estimate the effect of that
trait on disease. Existing workflows commonly analyze these legs separately
and combine significance decisions informally. Conversely, a single
multi-state model can force mediation and pleiotropy into mutually exclusive
classes even though they may coexist. Such classification is not identified
when a pleiotropic alternative is allowed to reproduce the same covariance as
the mediated path.

BMEDIATOR addresses this problem through instrument partitioning and
factorized inference. Risk-factor instruments outside the molecular cis region
and molecular cis instruments not associated with the risk factor have
different identifying roles. The two causal legs, residual risk-factor effect,
and correlated pleiotropy are estimated as separate components rather than as
exclusive biological states. The two-stage null is tested by an
intersection-union procedure, while fixed-prior Bayes factors quantify
evidence for each leg. A separate LD-aware regional analysis distinguishes
shared from distinct molecular and disease association signals. The resulting
mediation claim is explicitly conditional on instrument validity and the
absence of direct effects exactly aligned with molecular instrument strength.

## Overview of methods

Let `X` denote a risk factor, `M_j` molecular trait `j`, and `Y` an outcome.
The pathway coefficients are `beta_1j` for `X -> M_j`, `beta_2j` for
`M_j -> Y`, and `beta_3j` for the residual `X -> Y` path after accounting for
the molecular path. BMEDIATOR partitions variants separately for every
molecular trait:

- Set A contains genome-wide risk-factor instruments outside the molecular
  cis region. Exact molecular associations at these variants identify
  `beta_1j`.
- Set B contains cis molecular-QTL instruments that are not risk-factor
  instruments. These variants identify `beta_2j` independently of Set A.
- Set C contains variants satisfying both definitions. They inform regional
  configuration analysis but are excluded from confirmatory estimation of the
  two causal legs because their two roles cannot be separated without further
  assumptions.

Each causal leg is fit with an LD-aware errors-in-variables likelihood. A
zero-centered normal prior is integrated deterministically to obtain a Bayes
factor against a point null. In parallel, a conservative conditional score
test is calculated for each leg. The mediation null is

```text
H0_med,j: beta_1j = 0 OR beta_2j = 0,
```

so its valid intersection-union p-value is

```text
p_med,j = max(p_XM,j, p_MY,j).
```

Benjamini-Yekutieli correction is applied across molecular traits. Regional
H3/H4 probabilities are not treated as p-values; they form a separate
Bayesian gate for whether molecular and outcome signals are distinct or
shared. Confirmatory reporting requires evidence for both causal legs,
adequate independent instruments, a matched LD reference, independent GWAS
estimation errors, and a shared regional signal.

## Detailed methods

### Input data and harmonization

The analysis requires risk-factor GWAS summary statistics, outcome GWAS
summary statistics, full cis and risk-factor-variant association statistics
for each molecular trait, gene coordinates, and an ancestry- and
genome-build-matched PLINK reference panel. Effect estimates are aligned to the
same reference-panel allele. Variants with incompatible alleles, invalid
standard errors, disallowed allele frequencies, or absence from the reference
panel are removed.

Risk-factor instruments and cis molecular-QTL instruments are selected using
prespecified significance thresholds and LD clumping. Optional HEIDI and
Steiger filters are applied before factor inference. A variant belongs to only
one of Sets A, B, and C for a given molecular trait. Proxy-projected molecular
associations are excluded from the factorized first stage because their
reported molecular-QTL SE does not include uncertainty from estimating proxy
LD.

### Summary-statistic covariance

For one causal leg, let `x` and `y` be vectors of exposure and outcome
association estimates, with standard-error diagonal matrices `D_x` and `D_y`.
Let `R` be the signed LD correlation matrix from the reference panel. The
sampling covariance matrices are

```text
S_x = D_x R D_x,
S_y = D_y R D_y.
```

If the two GWAS have correlated estimation errors, a scalar error correlation
`c_xy` can be supplied and

```text
C_xy[i,k] = c_xy se_x[i] R[i,k] se_y[k].
```

The current confirmatory decision is withheld whenever a nonzero overlap
correlation is declared because selection-aware overlap inference is not yet
implemented. The covariance is nevertheless used for exploratory effect and
Bayes-factor calculations.

### Errors-in-variables likelihood

For candidate effect `b`, define

```text
r(b) = y - b x,
V(b) = S_y + b^2 S_x - b(C_xy + C_xy').
```

The limited-information likelihood is

```text
ell(b) = -1/2 {log|V(b)| + r(b)'V(b)^-1 r(b) + n log(2 pi)}.
```

BMEDIATOR maximizes this scalar likelihood using a deterministic grid followed
by golden-section refinement. The local observed curvature supplies the
reported standard error. Covariance calculations use Cholesky factorization
with a small escalating ridge only when reference-panel sampling noise makes a
matrix numerically non-positive-definite.

For a prespecified normal alternative `b ~ N(0,W)`, evidence against `b=0` is

```text
BF_10 = integral exp{ell(b)} N(b;0,W) db / exp{ell(0)}.
```

The one-dimensional integral is evaluated by deterministic Simpson
quadrature. `W` is fixed before analysis and written to the hyperparameter
file. The weaker of the two leg-specific log Bayes factors is reported as
`factor_min_log_BF`; this is a conservative evidence summary and is not
described as a Bayes factor for the composite union null.

### Risk-factor-to-molecular effect

The first leg applies the likelihood to exact molecular associations `alpha_A`
and risk-factor associations `gamma_A` in Set A:

```text
alpha_A = beta_1j gamma_A + residual_A.
```

Because Set A was selected in the risk-factor GWAS, the molecular association
test is conditional on the observed risk-factor associations and requires
independent risk-factor and molecular GWAS estimation errors for confirmatory
calibration.

### Molecular-to-outcome effect

The second leg uses Set B only:

```text
Gamma_B = beta_2j alpha_B + direct_B + sampling_error_B.
```

Using Set B prevents covariance among risk-factor instruments from generating
the evidence for `beta_2j`. Balanced direct effects are accommodated by the
overdispersion calculation in the score test. Causal interpretation still
requires an exclusion restriction: direct effects must not be systematically
proportional to `alpha_B`.

### Residual risk-factor effect and pleiotropy

After estimating `beta_2j`, BMEDIATOR constructs

```text
Gamma_A* = Gamma_A - beta_2j alpha_A
```

and estimates `beta_3j` from `Gamma_A*` and `gamma_A`, propagating molecular
and second-stage uncertainty into its marginal SE. Thus mediation and a
residual risk-factor effect can both be present.

Residual molecular and outcome effects at Set A are used to calculate a
de-noised correlation diagnostic. This quantity may indicate coexisting
pleiotropy, but it is not used in the confirmatory mediation p-value or treated
as a fully calibrated test.

### Conditional score tests

For each leg, the score under `b=0`, conditional on `x`, is

```text
U = x' S_y^-1 y,
I = x' S_y^-1 x,
z = U / sqrt(I phi).
```

The scale `phi` is the whitened residual sum of squares divided by `n-1`,
truncated below at one so heterogeneity can only increase uncertainty. A
Student t reference with `n-1` degrees of freedom is used. Under independent
GWAS errors and the proportional overdispersion working model, this gives a
small-sample test that is conservative when the reported sampling variance is
sufficient.

The intersection-union p-value `max(p_XM,j,p_MY,j)` rejects mediation only when
both required legs reject their respective nulls. Across traits, BMEDIATOR
computes Benjamini-Yekutieli q-values, which control FDR under arbitrary
dependence under valid marginal p-values. This calibration is analytical: no
simulated labels, empirical null fitted across proteins, or truth-dependent
thresholds enter the calculation.

### LD-aware regional hypotheses

Within each cis region, BMEDIATOR uses unpruned molecular and outcome summary
statistics and signed reference-panel LD to identify conditional signals.
Wakefield approximate Bayes factors and credible sets are calculated for every
retained molecular-outcome signal pair under:

- H0: neither trait has a regional association;
- H1: molecular association only;
- H2: outcome association only;
- H3: both traits are associated through distinct causal variants;
- H4: both traits are associated through a shared causal variant.

The shared prior defaults to the product of the two marginal association
priors. `regional_PP_shared`, `regional_PP_distinct`, and
`regional_shared_given_both` report H4, H3, and H4/(H3+H4), respectively. H3
support rejects a simple cis-mediated interpretation; H4 support is necessary
but not sufficient for mediation because a shared variant can act
pleiotropically.

### Decision rule and output

The default Bayesian status requires at least two exact Set A associations,
at least two independent Set B instruments, both leg-specific BFs at least 10,
an LD reference, zero declared sampling-error correlation for the causal legs,
and regional support for H4. The frequentist status applies the same
identification gates and requires the BY q-value to be no greater than the
prespecified alpha. Effect estimates, SEs, leg-specific BFs and p-values,
conjunction p and q-values, regional probabilities, instrument counts, and
explicit unresolved states are retained for every molecular trait.

The mediated effect is `beta_1j beta_2j`; its SE uses a first-order delta
method. A positive status is a conditional causal conclusion under the stated
instrument and regional assumptions, not assumption-free proof of biological
mediation.

### Identification boundary

Suppose direct Set B effects satisfy `direct_B = lambda alpha_B` for every
instrument. Then

```text
Gamma_B = beta_2j alpha_B + direct_B
        = (beta_2j + lambda) alpha_B.
```

No number of such instruments separates the mediated slope from the aligned
direct slope using these summary statistics alone. BMEDIATOR does not assign a
data-derived label at this boundary. Its mediation interpretation is
conditional on excluding exact systematic alignment, supported where possible
by heterogeneity, regional, and external replication analyses.

### Validation strategy

Validation must assess the factorized estimand rather than six-state argmax
classification. Required experiments include null calibration under M2, M4,
and pleiotropic alternatives; power under M1; exact aligned-pleiotropy and
no-Set-B boundaries; correlated instruments with matched and mismatched LD;
weak instruments; balanced and directional pleiotropy; sample overlap;
winner's curse; sparse missing molecular statistics; multiple regional
signals; and proteome-wide FDR. Comparisons should include two-step MR,
MR-Egger, weighted-median or robust MR, SMR/HEIDI, colocalization, and
multivariable MR where their estimands are applicable. Real-data evaluation
should include independent molecular resources, positive controls, negative
controls, and replication across outcome studies.
