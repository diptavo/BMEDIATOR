# BMEDIATOR: factorized Bayesian inference for two-stage molecular mediation using genetic instruments

> **Historical draft, not current BMEDIATOR methods.** This document describes
> the factorized two-leg development method retained for audit compatibility.
> The current joint Bayesian development approach is JG-0.2.8, documented in
> `docs/JOINT_GRAPH_MODEL_V0_2.md`. JG-0.2.4 failed its frozen numerical
> completion gate, so neither this draft nor the joint model is ready for a
> submission claim. See `docs/PRODUCTION_READINESS.md`.

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
intersection-union procedure, while fixed-prior Bayes factors and posterior
expected-FDR values quantify model-based evidence for each leg. A separate LD-aware regional analysis distinguishes
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

Disjoint variant identifiers alone do not establish distinct instrument roles.
We therefore cross-clump Set A and Set B using chromosome-wide reference-panel
LD. Every risk-factor instrument within the molecular cis window is assigned
to Set C or excluded from confirmatory estimation, irrespective of its
molecular-association p-value. A cis molecular-QTL above the cross-set LD
threshold with any risk-factor instrument is also assigned to Set C and the
corresponding risk-factor instrument is removed from Set A. The maximum
retained cross-set `r²` is reported and threshold violations fail the causal
identification gate.

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
estimation errors, and at least two independently matched shared regional
signals. Even then, the interpretation remains conditional on exclusion and
independence assumptions.

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
prespecified significance thresholds, within-set clumping, and cross-set LD
exclusion. Optional HEIDI and
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

BMEDIATOR augments the covariance by `tau^2 I`, where `tau` is a residual
heterogeneity SD. Point estimation uses a generalized adjusted profile score:
the quadratic term is optimized in `b` without the residual-variance
determinant that otherwise attenuates weak-instrument estimates, while `tau`
is estimated from the residual likelihood. Independent instruments use the
adjusted-profile sandwich variance and correlated instruments use profile
curvature. Covariance
calculations use Cholesky factorization with a small escalating ridge only when
reference-panel sampling noise makes a matrix numerically
non-positive-definite.

For Bayesian evidence, define `s_i=sign(x_i)` and add an allele-oriented
pleiotropic component `eta s`, with `eta ~ N(0,W_eta)`. The method compares
four models: neither slope nor directional component, slope only, directional
only, and both. With `b ~ N(0,W)` and `tau ~ HN(s_tau)`, the directional
coefficient is integrated analytically and `b,tau` are integrated by
prespecified quadrature.

```text
BF_slope = p(data | slope present, averaged over eta presence)
           / p(data | slope absent, averaged over eta presence).
```

The integrals are represented by normalized deterministic Simpson-weighted
grids, giving proper finite discrete mixture distributions. `W` and `s_tau`
are fixed before analysis and written to the hyperparameter file. The weaker
of the two leg-specific log Bayes factors is reported as
`factor_min_log_BF`; this is an evidence summary and is not described as a
Bayes factor for the composite union null.

Numerical validity of the integrated evidence, balanced effect estimate, and
joint directional sensitivity fit are assessed separately. The reported
confirmatory effect is the balanced/InSIDE adjusted-profile slope. The joint
slope/oriented-intercept fit supplies a separate directional-pleiotropy
diagnostic; its collinearity does not make the balanced slope unavailable. If
balanced profile curvature is invalid, the software retains well-defined
projected tests and integrated evidence but withholds the ordinary effect
estimate and assigns `UNRESOLVED_EFFECT_ESTIMATION`.

### Fixed-prior Bayesian posterior FDR

Let `pi_XM` and `pi_MY` be prespecified probabilities of nonzero slopes and
let `pi_eta` be the prespecified directional-component probability. The leg
posterior probabilities are obtained by normalizing the four model weights;
equivalently they are slope-present probability marginalized over directional
presence:

```text
PP_l = P(slope present | data).
```

After cross-clumping confirmatory Sets A and B, the factorized working-model
probability of two nonzero slopes is
`PP_two-stage = PP_XM PP_MY`. The posterior local FDR is one minus this
quantity. Traits are sorted by local FDR and their running mean estimates the
posterior expected false-discovery proportion. The development defaults
`pi_XM=0.50`, `pi_MY=0.25`, and `pi_eta=0.10` were frozen before confirmatory
validation and are recorded in every analysis. This is Bayesian calibration under the fixed-prior
likelihood model, not a frequentist guarantee under model or prior
misspecification. The product is exact only under the working likelihood
factorization. Because nonzero LD can correlate sampling errors even for
different variants, the maximum retained cross-set `r²` is reported and a
configured-threshold violation fails closed.

### Risk-factor-to-molecular effect

The first leg applies the likelihood to exact molecular associations `alpha_A`
and risk-factor associations `gamma_A` in Set A:

```text
alpha_A = beta_1j gamma_A + residual_A.
```

Because Set A was selected in the risk-factor GWAS, the molecular association
test is conditional on the observed risk-factor associations and requires
independent risk-factor and molecular GWAS estimation errors for confirmatory
calibration. Exposure-only selection may use these same exposure estimates;
an independent discovery statistic is preferred for effect estimation and
Bayesian evidence to reduce winner's-curse bias.

### Molecular-to-outcome effect

The second leg uses Set B only:

```text
Gamma_B = beta_2j alpha_B + direct_B + sampling_error_B.
```

Using Set B prevents the same risk-factor instruments from generating the
evidence for `beta_2j`. The residual-scaled score permits a common Gaussian
variance multiplier relative to the declared LD covariance; it is not valid
for arbitrary balanced, sparse, or variant-specific direct effects. Causal
interpretation still requires an exclusion restriction: direct effects must
not be systematically proportional to `alpha_B`.

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

For each leg let `s_i=sign(x_i)`. After whitening by the reported outcome
covariance, we project both the exposure and outcome vectors away from `s`.
Writing the residualized vectors as `x_perp` and `y_perp`, the score under
`b=0`, conditional on `x`, is

```text
U = x_perp' y_perp,
I = x_perp' x_perp,
z = U / sqrt(I phi).
```

The scale `phi` is the whitened residual sum of squares divided by `n-2`. A
Student t reference with `n-2` degrees of freedom is used. Conditional on `x`,
this is exact under independent GWAS errors when the declared outcome
covariance is correct up to a common positive Gaussian scale multiplier.

A second score uses the reported outcome covariance without estimating
overdispersion and has a standard-normal null reference. It is exact under the
declared Gaussian covariance while allowing an arbitrary allele-oriented
intercept, but it is not robust to additional residual heterogeneity.

The intersection-union p-value `max(p_XM,j,p_MY,j)` rejects mediation only when
both required legs reject their respective nulls. Across traits, BMEDIATOR
computes Benjamini-Yekutieli q-values, which control FDR under arbitrary
dependence under valid marginal p-values. This calibration is analytical: no
simulated labels, empirical null fitted across proteins, or truth-dependent
thresholds enter the calculation.

### Safe e-values and arbitrary-dependence FDR

Let `T` denote the residual-scaled statistic above and `nu=n-2`. Conditional on
the selected exposure associations, `T` has density `f_nu`, the central Student
t density, under the leg null for every unrestricted common allele-oriented
intercept and positive common Gaussian scale multiplier. We prespecify

```text
g_nu(t) = (1/6) [
  sum_(mu in {2,4,6}) {f_nu(t-mu) + f_nu(t+mu)}/2
  + sum_(s in {2,4,8}) f_nu(t/s)/s
],
E_leg = g_nu(T)/f_nu(T).
```

Every shifted or scaled component and therefore `g_nu` is a proper density.
It follows directly that `E_0(E_leg|x)=1`. The symmetric shift and scale grids
have equal fixed weights, are recorded in the hyperparameter output, and are
not estimated from outcomes or simulation labels. The guarantee requires
independent exposure and outcome GWAS errors and correct LD/error covariance
up to the common scalar; it excludes arbitrary heterogeneous, sparse, and
exactly slope-aligned direct effects. If a nonzero sampling-error correlation
is declared for either causal leg, the software reports exploratory leg
statistics but does not calculate BY or e-BH q-values and labels the result
unresolved for sample overlap.

The composite mediation null is `H_XM,0 union H_MY,0`. The statistic

```text
E_med,j = min(E_XM,j, E_MY,j)
```

is valid under this union because it is no greater than the valid leg e-value
for whichever null is true. Let `E_(1) >= ... >= E_(m)` be the ordered
mediation e-values. Base e-BH selects the largest `k` satisfying
`E_(k) >= m/(alpha k)`. This controls FDR under arbitrary cross-protein
dependence when the leg likelihood and independence assumptions hold. The
implementation reports adjusted e-BH values in `factor_e_q_EBH` (Wang and
Ramdas, 2022).

This FDR statement applies to the two-leg statistical hypotheses. Regional
H3/H4 classification is a subsequent causal-identification assessment. The
post-gate subset is not claimed to inherit the e-BH guarantee without a
separate structured e-value argument.

### Balanced-pleiotropy score and adaptive partial-conjunction testing

The unrestricted allele-oriented intercept protects against directional
pleiotropy but can be nearly collinear with the causal slope when instrument
magnitudes have a narrow range. We therefore report a separate higher-power
working model that assumes mean-zero pleiotropy conditional on instrument
strength. In whitened coordinates this model is `y=beta*x+epsilon`, without an
intercept. The through-origin score is divided by the residual scalar
dispersion after fitting the slope and is compared with a Student t reference
with `n-1` degrees of freedom. The 2-of-2 mediation p-value remains the maximum
of the two leg p-values. BY adjustment supplies an arbitrary-dependence result
when the base p-values are valid under this balanced/InSIDE model.

We additionally apply ordinary BH directly to these protein-level
partial-conjunction p-values. The maximum is super-uniform under the mediation
union null whenever either null leg p-value is valid, regardless of dependence
between the two legs. BH therefore controls FDR under independence or PRDS
across proteins, while BY remains the reported arbitrary-dependence analysis.
The BH decision is kept separate from AdaFilter because it requires no
data-dependent filtering and weaker assumptions about dependence between the
two causal legs. Correlation alone does not establish PRDS, particularly for
two-sided Gaussian tests, so the ordinary-BH result is not labeled valid under
unspecified cross-protein dependence.

The two-leg hypothesis is also a partial-conjunction hypothesis. For protein
`j`, define filtering and selection p-values

```text
F_j = min(p_XM,j, p_MY,j),
S_j = max(p_XM,j, p_MY,j).
```

After sorting by `S_(1)<=...<=S_(m)`, let
`m_AF,(j)=sum_h 1{F_h<=S_(j)}`. The AdaFilter-BH adjusted value is

```text
q_AF,(j) = min_{h>=j} min{1, S_(h) m_AF,(h) / h}.
```

This procedure is more efficient than direct adjustment when few proteins
have even one credible leg. Its finite-sample FDR theorem requires independent
base p-values; the asymptotic extension permits weak dependence within each
study while retaining independence between studies (Wang et al., 2022).
Because proteomic errors and shared GWAS can induce stronger dependence,
BMEDIATOR labels this output assumption-conditional and retains BY and e-BH as
the more conservative dependence-robust tracks.

As a dependence-robust alternative to AdaFilter, we transform the balanced
partial-conjunction p-value with the fixed calibrator mixture

```text
e(p) = (1/4) sum_{kappa in {0.10,0.25,0.50,0.75}}
       kappa p^(kappa-1).
```

Each component integrates to one under a uniform p-value and has expectation
at most one for a super-uniform p-value. The maximum of the two leg p-values is
super-uniform under the union null without requiring independence between the
legs. Applying e-BH across proteins therefore permits arbitrary cross-protein
dependence, conditional on validity of the balanced/InSIDE leg tests. The grid
and equal weights are prespecified and are not estimated from simulations or
the analysis data.

We additionally construct a more information-preserving e-value directly from
each balanced Student t statistic. Under the balanced leg null, the statistic
has `n-1` degrees of freedom after estimating the through-origin slope and
common residual scale. We divide the prespecified equal-weight mixture of
symmetric shifted Student t densities (shifts 2, 4, and 6) and scaled Student t
densities (scales 2, 4, and 8) by the central null density. This density ratio
has expectation one under the leg null for every common positive scale. The
minimum of the two leg e-values is valid for the mediation union null and is
passed to e-BH. This procedure does not require independence between legs or
proteins, but it remains conditional on balanced/InSIDE pleiotropy, correct
LD-aware covariance up to a scalar, independent selection, and no causal-leg
sample overlap.

### Information-adaptive strict e-value

Let `Z` be the strict standard-normal score after projecting both whitened
vectors away from the allele-oriented intercept. For a fixed
signal-to-noise variance `g`, the likelihood ratio obtained by mixing the
alternative score mean under a zero-centered normal prior is

```text
E_g(Z) = (1+g)^(-1/2) exp{g Z^2/[2(1+g)]}.
```

Each component has expectation one under `Z~N(0,1)`. We average the frozen
grid `g={0.25,1,4,16}`, take the minimum across causal legs, and apply e-BH.
The grid is independent of outcome data and simulation labels. This strict
track retains the oriented-intercept nuisance model and arbitrary-dependence
e-BH guarantee, but requires the reported sampling covariance to be correct
rather than merely correct up to an unknown scalar.

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

The primary development Bayesian family rule uses the posterior cumulative
FDR. A separate evidence status retains the historical requirement that both
leg-specific BFs exceed 10. Confirmatory interpretation additionally requires
at least three exact Set A associations, at least three independent Set B
instruments, an LD reference, at least two independently matched H4 signal
pairs, and no strong residual heterogeneity on either causal leg. The robust
frequentist and safe-e statuses apply the same interpretation gates after BY
or e-BH selection. Strict Gaussian statuses are explicitly labeled as
exclusion-model sensitivities. Effect estimates, heterogeneity estimates,
SEs, leg-specific BFs and posterior probabilities, e-values and p-values,
family-level decisions, regional probabilities, instrument counts, and
explicit unresolved states are retained for every molecular trait.

The mediated effect is `beta_1j beta_2j`; its SE includes the second-order
product-of-variance term. Its confidence interval is the range of the four
corner products from Bonferroni simultaneous t intervals for the two legs,
using `n-2` degrees of freedom on each leg. This interval is deliberately
conservative.

For confidence statements restricted to proteins selected by balanced
partial-conjunction testing, let `R` be the number selected among `m` tested
proteins at family level `q`. Under independent protein vectors, BMEDIATOR
constructs each selected protein's simultaneous two-leg confidence set with
total noncoverage `qR/m`. Each leg receives half this budget by Bonferroni, and
the indirect-effect set is the corresponding corner-product range. The
arbitrary-dependence version replaces this budget by `qR/(mH_m)`, where `H_m`
is the `m`th harmonic number. These are false-coverage-rate adjustments in the
sense of Benjamini and Yekutieli. A selected protein with unstable effect
curvature receives an unbounded set and remains in the FCR denominator. The
claim is conditional on valid marginal adjusted-profile intervals; these are
not weak-instrument-exact confidence sets.

A positive status is a conditional causal conclusion under the stated
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

### Frozen analytical-calibration result

We froze balanced partial-conjunction BH and its decision criteria before
examining results, then evaluated a new simulation seed in 3,600 families of
500 proteins each. The design comprised 200 replicates of each mixed or stress
cell and 1,000 replicates of each global-null and least-favorable-null cell. At
a 5% target, mean FDP was 0.0376 in the broad balanced cell and 0.0298 in the
narrow balanced cell, with mean power 0.9720 and 0.5945. The least-favorable
and global-null family error rates were 0.009 and 0, with one-sided 95% Wilson
upper bounds 0.0154 and 0.0027. Mean FDP was 0.0043 under balanced
heterogeneity, 0.0082 with sparse outliers, and 0.0232 in the dense-dependence
stress cell. All prespecified criteria passed.

These results establish finite-grid operating characteristics, not the truth
of the causal assumptions. The dense-dependence cell does not prove PRDS or
arbitrary-dependence validity. Under directional pleiotropy, which violates
the balanced/InSIDE model, mean FDP was 0.5753. At the exactly proportional M5
boundary, mean FDP was 0.5008 because the observed summaries cannot identify
mediation separately from a direct effect proportional to pQTL strength. The
BH rule is therefore the leading confirmatory statistical procedure under its
stated assumptions, while regional H3/H4 evidence, heterogeneity, and external
replication remain separate causal-identification checks.

### Frozen post-selection effect-estimation result

We separately froze a validation of first-stage, second-stage, and indirect
effect inference after balanced partial-conjunction selection. The held-out
seed comprised 2,700 simulated proteomes of 500 proteins each, for 1.35 million
protein analyses. All 87 prespecified component criteria passed. For BH-selected
families, mean FCP was 0.01157 with strong instruments, 0.00165 with moderate
instruments, 0 with four instruments, and 0.02801 under balanced heterogeneity.
For BY-selected families, mean FCP was 0.00027 with strong instruments, 0.00114
under balanced heterogeneity, and 0.00033 under dense cross-protein dependence.
All selected confidence-set families were complete and none was unbounded.

Among true-M1 proteins in the four independent confirmatory cells, estimates
of both causal legs and their product were available in every analysis. The
minimum nominal causal-leg coverage was 0.94487, and the largest absolute bias
relative to empirical standard deviation among BH-selected M1 effects was
0.18997. In the declared-overlap cell, no finite selected confidence set was
emitted for any of 39,923 tested proteins, as prespecified.

This effect-validation result is conditional on valid balanced/InSIDE marginal
profile intervals. Under directional pleiotropy, BH mean FCP was 0.29491. The
near-threshold weak-instrument cell made no selections, and the current
finite-instrument Student-t intervals are not weak-instrument-exact. BY was
highly conservative in moderate and four-instrument cells, while same-sample
selection produced too few calls for confirmatory evaluation. These outcomes
define the scope of the effect-calibration claim rather than extending it to
violated assumptions.

## Methodological references

Wasserman L, Ramdas A, Balakrishnan S. Universal inference. *Proceedings of
the National Academy of Sciences* 117, 16880-16890 (2020).

Wang R, Ramdas A. False discovery rate control with e-values. *Journal of the
Royal Statistical Society: Series B* 84, 822-852 (2022).

Wang J, Gui L, Su WJ, Sabatti C, Owen AB. Detecting multiple replicating
signals using adaptive filtering procedures. *Annals of Statistics* 50,
1890-1909 (2022).

Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and
powerful approach to multiple testing. *Journal of the Royal Statistical
Society: Series B* 57, 289-300 (1995).

Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple
testing under dependency. *Annals of Statistics* 29, 1165-1188 (2001).

Benjamini Y, Yekutieli D. False discovery rate-adjusted multiple confidence
intervals for selected parameters. *Journal of the American Statistical
Association* 100, 71-81 (2005).

Wang S, Kang H. Weak-instrument robust tests in two-sample summary-data
Mendelian randomization. *Biometrics* 78, 1699-1713 (2022).
