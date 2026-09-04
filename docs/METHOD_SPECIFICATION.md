# Confirmatory Method Specification

## Scope

This document defines the estimand, assumptions, evidence layers, and release
criteria for the factorized BMEDIATOR method. It is normative for development:
software output must not make a stronger claim than this specification permits.

For risk factor `X`, molecular trait `M_j`, and outcome `Y`, the target path
coefficients are

```text
beta_XM,j: X -> M_j
beta_MY,j: M_j -> Y
indirect_j = beta_XM,j * beta_MY,j.
```

The product is an interventional indirect effect only under linearity, no
`X`-by-`M_j` interaction on the analyzed scale, valid instruments for both
legs, and compatible populations. Otherwise it is a two-stage genetic path
coefficient.

## Instrument Roles

- Set A contains risk-factor instruments outside the molecular cis region.
  Exact Set A molecular associations identify `beta_XM`.
- Set B contains cis molecular instruments that are not selected risk-factor
  instruments. Set B identifies `beta_MY` independently of Set A.
- Set C contains variants with both roles. Set C is excluded from confirmatory
  leg estimation and retained for regional diagnostics.

Factorized construction cross-clumps Set A and Set B using the configured LD
threshold. Every RF instrument inside the cis window belongs to Set C or is
excluded from confirmatory estimation, even if its molecular association does
not pass the cis-pQTL threshold. A cis-pQTL in LD above the threshold with any
RF instrument is also assigned to Set C, and the corresponding RF instrument
is removed from Set A. The largest retained cross-set `r²` is reported; a
threshold violation fails the identification gate.

At least three usable Set A and three usable Set B variants are required to fit
the confirmatory leg models. Three is the minimum that leaves residual
information after estimating a slope and an allele-oriented intercept. It does
not by itself prove the exclusion restriction.

## Summary-Statistic Model

For either leg, let `x` and `y` denote exposure and outcome association
estimates at the selected variants. Let the latent exposure associations be
`g`, signed LD be `R`, sampling covariance matrices be

```text
Sx = Dx R Dx,
Sy = Dy R Dy,
C  = Cov(error_x, error_y),
```

and let `a` denote direct effects on the leg outcome. The observation model is

```text
x = g + error_x,
y = beta g + a + error_y.
```

With a flat prior for `g` and `a ~ N(0, tau^2 I)`, integrating `g` conditional
on the observed `x` gives

```text
y | x, beta, tau ~ N(
    beta x,
    Sy + beta^2 Sx - beta(C + C') + tau^2 I
).
```

This is the interpretation of the errors-in-variables evidence likelihood. The
point estimator uses a generalized adjusted profile score, with an exact
adjusted-profile sandwich variance for independent instruments and profile
curvature for correlated instruments. Independent discovery statistics remain
required to avoid winner's-curse reuse in confirmatory inference. The release
analysis must report instrument-strength strata and three-sample sensitivity.

The additive `tau` model is a working model for balanced residual
heterogeneity. Directional pleiotropy is represented separately by `eta s`,
where `s_i=sign(x_i)`. This is an intercept after orienting exposure effects
positive. Exactly proportional pleiotropy remains nonidentified.

## Hypotheses

The two leg nulls are

```text
H_XM: beta_XM = 0,
H_MY: beta_MY = 0.
```

The statistical two-stage null is the union

```text
H_two-stage = H_XM union H_MY.
```

A valid intersection-union p-value is

```text
p_two-stage = max(p_XM, p_MY).
```

Rejecting this null establishes evidence for two nonzero genetic slopes. It
does not, by itself, establish biological mediation.

## Analytical Calibration

Each leg has two p-values calculated conditionally on selected exposure
associations. Exposure-only selection may use the same exposure estimates,
provided outcome information is not used and exposure/outcome errors are
independent. Both scores first GLS-project the
exposure design away from the allele-oriented intercept. Under independent
exposure/outcome errors and the declared Gaussian covariance, the
known-variance score has an exact standard-normal null reference while allowing
an arbitrary oriented intercept. A residual-scaled score estimates the ordinary
residual scale and uses a Student t reference with `n-2` degrees of freedom; its exact
reference additionally requires the declared covariance to be correct up to a
common scalar multiplier. It is not calibrated for arbitrary heterogeneous or
sparse direct effects. Neither score identifies direct effects exactly
proportional to instrument strength.

Three family-level procedures are reported separately:

1. BY adjustment of the residual-scaled `p_two-stage`, valid under arbitrary
   cross-protein dependence when every input p-value is valid under the scalar
   dispersion model.
2. BY adjustment of the strict Gaussian `p_two-stage`; this has the same
   dependence guarantee under the declared Gaussian covariance and
   allele-oriented-intercept mean model.
3. e-BH applied to fixed Student t density-ratio e-values. The null t statistic
   eliminates the common scalar dispersion and unrestricted allele-oriented
   intercept; this is also valid under arbitrary cross-protein dependence.

Two experimental analytical tracks address the observed power loss without
changing the preceding outputs:

4. A through-origin, residual-scaled score assumes mean-zero pleiotropy
   conditional on instrument strength (balanced pleiotropy/InSIDE). It has a
   Student t reference with `n-1` degrees of freedom under covariance
   proportional to the declared LD-aware covariance. BMEDIATOR applies BY to
   `max(p_XM_balanced,p_MY_balanced)` for arbitrary cross-protein dependence.
5. Ordinary BH is also applied to the balanced partial-conjunction p-value.
   It does not require independence between the leg p-values and controls FDR
   under independence or PRDS across the protein-level p-values.
6. AdaFilter-BH treats the two legs as a 2-of-2 partial conjunction, with
   `F=min(p_XM_balanced,p_MY_balanced)` and
   `S=max(p_XM_balanced,p_MY_balanced)`. This can reduce the multiplicity cost
   of proteins for which neither leg is plausible, but its guarantee requires
   independent base studies and independence or weak within-study dependence
   across proteins. It is reported as assumption-conditional.
7. A fixed p-to-e calibrator is applied directly to the balanced 2-of-2
   p-value, followed by e-BH. This retains the balanced/InSIDE and scalar
   covariance assumptions of the base p-value, but does not require
   independence between legs or across proteins.
8. A fixed Student t density-ratio e-value is also computed directly from each
   balanced residual-scaled leg statistic. The minimum leg e-value is passed to
   e-BH. This preserves the same mean-zero/InSIDE and scalar-dispersion model
   while permitting arbitrary dependence between legs and across proteins.

An additional strict e-value residualizes the oriented intercept and forms a
fixed mixture of Gaussian likelihood ratios. For `Z~N(0,1)` under the leg null,

```text
E_g(Z) = (1+g)^(-1/2) exp{g Z^2/[2(1+g)]},
g in {0.25, 1, 4, 16}.
```

The equally weighted mixture has null expectation one. The mediation e-value
is the minimum leg e-value and is passed to e-BH. This path permits arbitrary
cross-protein dependence but assumes the declared covariance is known; it does
not absorb residual overdispersion. See `docs/ANALYTICAL_CALIBRATION.md`.

The strict p-to-e sensitivity uses a fixed mixture of calibrators

```text
e_kappa(p) = kappa * p^(kappa - 1),  0 < kappa < 1,
e_mix(p) = sum_k w_k e_kappa_k(p),   sum_k w_k = 1.
```

For a super-uniform null p-value, each component and any prespecified mixture
have expectation at most one. The `kappa` values and weights must be frozen
before final validation. They may not be selected from the analysis p-value
distribution or from the final simulation labels.

No confirmatory p-value, e-value, or adjusted value is reported when causal-leg
GWAS sampling errors overlap unless a selection-aware overlap derivation is
implemented and validated.

## Bayesian Evidence

For each leg, a zero-centered normal prior on `beta`, a half-normal prior on
`tau`, and a zero-centered normal prior on the oriented intercept `eta` define
a four-model comparison: null, slope only, directional only, and both. The same
heterogeneity prior is integrated under the null and alternative. Fixed leg
and directional prior probabilities yield slope probabilities marginalized
over directional presence. After cross-clumping Set A and Set B, the
working-model two-stage probability is
`PP_two_stage = PP_XM * PP_MY`, and its posterior local FDR is
`1 - PP_two_stage`. Proteins are ordered by local FDR and the running mean is
the posterior expected FDP used for selection.

This product is exact only under the specified factorized likelihood. Disjoint
variant names do not establish independence when residual cross-set LD induces
correlated sampling errors. The reported maximum cross-set `r²` makes that
approximation auditable, and a threshold violation fails closed.

The default development priors are `P(XM)=0.50` and `P(MY)=0.25`; both are
recorded in `.hyp` and can be changed explicitly. They and the Student-t
alternative-density grid were frozen before the held-out validation seeded
`20261603`. Posterior FDR is a model-based
expectation, not a frequentist guarantee under prior or likelihood
misspecification. Results must include prior sensitivity; the Bayesian family
rule remains provisional because the held-out success applies to balanced
partial-conjunction BH, not to posterior FDR.

The heterogeneity Bayes factor compares the effect model integrated over
`tau > 0` with the corresponding `tau = 0` model. Strong Set B heterogeneity
withholds a mediation interpretation but does not negate evidence for two
genetic slopes.

Integrated evidence, balanced effect estimation, and directional sensitivity
analysis have separate numerical validity checks. Confirmatory effect estimates
use the balanced/InSIDE adjusted-profile slope and do not require an oriented
directional intercept to be identifiable. The joint slope/intercept fit is
reported separately as a sensitivity diagnostic. If balanced profile curvature
is invalid, the ordinary effect estimate is unavailable and the identification
gate returns `UNRESOLVED_EFFECT_ESTIMATION`.

For effect statements made after balanced partial-conjunction selection, the
software reports FCR-adjusted joint confidence sets. With `R` selections among
`m` tested proteins, the BH version assigns total noncoverage `qR/m` to each
selected protein; the BY version uses `qR/(mH_m)`. That budget is split between
the two leg intervals, and the indirect-effect interval is the range of their
four endpoint products. An unstable selected estimate receives an explicitly
unbounded set instead of being removed from the selected-family denominator.
These constructions are conditional on validity of the underlying marginal
profile intervals. BH-FCR requires independent protein vectors; the harmonic
version is the arbitrary-dependence procedure.

## LD and Causal Identification

Regional analysis asks whether molecular and outcome associations are
compatible with shared rather than distinct causal variants. H3 supports
distinct variants; H4 supports a shared variant. H4 is necessary compatibility
evidence but is not a mediation test because same-variant horizontal
pleiotropy has the same H4 signature.

Conditional protein and outcome signals are paired using their H3/H4 evidence.
Supported H4 edges are reduced to a maximum bipartite matching, so every
conditional signal is counted at most once. The evidence tiers are:

- one matched H4 signal: `UNRESOLVED_SINGLE_SHARED_SIGNAL`;
- at least two matched H4 signals with a stable Set B slope:
  `SUPPORTED_EXCLUSION_RESTRICTION_CONDITIONAL`;
- supported H3: `REJECTED_DISTINCT_REGIONAL_SIGNALS`;
- strong Set A or Set B heterogeneity: an explicit `UNRESOLVED_*_HETEROGENEITY`
  state;
- weak or conflicting regional evidence: unresolved.

Even infinitely many instruments cannot distinguish

```text
y_k = beta_MY x_k
```

from direct effects `a_k = beta_MY x_k` with causal `beta_MY = 0`. This aligned
pleiotropy boundary is observationally nonidentified. The strongest reportable
claim is therefore conditional on the Set B exclusion restriction; the method
must never label this boundary as empirically resolved.

## Required Outputs

Every molecular trait must report, separately:

- both leg estimates, uncertainty, projected p-values, Bayes factors, and e-values;
- indirect product estimate and uncertainty;
- Set A and Set B instrument counts and strength diagnostics;
- residual heterogeneity estimates and heterogeneity Bayes factors;
- directional-component BFs/posteriors and slope-intercept collinearity;
- maximum retained Set A/Set B cross-LD `r²`;
- balanced/InSIDE adjusted-profile effect estimates and uncertainty, with the
  joint directional-intercept fit reported separately as a diagnostic;
- H3/H4 evidence and independent matched-shared-signal count;
- `factor_two_stage_status`, which makes no mediation claim;
- `factor_mediation_status`, which applies identification gates and names the
  exclusion-restriction condition;
- factorized leg and two-stage posterior probabilities, posterior local FDR,
  cumulative posterior FDR, rank, and `factor_posterior_status`;
- family-level BY and e-BH decisions when their assumptions are met.
- balanced partial-conjunction BH q-values and statuses when cross-protein
  independence or PRDS is defensible.
- BH- and BY-adjusted selected-protein confidence sets, their per-protein
  noncoverage level, and explicit bounded/unbounded status.

Legacy M0-M5 probabilities may be emitted for compatibility but are not part
of confirmatory factorized inference.

## Release Criteria

The development suffix must remain until an untouched validation program
shows all of the following:

1. Null size and family FDR at or below the declared level under LD, weak
   instruments, cross-protein dependence, and valid sample configurations.
2. Explicit failure, not false certainty, at the aligned-pleiotropy boundary.
3. Controlled false mediation calls under distinct variants in LD, directional
   pleiotropy, sparse outliers, and LD-reference mismatch.
4. Bias and interval coverage reported for both legs and the indirect product.
5. Stable conclusions across prespecified effect and heterogeneity priors.
6. Agreement or explained disagreement with MR-RAPS, CAUSE, SMR/HEIDI,
   coloc-SuSiE, conventional two-step MR, and multivariable MR benchmarks.
7. Reproducible proteome-scale runtime and bounded memory use.
8. Replication across at least two molecular-QTL resources and multiple disease
   outcomes with ancestry-matched LD sensitivity analyses.

Passing simulations is necessary for release but cannot prove the exclusion
restriction in real data. Real-data conclusions must retain the
assumption-conditional label.

The frozen 1.8-million-analysis validation from commit `dc6eb00` satisfied the
prespecified FDR and power criteria for balanced partial-conjunction BH. This
completes the statistical-family criterion only under balanced/InSIDE leg
validity and cross-protein independence/PRDS. It does not complete the other
release criteria above, establish arbitrary-dependence control, or identify
exactly proportional pleiotropy.

A second held-out validation at commit `bc5d8f1` evaluated the reported
`beta1`, `beta2`, and indirect-effect confidence sets after proteome-wide
selection. The frozen design contained 2,700 families of 500 proteins (1.35
million analyses), and all 87 component criteria passed. BH mean FCP was
0.01157, 0.00165, 0, and 0.02801 in the strong, moderate, four-instrument, and
balanced-heterogeneity confirmatory cells. BY mean FCP was 0.00027 under strong
instruments, 0.00114 under balanced heterogeneity, and 0.00033 under dense
cross-protein dependence. Selected interval families were complete with no
unbounded sets. True-M1 marginal effect availability was 100%; causal-leg
coverage was at least 0.94487; and the largest absolute standardized bias among
BH-selected M1 effects was 0.18997. Declared overlap emitted no finite selected
sets for 39,923 tested proteins.

This result does not remove the release requirements above. Directional
pleiotropy produced BH mean FCP 0.29491, demonstrating failure when the
balanced/InSIDE condition is violated. The near-threshold cell made no
discoveries, so weak-instrument-exact inference remains unresolved. BY was too
conservative for a confirmatory claim in the moderate and four-instrument
cells, and same-sample selection yielded too few calls for a strong conclusion.
The complete frozen decision and provenance are in
`sim/results/factorized_effect_validation_20260906_v1/`.
