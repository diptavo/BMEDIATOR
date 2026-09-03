# Factorized two-stage mediation model

## Status

The factorized engine is experimental. Invoke it with
`--structural-method factorized`. It does not replace the legacy M0-M5 output
yet. The legacy columns remain in the output for comparison, while all
confirmatory conclusions from this engine are in columns prefixed `factor_`.

The design fixes a structural defect in a mutually exclusive M1/M5 model:
mediation, a residual RF-to-outcome effect, and correlated pleiotropy can
coexist. The factorized model estimates these components separately.

## Identification

For protein `j`, BMEDIATOR retains three disjoint instrument sets:

- Set A: RF instruments outside the protein cis region. Only variants with an
  exact observed protein association are used by the factorized first stage;
  proxy-projected associations are excluded because their reported SE does not
  include LD-estimation uncertainty.
- Set B: cis-pQTL instruments that are not RF instruments. These identify the
  protein-to-outcome effect independently of the RF instruments.
- Set C: variants that are both RF instruments and cis variants. These are not
  used to estimate either confirmatory causal leg because their two roles are
  not separable without additional assumptions. They remain available for
  regional H3/H4 analysis.

The confirmatory interpretation requires at least two usable Set A variants,
at least two Set B variants, independent estimation errors for the GWAS used
in the two legs, and a regional shared-signal result. The instrument counts can
be changed, but lowering them weakens overidentification.

## Leg likelihoods

For either causal leg, let `x` be the vector of exposure associations, `y` the
outcome associations, and `R` the signed LD matrix from the reference panel.
For a candidate causal effect `b`, the errors-in-variables residual likelihood
is

```text
r(b) = y - b x
V(b) = S_y + b^2 S_x - b(C_xy + C_xy')
r(b) ~ N(0, V(b))
```

where `S_x = D_x R D_x`, `S_y = D_y R D_y`, and `D_x` and `D_y`
contain reported standard errors. If sampling-error correlation is supplied,
`C_xy[i,j] = corr_xy * se_x[i] * R[i,j] * se_y[j]`.

The maximum-likelihood estimate is obtained by deterministic one-dimensional
optimization. Its standard error is the inverse square root of the observed
likelihood curvature. This is used for effect estimation, not for the
confirmatory p-value.

For a fixed zero-centered normal effect prior with variance `W`, the Bayes
factor is

```text
BF_10 = integral L(b) Normal(b; 0, W) db / L(0).
```

BMEDIATOR evaluates this integral by deterministic Simpson quadrature. The
prior variances are the existing `--sigma2-beta1`, `--sigma2-beta2`, and
`--sigma2-beta3` values. They are fixed by default and recorded in `.hyp`; the
analyzed proteins do not calibrate their own priors.

The first leg fits RF-to-protein using Set A. The second leg fits
protein-to-outcome using Set B only. The residual RF-to-outcome effect is fit
after subtracting the estimated protein path from outcome associations at Set
A. This separation prevents correlated effects among RF instruments from
being treated as evidence for the protein-to-outcome leg.

## Analytical frequentist calibration

Each leg also has a null score test conditional on its selected exposure
associations. Under independent GWAS estimation errors and the Gaussian
summary-statistic model,

```text
z = x' S_y^-1 y / sqrt(x' S_y^-1 x).
```

This null test remains valid with weak exposure instruments because the
selected `x` values are conditioned on and selection occurred in a separate
GWAS. Residual overdispersion is estimated from the whitened regression and is
only allowed to increase the variance. A Student t reference with `n-1`
degrees of freedom is used, making the small-instrument test conservative.

Two-stage mediation has the composite null

```text
H0: beta1 = 0 OR beta2 = 0.
```

The intersection-union p-value is therefore

```text
p_conjunction = max(p_XM, p_MY).
```

Across proteins, BMEDIATOR applies the Benjamini-Yekutieli procedure to these
conjunction p-values. BY is valid under arbitrary dependence, which is useful
when proteins share instruments or correlated molecular measurements.

The p-value calibration does not use simulated truth, fitted empirical nulls,
or the observed cross-protein p-value distribution. Simulations are still
required to verify finite-sample behavior and power under departures from the
working assumptions.

## Regional LD gate

Regional H3/H4 evidence is not a p-value and is not included in the
intersection-union calculation. It is a separate Bayesian identification
gate:

- H4/shared supported: the protein and outcome regional signals are compatible
  with a shared causal variant.
- H3/distinct supported: distinct regional signals explain the associations,
  so mediation is rejected for confirmatory reporting.
- Ambiguous: two-stage evidence is reported but mediation remains unresolved.

Signed LD is used both in the leg likelihoods and in regional multi-signal
fine-mapping. Without an LD reference, the engine reports estimates and BFs but
cannot produce a conditionally identified mediation call.

## Interpretation

`factor_mediation_status` is the Bayesian evidence status. By default it
requires both leg-specific Bayes factors to be at least 10, adequate Set A and
Set B counts, an LD reference, no declared sample-overlap correlation, and a
shared regional configuration. `factor_min_log_BF` is the weaker of the two
log Bayes factors; it is an evidence summary, not a Bayes factor for the union
null.

`factor_frequentist_status` is based on the conjunction p-value and its BY
q-value. It uses the same instrument and regional gates. The Bayesian and
frequentist statuses may differ because the fixed prior and conservative
small-sample t test answer different questions.

`factor_pattern` reports which individual effects or residual-correlation
diagnostics cross the nominal `--factor-alpha` threshold. It is descriptive
and is not a mutually exclusive biological state.

## Fundamental boundary

If direct effects at Set B are exactly proportional to the protein effects,
then horizontal pleiotropy satisfies the same observable equation as a causal
protein effect. Multiple instruments do not solve this when every direct
effect is aligned with the same slope. BMEDIATOR therefore makes a conditional
mediation claim under the Set B exclusion restriction; it does not claim to
identify this exact aligned-pleiotropy boundary from summary statistics.

## Current limitations

- Nonzero overlap correlations are incorporated into effect likelihoods, but
  confirmatory status is `UNRESOLVED_SAMPLE_OVERLAP`. Selection-aware overlap
  inference is not implemented yet.
- The residual pleiotropy correlation is a de-noised diagnostic, not a
  confirmatory test and not part of the mediation p-value.
- Factorized mode skips legacy CAVI, but runtime and memory benchmarking on a
  complete proteome must be completed before it becomes the default.
- Broad simulations with realistic LD, instrument selection, winner's curse,
  sparse and directional pleiotropy, ancestry mismatch, and missingness remain
  required.
