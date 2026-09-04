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

The factorized partition uses both position and LD. Every RF instrument inside
the cis window is assigned to Set C or excluded from confirmatory estimation,
regardless of whether its protein association passes the cis-pQTL threshold.
A cis-pQTL with `r²` at or above `--clump-r2` to any RF instrument is also
assigned to Set C, and that RF instrument is removed from Set A. This check is
chromosome-wide in the reference panel rather than limited to the ordinary
clumping window. `factor_cross_set_max_r2` records the maximum retained Set
A/Set B cross-LD; a value above the configured threshold fails closed as
`UNRESOLVED_CROSS_SET_LD`.

The confirmatory interpretation requires at least three usable Set A variants,
at least three Set B variants, independent estimation errors for the GWAS used
in the two legs, and at least two conditionally independent regional
protein-outcome signal pairs classified as shared. A single H4 signal is
compatible with mediation but cannot distinguish mediation from same-variant
horizontal pleiotropy. The instrument counts can be changed, but lowering them
does not remove this identification boundary.

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

The point estimate uses a generalized adjusted profile score. Omitting the
residual-variance determinant while optimizing the slope removes the
weak-instrument attenuation induced by treating the noisy exposure estimate
as fixed. For independent instruments, BMEDIATOR uses the adjusted-profile
sandwich variance; correlated instruments use profile curvature. Effect
estimation is separate from the confirmatory null score.

To prevent heterogeneity from being absorbed into the causal slope, each leg
also contains an independent random-effect SD `tau`:

```text
V(b, tau) = V(b) + tau^2 I.
```

The prior on `tau` is half-normal and is fixed before analysis. Separate prior
scales are available for the XM, MY, and residual XY legs. The same nuisance
model is used under the null and alternative. For a fixed zero-centered normal
effect prior with variance `W`, the robust Bayes factor is

```text
BF_10 = integral integral L(b,tau) p(b) p(tau) db dtau
        / integral L(0,tau) p(tau) dtau.
```

BMEDIATOR evaluates a normalized, deterministic Simpson-weighted grid. This is
a proper finite discrete mixture approximation, so its weights sum exactly to
one. The effect prior variances are the existing `--sigma2-beta1`,
`--sigma2-beta2`, and `--sigma2-beta3` values. Heterogeneity prior scales are
set by `--factor-pleio-sd-xm`, `--factor-pleio-sd-my`, and
`--factor-pleio-sd-xy`. All are fixed and recorded in `.hyp`; the analyzed
proteins do not calibrate their own priors.

The first leg fits RF-to-protein using Set A. The second leg fits
protein-to-outcome using Set B only. The residual RF-to-outcome effect is fit
after subtracting the estimated protein path from outcome associations at Set
A. This separation prevents correlated effects among RF instruments from
being treated as evidence for the protein-to-outcome leg.

For each leg, BMEDIATOR also compares the effect model integrated over the
half-normal heterogeneity prior with the same effect model at `tau = 0`:

```text
BF_heterogeneity = integral integral L(b,tau) p(b) p(tau) db dtau
                   / integral L(b,0) p(b) db.
```

The XM and MY heterogeneity Bayes factors are overidentification diagnostics
for both causal legs. Strong evidence on either leg changes the mediation
status to an explicit unresolved-heterogeneity state; it does not erase the
reported two-stage association or prove a particular pleiotropic mechanism.

Each leg additionally compares four prespecified models: neither slope nor
directional component, slope only, an allele-oriented intercept only, and both.
If `s_i = sign(x_i)`, the directional component is `eta s`; after orienting all
exposure associations positive this is an MR-Egger-type intercept. `eta` has a
fixed zero-centered normal prior and is integrated analytically. Leg slope BFs
and posterior probabilities are averaged over absence or presence of this
competitor. Thus a common allele-oriented shift is not forced into the causal
slope. `factor_directional_collinearity_*` reports the GLS cosine between `x`
and `s`; a value near one means that the data contain little information for
separating the two mechanisms. The reported causal slope and directional
intercept are estimated jointly by adjusted profiling, so the effect estimate
does not silently revert to a zero-intercept model.
If the joint curvature is not positive definite, BMEDIATOR retains any
well-defined projected p-values, e-values, and integrated Bayesian evidence but
reports the effect estimate as unavailable. Such a result cannot receive a
supported mediation status and is labeled `UNRESOLVED_EFFECT_ESTIMATION` after
the regional gate. A zero-intercept estimate is never substituted silently.
The three component-model BFs versus the neither-component model are retained
in the output, allowing all four posterior weights to be recomputed under
alternative prespecified slope and directional priors.

## Bayesian posterior-FDR decision

The robust leg Bayes factors are converted using fixed causal-leg priors:

```text
PP_leg = P(slope present | data), marginalized over directional presence
PP_two_stage = PP_XM * PP_MY
local_FDR = 1 - PP_two_stage.
```

The product uses a working likelihood factorization across Sets A and B.
Different variant identifiers are not sufficient for exact independence:
sampling errors can remain correlated when variants are in LD. BMEDIATOR
therefore cross-clumps the sets and reports their maximum residual `r²`; the
product approximation is most defensible when that value is negligible. The
default development priors are `P(XM)=0.50` and `P(MY)=0.25`, controlled by
`--factor-prior-xm` and `--factor-prior-my`. They are fixed before analysis and
recorded in `.hyp`. Sorting local FDR values and taking their running mean gives
`factor_posterior_cum_fdr`; a protein passes the Bayesian family rule when this
is at most `--factor-alpha`.

This is analytical Bayesian calibration under the stated likelihood and fixed
priors. It does not use simulated labels or fit an empirical null. It is not a
frequentist guarantee under misspecified pleiotropy or priors. Regional LD,
shared-signal, and heterogeneity checks are applied after statistical
selection, and the strongest status remains assumption-conditional.

## Analytical frequentist calibration

Each leg also has a null score test conditional on the selected exposure
associations. Selection may use `x` itself, provided it uses no outcome
information and the exposure/outcome estimation errors are independent. Let
`s_i=sign(x_i)`, whiten `x`, `y`, and
`s` by `S_y`, and let `x_perp` be whitened `x` after GLS projection on `s`.
Under independent GWAS estimation errors and the Gaussian summary-statistic
model,

```text
z = x_perp' y / sqrt(x_perp' x_perp).
```

The nuisance projection makes the score invariant to any common
allele-oriented intercept. This null test remains valid with weak exposure
instruments because the selected `x` values are conditioned on. Independent
discovery `P_SELECT` statistics are nevertheless preferred for effect and
Bayes-factor interpretation because they reduce winner's-curse bias. The
residual scale is estimated from the whitened regression. A Student t reference
with `n-2` degrees of freedom is used for this residual-scaled version. Its
exact reference requires the declared covariance to be correct up to a common
positive Gaussian scale multiplier; arbitrary heterogeneous or sparse direct
effects are outside that guarantee.

BMEDIATOR also reports the known-variance Gaussian version of the projected
score. It is exact under the stated covariance model and permits an arbitrary
oriented intercept, but not residual heterogeneity beyond that covariance.
The residual-scaled t version permits common variance inflation but is not a
general pleiotropy-robust test.

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

For the separate balanced/InSIDE track, BMEDIATOR forms the 2-of-2 p-value
`p_balanced=max(p_XM_balanced,p_MY_balanced)` and applies the same frozen
calibrator mixture directly to that valid union-null p-value. The resulting
`factor_log_e_p2e_balanced_mediation` is passed to e-BH across proteins and
reported as `factor_e_q_p2e_balanced_EBH`. This route does not require
independence between the two legs or across proteins, but it does require the
balanced/InSIDE and scalar-dispersion assumptions used to obtain the base leg
p-values. See `docs/ANALYTICAL_CALIBRATION.md` for the full derivation and the
AdaFilter comparator.

The direct balanced e-value applies the prespecified Student-t
alternative-density mixture to each residual-scaled balanced leg statistic,
uses `min(E_XM_balanced,E_MY_balanced)` for the mediation union null, and
applies e-BH. Unlike the p-to-e conversion, this retains more of the score
magnitude; unlike AdaFilter, its FDR theorem permits arbitrary dependence. It
retains the balanced/InSIDE and common scalar-dispersion assumptions.

## Safe e-values and e-BH

Let `T` be the residual-scaled statistic above and `nu=n-2`. Under the leg null,
conditional on the selected exposure associations,

```text
T ~ t_nu
```

for every unrestricted common allele-oriented intercept and every positive
common Gaussian scale multiplier. Let `f_nu` be the central Student t density.
BMEDIATOR fixes a proper symmetric alternative density before analysis:

```text
g_nu(t) = (1/6) [
  sum_(mu in {2,4,6}) {f_nu(t-mu) + f_nu(t+mu)}/2
  + sum_(s in {2,4,8}) f_nu(t/s)/s
].
E_leg = g_nu(T) / f_nu(T).
```

Each shifted or scaled component integrates to one, hence `g_nu` is a density
and `E_null(E_leg | x)=1` exactly under the stated null model. The symmetric
shift components target moderate standardized effects and the scale components
retain evidence in the tails. The grids and equal weights are fixed in the
software and recorded in `.hyp`; no analyzed outcome or simulation label
chooses them. The guarantee requires independent exposure and outcome GWAS
estimation errors and correct summary-error/LD covariance up to a common
positive scalar. It does not cover arbitrary heterogeneous, sparse, or exactly
slope-aligned direct effects.
When either causal leg has declared sampling-error overlap, BMEDIATOR retains
the leg statistics for sensitivity analysis but reports both BY and e-BH
q-values as unavailable and sets the confirmatory status to
`UNRESOLVED_SAMPLE_OVERLAP`.

The mediation null is a union: either the XM null or the MY null is true.
Therefore

```text
E_mediation = min(E_XM, E_MY)
```

is an e-value under the complete mediation null, because it is bounded by the
valid e-value for whichever leg is null. BMEDIATOR applies base e-BH to these
mediation e-values. If `E_(1) >= ... >= E_(m)`, e-BH selects the largest `k`
such that

```text
E_(k) >= m / (alpha k).
```

This controls FDR under arbitrary dependence among valid e-values. The output
includes leg-specific log e-values, `factor_log_e_mediation`, and
`factor_e_q_EBH`. No simulated truth, empirical-null fit, or cross-protein
prior estimation enters this calculation.

The formal e-BH guarantee applies to the two-leg hypothesis family. The
regional H3/H4 assessment is reported as a subsequent identification gate; a
post hoc regional subset must not itself be described as retaining the e-BH
FDR guarantee without a dedicated structured e-value procedure.

## Regional LD gate

Regional H3/H4 evidence is not a p-value and is not included in the
intersection-union calculation. It is a separate Bayesian identification
gate:

- H4/shared supported: the protein and outcome regional signals are compatible
  with a shared causal variant. H4 alone is not evidence that the protein is
  the causal route from that variant to the outcome.
- H3/distinct supported: distinct regional signals explain the associations,
  so mediation is rejected for confirmatory reporting.
- Ambiguous: two-stage evidence is reported but mediation remains unresolved.

Shared signal-pair edges are reduced to a maximum bipartite matching so that
each conditional protein and outcome signal is counted at most once. This
count is reported as `regional_independent_shared_signals`. One matched signal
is reported as `UNRESOLVED_SINGLE_SHARED_SIGNAL`; at least two are required
for an overidentified, assumption-conditional call.

Signed LD is used both in the leg likelihoods and in regional multi-signal
fine-mapping. Without an LD reference, the engine reports estimates and BFs but
cannot produce a conditionally identified mediation call.

## Interpretation

`factor_two_stage_status` records whether both leg-specific Bayes factors pass
the evidence threshold after instrument and overlap checks. It deliberately
does not claim mediation. `factor_mediation_status` then applies the LD,
independent-shared-signal, and Set B heterogeneity gates. Its strongest state,
`SUPPORTED_EXCLUSION_RESTRICTION_CONDITIONAL`, remains conditional on the
absence of direct outcome effects proportional to the protein effects.
`factor_min_log_BF` is the weaker of the two log Bayes factors; it is an
evidence summary, not a Bayes factor for the union null.

`factor_frequentist_status` is based on the conjunction p-value and its BY
q-value. It uses the same instrument and regional gates. The Bayesian and
frequentist statuses may differ because the fixed prior and conservative
small-sample t test answer different questions.

`factor_ebh_status` reports whether the safe mediation e-value passed e-BH and
then records the same identification conditions. It is deliberately separate
from the Bayes-factor evidence status.

`factor_balanced_p2e_status` analogously combines balanced p-to-e/e-BH evidence
with the identification gates. The e-BH theorem applies to the pre-gate
two-leg family selected through `factor_e_q_p2e_balanced_EBH`; it does not imply
a separate FDR theorem for the post-gate status subset.

`factor_balanced_ebh_status` applies the same interpretation to the direct
balanced Student-t e-value. Its formal guarantee likewise belongs to the
pre-gate family selected by `factor_e_q_balanced_EBH`.

`factor_posterior_status` reports the fixed-prior posterior expected-FDR rule
and then applies the same regional and heterogeneity interpretation gates.
`factor_PP_two_stage` is evidence for two slopes under the working model; it is
not itself proof of biological mediation.

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
This remains true with multiple independent instruments: multiple instruments
detect departures from a common slope, not an alternative mechanism with the
same common slope.

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
- The default heterogeneity prior is intentionally broad. Bayes-factor power
  can be low with only two or three instruments, and prior-sensitivity analyses
  are required before biological interpretation.

## References

- Wasserman L, Ramdas A, Balakrishnan S. Universal inference. *PNAS* 2020.
  https://doi.org/10.1073/pnas.1922664117
- Wang R, Ramdas A. False discovery rate control with e-values. *JRSS B* 2022.
  https://doi.org/10.1111/rssb.12489
