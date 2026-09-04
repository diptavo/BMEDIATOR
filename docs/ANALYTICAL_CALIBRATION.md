# Analytical calibration tracks

## Status

These outputs are experimental in version 1.2.0-dev. They use fixed
mathematical reference distributions and do not estimate a truth-dependent
calibration map from simulations or real data. The direct balanced
Student-t/e-BH and balanced p-to-e/e-BH rules passed the frozen FDR checks in
the tested identifiable cells, but both failed the prespecified power targets.
AdaFilter retained substantially more power but failed its all-null FDR
endpoint. Ordinary balanced partial-conjunction BH is undergoing a frozen
new-seed validation. None is designated as a production primary method.

## Hypothesis

For each protein, BMEDIATOR tests the two causal legs

\[
H_{XM}: \beta_{XM}=0, \qquad H_{MY}: \beta_{MY}=0.
\]

The mediation null is the union

\[
H_{med}=H_{XM}\cup H_{MY},
\]

so mediation requires rejection of both leg nulls. This is a 2-of-2 partial
conjunction problem. It is distinct from the regional shared-versus-distinct
signal analysis and does not, by itself, prove that the exclusion restrictions
hold.

## Balanced/InSIDE score

For one leg, let `x` be the independently selected exposure associations, `y`
the outcome associations, and `V` their LD-aware sampling covariance. After
whitening by `V`, the balanced working model is

\[
V^{-1/2}y = \beta V^{-1/2}x + \epsilon,
\qquad E(\epsilon\mid x)=0.
\]

The through-origin score is studentized by the residual scalar dispersion after
fitting the slope. Under Gaussian errors with covariance proportional to `V`,
the resulting statistic has a Student t distribution with `n-1` degrees of
freedom. This track assumes mean-zero pleiotropy conditional on instrument
strength, usually described as a balanced-pleiotropy/InSIDE assumption.

The leg outputs are `factor_p_XM_balanced` and `factor_p_MY_balanced`. Their
2-of-2 partial-conjunction p-value is

\[
p_{med}=\max(p_{XM},p_{MY}).
\]

`factor_balanced_conjunction_q_BY` applies Benjamini-Yekutieli correction and
therefore targets FDR control under arbitrary cross-protein dependence, subject
to validity of the two base p-values.

## Balanced Student-t e-value

The residual-scaled balanced statistic has a central Student t distribution
with `n-1` degrees of freedom under its leg null. BMEDIATOR evaluates the same
prespecified proper alternative-density mixture used by the directional safe
e-value: symmetric location shifts `{2,4,6}` and scale factors `{2,4,8}` with
equal weights. The density ratio has null expectation one for every positive
common scalar dispersion.

For mediation, BMEDIATOR takes the minimum of the two leg e-values. This is a
valid e-value for the union null because at least one leg is null and the
minimum is bounded by that null leg's e-value. e-BH then controls FDR under
arbitrary dependence. Outputs are `factor_log_e_XM_balanced`,
`factor_log_e_MY_balanced`, `factor_log_e_mediation_balanced`,
`factor_e_q_balanced_EBH`, and `factor_balanced_ebh_status`.

This construction is fully prespecified and does not require independence
between studies or proteins. It does require independent instrument selection,
no causal-leg sample overlap, valid LD-aware covariance up to a common scalar,
and mean-zero/InSIDE direct effects conditional on instrument strength.

## Balanced p-to-e calibration

The balanced conjunction p-value is also converted using the prespecified
mixture

\[
e(p)=\frac{1}{4}\sum_{\kappa\in\{0.10,0.25,0.50,0.75\}}
\kappa p^{\kappa-1}.
\]

Every component has expectation at most one for a super-uniform null p-value.
Because `max(p_XM,p_MY)` is super-uniform under the union null whenever either
leg p-value is valid, this construction does not require independence between
the two legs. BMEDIATOR applies e-BH to these e-values across proteins; e-BH
controls FDR under arbitrary dependence of valid e-values. The outputs are
`factor_log_e_p2e_balanced_mediation`, `factor_e_q_p2e_balanced_EBH`, and
`factor_balanced_p2e_status`.

This is a dependence-robust analytical track when the balanced/InSIDE and
scalar-dispersion assumptions are defensible. In the frozen one-million-
analysis run it passed the tested FDR criteria but had 0.0969 broad-signal
power and zero narrow-signal power at 0.05. Simulations test the assumptions
and finite-sample behavior; they do not estimate or tune the calibrator.

The formal e-BH guarantee applies to selection by
`factor_e_q_p2e_balanced_EBH` in the two-leg hypothesis family. The subsequent
regional identification gate adds biological interpretation but is not itself
part of the e-BH theorem; `factor_balanced_p2e_status` must therefore be read as
an exclusion-restriction-conditional status rather than as a separately proven
FDR-controlled regional subset.

## Balanced partial-conjunction BH

The same valid union-null p-value `max(p_XM,p_MY)` can be entered directly into
the Benjamini-Hochberg procedure across proteins. This does not require the two
leg p-values to be independent: under the mediation union null, at least one
leg p-value is super-uniform, and their maximum is therefore super-uniform.
BH controls FDR when the resulting protein-level p-values are independent or
satisfy positive regression dependency on each null statistic (PRDS).

The output is `factor_balanced_conjunction_q_BH`; the corresponding gated
interpretation is `factor_balanced_bh_status`. The BH theorem concerns the
pre-gate family of two-leg hypotheses. The later regional H3/H4 gate is a
causal-identification diagnostic, not part of that FDR theorem. For arbitrary
cross-protein dependence, use `factor_balanced_conjunction_q_BY`; for designs
that justify stronger leg and cross-protein independence assumptions,
AdaFilter is retained as a separate comparison.

## AdaFilter partial-conjunction FDR

For each protein define

\[
F=\min(p_{XM},p_{MY}), \qquad S=\max(p_{XM},p_{MY}).
\]

After ordering proteins by `S`, the AdaFilter-BH adjustment count for rank `j`
is the number of proteins with `F <= S[j]`. The adjusted value is the reverse
cumulative minimum of `S[j] * adjustment_count[j] / j`. The implementation
follows Wang et al., *Annals of Statistics* (2022),
[doi:10.1214/21-AOS2139](https://doi.org/10.1214/21-AOS2139).

AdaFilter can be substantially more powerful than applying BH or BY directly
to `max(p_XM,p_MY)`, but its finite-sample FDR theorem assumes independence of
all base p-values. Its asymptotic theorem permits weak dependence within a
study, while retaining independence between studies. BMEDIATOR therefore
reports `factor_balanced_conjunction_q_AdaFilter` and
`factor_adafilter_status` as assumption-conditional outputs. Dense correlated
protein measurement error is an explicit validation stress case.

## Information-adaptive e-value

The directional sensitivity model first projects out the allele-oriented
intercept `s=sign(x)`. With whitened vectors `wx` and `wy`, define

\[
r_x=(I-P_s)w_x, \qquad r_y=(I-P_s)w_y,
\qquad Z=\frac{r_x^T r_y}{\sqrt{r_x^T r_x}}.
\]

Under the strict null, conditional on independently selected `x`, `Z` is
standard normal when `V` is the correct outcome covariance and sampling errors
do not overlap across the exposure and outcome studies. For a fixed
signal-to-noise prior variance `g`,

\[
E_g(Z)=(1+g)^{-1/2}
       \exp\left\{\frac{g Z^2}{2(1+g)}\right\}
\]

is the likelihood ratio of a `N(0,1+g)` alternative to the `N(0,1)` null, so
its null expectation is exactly one. BMEDIATOR averages the prespecified grid
`g={0.25,1,4,16}`. The prior scale depends on exposure-side information but not
on outcome evidence, preserving conditional validity.

For two legs, `E_med=min(E_XM,E_MY)` is an e-value for the union mediation null:
whichever leg is null bounds the minimum. e-BH is then applied across proteins,
which controls FDR under arbitrary dependence of valid e-values. Outputs are
`factor_log_e_XM_adaptive`, `factor_log_e_MY_adaptive`,
`factor_log_e_mediation_adaptive`, and `factor_e_q_adaptive_EBH`.

This track permits an unrestricted allele-oriented intercept but assumes known
sampling covariance rather than an unknown residual overdispersion. It is
expected to be conservative when instrument magnitudes are nearly collinear
with their signs because projecting the intercept then removes most slope
information.

## Identification boundary

Analytical calibration cannot identify a causal slope against direct effects
that are exactly proportional to instrument strength. If

\[
y=\beta x+\epsilon
\]

and an alternative pleiotropic mechanism produces the same mean vector
`beta*x`, the summary-data likelihoods are identical. Multiple instruments
resolve mediation only when they provide overidentifying variation: distinct
instrument sets, sufficiently varied strengths, nonproportional residual
patterns, or independently supported regional signals. Exact proportional
pleiotropy remains an explicitly unresolved boundary.
