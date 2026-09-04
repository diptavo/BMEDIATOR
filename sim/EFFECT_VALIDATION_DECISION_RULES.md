# Frozen effect-estimation validation rules

These rules are specified before the held-out simulation is run. The validation
targets BMEDIATOR's reported first-stage effect (`beta1`), second-stage effect
(`beta2`), and indirect effect (`beta1 * beta2`) after proteome-wide selection.

## Estimands and selected family

- One protein is one confidence-set family containing `beta1`, `beta2`, and
  their product.
- Proteins are selected by the balanced partial-conjunction p-value
  `max(p_XM, p_MY)` followed by BH or BY at `q = 0.05`.
- The false coverage proportion (FCP) in one simulated proteome is the number
  of selected protein families for which at least one of the three reported
  intervals misses its truth, divided by the number selected. It is zero when
  no protein is selected. FCR is the expectation of this replicate-level FCP.
- A missing interval for a selected protein counts as noncoverage. It is never
  silently removed from the denominator. If an effect estimate is unstable,
  BMEDIATOR emits an explicitly unbounded confidence set; this is complete and
  covering but is recorded separately as uninformative.

For `R` selected proteins among `m` tested proteins, the BH confidence-set
family uses total per-protein noncoverage `0.05 R/m`. The two causal legs split
that budget by Bonferroni; the indirect-effect interval is the range of the
four endpoint products. The BY confidence-set family additionally divides the
per-protein budget by `H_m = sum(1/i, i=1,...,m)`. This follows the selected
confidence-interval construction of Benjamini and Yekutieli (2005), applied to
a simultaneous two-leg confidence set.

The guarantee is conditional on the marginal profile intervals having their
stated coverage. The current profile intervals are finite-instrument Student-t
approximations, not weak-instrument-exact confidence sets. Near-threshold
instruments are therefore a diagnostic stress test rather than part of the
confirmatory claim.

## Confirmatory cells

The independent-protein validation cells are:

- `strong_independent`
- `moderate_independent`
- `few_instruments_independent`
- `balanced_heterogeneity_independent`

For BH-FCR in each cell, all of the following must hold:

1. Mean replicate FCP is at most 0.05.
2. Mean FCP plus 1.96 Monte Carlo standard errors is at most 0.07.
3. Every selected protein has a complete interval family.
4. At most 1% of selected families are unbounded because effect curvature could
   not be estimated.
5. At least 90% of replicate proteomes make one or more discoveries, preventing
   a vacuous no-selection pass.

The same criteria apply to BY-FCR in `strong_independent` and
`balanced_heterogeneity_independent`. The moderate and four-instrument settings
may have no BY discoveries because the harmonic correction is intentionally
severe; their BY results are reported but cannot pass vacuously. Under
`dense_cross_protein_dependence`, the same criteria apply to BY-FCR. BH is
reported descriptively because its theorem requires independence or a justified
positive-dependence condition.

## Marginal estimation criteria

Among all true-M1 proteins in each confirmatory independent cell:

- At least 99% must have each effect estimate.
- Each nominal interval must have coverage of at least 0.925.
- Absolute bias divided by the empirical standard deviation must not exceed
  0.15 (strong), 0.20 (moderate), 0.25 (few instruments), or 0.20 (balanced
  heterogeneity).

Among BH-selected true-M1 proteins, the absolute standardized bias of each
effect must not exceed 0.25. This explicitly checks winner's-curse distortion
in point estimates even though FCR adjustment targets intervals rather than
point-estimate bias.

## Fail-closed and diagnostic cells

When causal-leg overlap is correctly declared (`declared_overlap`), BH/BY
selection and FCR intervals must not be emitted, and an
`UNRESOLVED_SAMPLE_OVERLAP` status must be present.

The following cells are descriptive and cannot make the confirmatory result
pass:

- `near_threshold_independent`: weak-instrument behavior.
- `same_sample_selection`: instrument selection and estimation use the same
  summary statistics.
- `unreported_overlap`: overlap exists in the data-generating process but is
  supplied to BMEDIATOR as zero.
- `directional_pleiotropy`: the balanced/InSIDE exclusion assumption is
  violated.

Failure in a diagnostic cell defines a limitation or a required methodological
repair; it is not evidence against a theorem whose assumptions exclude that
cell.

## References

- Benjamini Y, Yekutieli D. False discovery rate-adjusted multiple confidence
  intervals for selected parameters. *JASA*. 2005;100:71-81.
  doi:10.1198/016214504000001907.
- Wang S, Kang H. Weak-instrument robust tests in two-sample summary-data
  Mendelian randomization. *Biometrics*. 2022;78:1699-1713.
  doi:10.1111/biom.13524.
