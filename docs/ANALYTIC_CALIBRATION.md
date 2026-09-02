# Analytical sensitivity calibration

BMEDIATOR's posterior local FDR is model-based. It is not guaranteed to be a
frequentist FDR when LD leakage, horizontal pleiotropy, or sample overlap lies
outside the fitted likelihood. `scripts/calibrate_bmediator_analytic.py`
provides a label-free sensitivity analysis; it does not learn a mapping from
simulation truth.

## Null and p-value

Let `alpha` denote the RF-to-protein effect and `beta` the
protein-to-outcome effect. For a prespecified bias bound `Delta >= 0`, define

```
H0: alpha = 0 OR |beta| <= Delta
H1: alpha != 0 AND |beta| > Delta.
```

The second component treats any effect no larger than `Delta` as potentially
explained by LD leakage or pleiotropy. A conservative Wald p-value is

```
p_beta(Delta) = 2 Phi(-max(0, (|beta_hat| - Delta) / se_beta)).
p_mediation   = max(p_alpha, p_beta(Delta)).
```

The maximum is an intersection-union p-value for the composite mediation
null. The script adjusts these p-values with Benjamini-Yekutieli (BY) by
default, which permits arbitrary dependence among proteins. BH is available
when its positive-dependence assumptions are defensible.

## Choosing the bias bound

`Delta` is a scientific sensitivity parameter, not a tuning parameter to be
chosen from simulated labels. Prefer protein-specific bounds from regional
conditional analysis. For cis-instrument effects `a`, outcome effects `y`,
and precision matrix `W`, suppose the unmodelled outcome component satisfies
`|u_k| <= U_k`. The induced slope bias is bounded by

```
Delta = sum_k |(W a)_k| U_k / (a' W a).
```

`U_k` can be bounded using an ancestry-matched LD matrix and conditional
outcome associations from variants outside the protein credible set. Report
results across several defensible bounds when a single upper bound cannot be
justified. Colocalization alone is not enough: a shared causal variant can
represent either vertical mediation or horizontal pleiotropy.

Example with one global sensitivity bound:

```bash
python3 scripts/calibrate_bmediator_analytic.py \
  --input analysis.mediation \
  --output analysis.analytic.tsv \
  --bias-bound 0.20 \
  --adjustment BY
```

Example with locus-specific bounds:

```bash
python3 scripts/calibrate_bmediator_analytic.py \
  --input analysis.mediation \
  --output analysis.analytic.tsv \
  --bounds-file protein_bias_bounds.tsv \
  --adjustment BY
```

The bounds file is tab-delimited and contains `Protein` and `bias_bound`.

## Stress-test audit

The sensitivity test was evaluated on the 7,000-replicate technical-report
stress suite. At a nominal 5% threshold, the standalone BY procedure gave:

| Delta | Overall FDR | Worst-cell FDR | Mean cell power |
|------:|------------:|---------------:|----------------:|
| 0.20 | 0.0273 | 0.0632 | 0.2972 |
| 0.21 | 0.0209 | 0.0442 | 0.2653 |
| 0.22 | 0.0157 | 0.0308 | 0.2387 |
| 0.23 | 0.0114 | 0.0200 | 0.2161 |
| 0.24 | 0.0085 | 0.0137 | 0.1981 |
| 0.25 | 0.0061 | 0.0102 | 0.1835 |

`Delta=0.21` was the smallest tested bound that kept every stress cell below
5% under BY. This validates the procedure's behavior for that simulation
scale; it does not justify 0.21 as a universal real-data bound.

## Scope and limitations

The current implementation uses the IVW estimates already present in a
`.mediation` file. It assumes their normal approximation and standard errors
are valid. For confirmatory analyses, use independent instrument-selection
and estimation samples and account for sample overlap.

The intended confirmatory extension is a SNP-level Anderson-Rubin/Fieller
test using pQTL and outcome covariance plus the LD matrix. That extension is
needed for formal weak-instrument robustness. No method using only the same
marginal RF, protein, and outcome associations can distinguish mediation from
an unrestricted proportional pleiotropic path; an explicit bound, exclusion
restriction, or additional data is mathematically necessary.
