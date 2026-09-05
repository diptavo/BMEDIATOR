# JG-0.1 development results

## Scope

These results cover the first four joint-core development steps only:

1. freeze the causal graphs, likelihood, priors, and outputs;
2. state the identification conditions and prove the aligned-pleiotropy
   equivalence;
3. implement a transparent base-R reference model and an independent C++
   evaluator;
4. test whether the frozen model separates clean partial mediation from
   non-aligned subset pleiotropy under its own data-generating assumptions.

They are implementation and model-behavior tests, not calibration evidence for
real GWAS analysis.

## Locked fixture results

Each fixture contained 60 independent variants in each of Sets A, B, and C.
The residual RF-to-outcome path was nonzero in every fixture. The R and C++
implementations compared every state posterior, aggregate posterior, and
posterior mean.

| Fixture | PP XM | PP global MY | PP non-aligned pleiotropy | PP two-path |
| --- | ---: | ---: | ---: | ---: |
| Null | 0.13223 | 0.04834 | 0.04835 | 0.00639 |
| Partial mediation | 1.00000 | 1.00000 | 0.04835 | 1.00000 |
| Non-aligned pleiotropic mimic | 1.00000 | 0.04980 | 1.00000 | 0.04980 |
| Mediation plus pleiotropy | 1.00000 | 1.00000 | 1.00000 | 1.00000 |

The maximum R/C++ absolute difference was `8.53e-13`, below the frozen `1e-8`
tolerance. A direct numerical test also confirmed that the clean global-slope
parameterization and its exactly aligned no-mediation reparameterization have
equal likelihood. The output consequently carries
`CONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY` rather than claiming to resolve
that boundary.

## New-seed development sweep

Before examining the sweep, the design was fixed at 50 replicates per
mechanism, 180 variants per replicate, and at least 80% correct posterior
patterns in every mechanism. Results were:

| Mechanism | Replicates passing | Rate | Mean PP XM | Mean PP global MY | Mean PP pleiotropy |
| --- | ---: | ---: | ---: | ---: | ---: |
| Null | 50/50 | 1.00 | 0.13223 | 0.04834 | 0.04836 |
| Partial mediation | 50/50 | 1.00 | 1.00000 | 1.00000 | 0.04843 |
| Non-aligned pleiotropy | 50/50 | 1.00 | 1.00000 | 0.05193 | 1.00000 |
| Coexistence | 50/50 | 1.00 | 1.00000 | 1.00000 | 1.00000 |

The standalone C++ process had a median elapsed time of approximately 0.07
seconds per 180-variant protein on the development machine. This is encouraging
but is not a formal proteome-scale runtime benchmark.

## What has been established

- The model is genuinely joint at the variant likelihood level: RF, molecular,
  and outcome associations enter one trivariate density.
- Sets A, B, and C all contribute to every graph comparison.
- The residual RF path coexists with every graph, so the target is partial
  rather than full mediation.
- Factorial states permit mediation and non-aligned pleiotropy to coexist.
- Under matched strong-signal simulations, the subset-mixture pattern is
  distinguishable from the global molecular-outcome slope.
- The efficient C++ calculation reproduces the transparent R calculation.

## What has not been established

- No type-I-error, FDR, coverage, or posterior-calibration claim has been made.
- The fixed role variances are reference values on standardized scales; they
  have not been estimated from real summary statistics.
- The prototype assumes LD-pruned variants and independent GWAS sampling
  errors.
- Weak instruments, winner's curse, sparse outliers, directional pleiotropy,
  LD mismatch, and sample overlap have not been added to this likelihood.
- The strong matched simulations are deliberately easy. Performance near the
  identification boundaries is unknown.
- Exact aligned pleiotropy remains observationally indistinguishable from a
  global molecular-outcome slope.
- The prototype is not wired into the main `bmediator` executable.

The next methodological phase is therefore model extension and adversarial
validation, not application to real proteins.
