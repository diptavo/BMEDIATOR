# BMEDIATOR joint graph model JG-0.2 series

## Status

`JG-0.2` is the first repair of the failed `JG-0.1` numerical prototype. It is
kept separate from the production `bmediator` command and from the historical
factorized method. Passing its development tests is necessary but not
sufficient for a production claim.

Patch `JG-0.2.1` adds exact within-block averaging over uncertain independent
orientation calls and applies the adaptive-evidence diagnostic only to graph
states with posterior probability above `1e-6`. The original `JG-0.2` run and
its failures remain reproducible at commit `efdeee8`.

The model targets candidate-specific partial mediation of risk factor `X`
through molecular trait `M` to outcome `Y`. A residual `X -> Y` path remains
free in every graph.

## Required inputs

For every selected variant, the analysis requires aligned effect estimates and
standard errors for all three traits, an A/B/C role, an LD-block label, and an
orientation fixed using independent discovery data and its probability of
being correct. It also requires:

- a signed LD matrix covering the variants in identical allele orientation;
- externally estimated latent variances `v_x`, `v_m`, and `v_y` for every row;
- declared sampling-error correlations `rho_xm`, `rho_xy`, and `rho_my`.

The supplied sampling correlations must be constant within one protein
analysis and form a positive-definite correlation matrix. Role variances should
be estimated using independent, LD-pruned summary statistics. The included R
estimator requires at least 30 variants per role; substantially larger panels
are preferable.

## Structural model

Within LD block `j`, let the independent causal-scale random effects be

```text
g_j ~ N(0, diag(v_x))
d_j ~ N(0, diag(v_m))
e_j ~ N(0, diag(v_y)).
```

For causal variants in the block,

```text
x_j = g_j
m_j = a g_j + d_j
y_j = c g_j + b m_j + h_j lambda d_j + eta s_j + e_j,
```

where `s_j` is the latent orientation vector informed by independent sign
probabilities. `JG-0.2.1` exactly averages over possible signs within each LD
block, with a maximum of 12 uncertain signs per block. The block-level
indicator `h_j ~ Bernoulli(q)` identifies non-aligned pleiotropic LD signals.
Moving `h` from SNPs to independent LD blocks makes the likelihood coherent
under LD and avoids an exponential sum over SNP-level allocations.

For signed block LD matrix `R_j`, the expected marginal GWAS associations are

```text
E(bhat_x) = 0
E(bhat_m) = 0
E(bhat_y | eta) = R_j eta s_j.
```

Define

```text
K_g = R_j diag(v_x) R_j'
K_d = R_j diag(v_m) R_j'
K_e = R_j diag(v_y) R_j'
u   = c + a b
w_h = b + h lambda.
```

Conditional on `h`, the genetic covariance blocks are

```text
Cov(x,x) = K_g
Cov(x,m) = a K_g
Cov(x,y) = u K_g
Cov(m,m) = a^2 K_g + K_d
Cov(m,y) = a u K_g + w_h K_d
Cov(y,y) = u^2 K_g + w_h^2 K_d + K_e.
```

Sampling covariance is added using the reported standard errors, signed LD,
and cross-GWAS correlations. For example,

```text
Cov(error_xi,error_mk) = rho_xm R_ik se_xi se_mk.
```

The block likelihood is the exact two-component normal mixture

```text
(1-q) N_3n(mu_0, Sigma_0) + q N_3n(mu_1, Sigma_1).
```

Independent block likelihoods are multiplied. Cross-block absolute LD above
`0.05` causes a hard failure; smaller ignored correlations are reported.

## Factorial states

Four binary indicators produce 16 states:

```text
zXM:          a is active
zMY:          b is active
zSparse:      lambda and q are active
zDirectional: eta is active.
```

The state name `Sabcd` records these indicators in that order. Primary outputs
are `PP_XM`, `PP_global_MY`, `PP_sparse_P`, `PP_directional_P`, `PP_any_P`, and
`PP_two_path`. Mediation-compatible evidence is `PP_two_path`, conditional on
the unchanged exact-alignment exclusion.

## Priors and adaptive evidence

Default inclusion probabilities are `0.25`, `0.10`, `0.10`, and `0.10`.
Normal prior standard deviations are `0.70` for `a`, `b`, `lambda`, and `eta`,
and `0.175` for `c`. Active `q` has a `Beta(2,2)` prior and is optimized on the
logit scale with the Jacobian included.

Every state is optimized from multiple deterministic starts. Its marginal
likelihood is calculated with mode-centered three-node adaptive
Gauss-Hermite quadrature using the numerical posterior Hessian. This fixes the
off-grid failure of fixed prior-centered quadrature. The executable fails
without a posterior if:

- any of the 16 state optimizations does not converge;
- the posterior Hessian requires a ridge;
- adaptive and ordinary Laplace log evidence differ by more than one log unit
  in any state with posterior probability above `1e-6`;
- the LD, scale, orientation, or sampling-correlation inputs are invalid.

The frozen optimizer limit is 1,500 iterations with a `1e-6` convergence
tolerance. These values were set after a two-replicate runner smoke exposed a
nonconverged coexistence state under the earlier `500` and `1e-7` development
defaults; the fitted posterior was unchanged when convergence was attained.

All priors can be supplied in a tab-delimited key/value configuration file and
are written into the result. This supports prespecified prior-sensitivity
analysis without recompilation.

## External scale estimation

`research/estimate_joint_graph_v02_scales.R` estimates role-specific `v_x` and
`v_m`, a shared `v_y`, and nuisance path coefficients by maximum likelihood.
The input must come from an independent estimation panel and must be LD-pruned.
Using the selected discovery estimates would reintroduce winner's curse and is
not permitted for confirmatory interpretation.

## Identification boundary

Exact aligned pleiotropy remains observationally equivalent to a global
protein-outcome slope. `JG-0.2` does not claim otherwise, and every output
retains

```text
CONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY.
```

The directional state protects against a declared orientation-aligned mean; it
does not identify arbitrary pleiotropic functions. Sparse-state information is
carried by the number of independent LD blocks, not the raw number of SNPs.

## Current production boundary

`JG-0.2` is not production-ready until it passes the frozen new-seed matrix,
larger held-out calibration, prior sensitivity, LD/reference mismatch,
winner's-curse, missing-data, and competitor benchmarks. It must not be wired
into the main executable before those gates are passed.
