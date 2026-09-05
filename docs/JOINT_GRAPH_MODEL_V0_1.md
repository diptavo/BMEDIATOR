# BMEDIATOR joint graph model JG-0.1

## Status and purpose

This document freezes the first reference specification for the new joint
Bayesian BMEDIATOR core. JG-0.1 is a research prototype used to establish the
likelihood, identification boundary, and agreement between independent R and
C++ implementations. It is not a production analysis mode and must not replace
the validated factorized sensitivity analysis until later validation is
complete.

The target is candidate-specific partial mediation of a risk factor `X`
through one molecular trait `M` to an outcome `Y`. The residual `X -> Y` path
is always allowed. It includes a biologically direct effect and pathways
through every mediator omitted from the candidate-specific model.

## Observed data and instrument roles

For variant `k`, the observed vector is

```text
bhat_k = (xhat_k, mhat_k, yhat_k)'.
```

Its three entries are GWAS effect estimates for `X`, `M`, and `Y`. JG-0.1 uses
three prespecified instrument roles:

- `A`: risk-factor anchor; large latent `X` variance and small residual
  molecular variance;
- `B`: molecular anchor; small latent `X` variance and large residual molecular
  variance;
- `C`: joint anchor; both latent variances are large.

All three roles enter one likelihood. They are not fitted as separate MR
regressions. The role labels determine reference prior variances, not observed
zero constraints.

The JG-0.1 reference likelihood assumes variants have already been LD-pruned,
that the three GWAS sampling errors are independent, and that effects use
compatible standardized scales. These restrictions isolate the graph-mixture
question for the first implementation. An LD- and sample-overlap-aware stacked
likelihood is required before production use.

## Structural model

For each variant, let

```text
g_k ~ N(0, vX_role)
d_k ~ N(0, vM_role)
e_k ~ N(0, vY)
```

be mutually independent latent direct genetic components. The structural
effects are

```text
x_k = g_k
m_k = a g_k + d_k
y_k = c g_k + b m_k + h_k lambda d_k + e_k.
```

The parameters are:

- `a`: `X -> M` path;
- `b`: global `M -> Y` path compatible with mediation;
- `c`: residual `X -> Y` path, always present;
- `lambda`: non-aligned shared-pleiotropy loading;
- `h_k`: a Bernoulli indicator with probability `q`.

When `h_k=0`, variant `k` follows the global structural graph. When `h_k=1`,
its residual molecular component additionally loads on `Y`. For `0 < q < 1`,
this creates a subset-specific non-Gaussian mixture pattern. It is the first
reference representation of non-aligned correlated pleiotropy.

The candidate indirect path coefficient is

```text
indirect = a b.
```

It is an interventional indirect effect only under the usual linear structural
assumptions and valid-instrument conditions. Otherwise it is a genetic path
coefficient.

## Marginal variant likelihood

For component `h` in `{0, 1}`, define

```text
T_h = [ 1          0              0
        a          1              0
        c + a b    b + h lambda   1 ].
```

With

```text
V_role = diag(vX_role, vM_role, vY)
E_k    = diag(seX_k^2, seM_k^2, seY_k^2),
```

integrating `(g_k, d_k, e_k)` gives

```text
bhat_k | h, theta ~ N(0, T_h V_role T_h' + E_k).
```

The per-variant likelihood is

```text
L_k(theta) = (1-q) phi_3(bhat_k; 0, Sigma_0k)
             + q phi_3(bhat_k; 0, Sigma_1k).
```

The complete reference likelihood is the product over LD-pruned variants.
This likelihood is trivariate and joint: every variant contributes its RF,
molecular, and outcome association simultaneously.

## Factorial graph states

JG-0.1 uses three protein-level indicators:

```text
zXM: a is free rather than fixed at zero
zMY: b is free rather than fixed at zero
zP:  lambda and q are active rather than fixed at zero
```

This produces eight states:

| State | zXM | zMY | zP | Interpretation |
| --- | ---: | ---: | ---: | --- |
| `S000` | 0 | 0 | 0 | no candidate molecular path |
| `S100` | 1 | 0 | 0 | RF-responsive molecular trait only |
| `S010` | 0 | 1 | 0 | molecular-outcome path not downstream of RF |
| `S110` | 1 | 1 | 0 | partial-mediation-compatible global paths |
| `S001` | 0 | 0 | 1 | subset shared pleiotropy only |
| `S101` | 1 | 0 | 1 | RF-responsive pleiotropic mimic |
| `S011` | 0 | 1 | 1 | molecular-outcome path plus pleiotropy |
| `S111` | 1 | 1 | 1 | partial mediation plus pleiotropy |

The residual path `c` is free in every state. Mediation never requires full
mediation or `c=0`.

The primary graph summaries are

```text
PP_XM              = P(zXM = 1 | data)
PP_global_MY       = P(zMY = 1 | data)
PP_nonaligned_P    = P(zP  = 1 | data)
PP_two_path        = P(zXM = 1, zMY = 1 | data)
PP_two_path_plus_P = P(S111 | data).
```

`PP_two_path` is mediation-compatible evidence conditional on excluding the
aligned-pleiotropy equivalence below. It is not an assumption-free probability
that biological mediation is true.

## Frozen reference priors and quadrature

The first reference implementation uses independent inclusion priors

```text
P(zXM=1) = 0.25
P(zMY=1) = 0.10
P(zP =1) = 0.10.
```

Active continuous parameters have independent priors

```text
a      ~ N(0, 0.70^2)
b      ~ N(0, 0.70^2)
c      ~ N(0, 0.175^2)
lambda ~ N(0, 0.70^2).
```

For `zP=1`, `q` has equal prior mass on

```text
q in {0.15, 0.35, 0.60}.
```

The reference latent variances are

| Role | vX | vM |
| --- | ---: | ---: |
| A | 0.0400 | 0.0025 |
| B | 0.0004 | 0.0400 |
| C | 0.0400 | 0.0400 |

and `vY=0.0025`. These fixed variances define the standardized reference
experiment; they are not acceptable production defaults for arbitrary GWAS
scales.

Marginal likelihoods use seven-node Gauss-Hermite quadrature for every active
normal parameter and exact summation over the three `q` values. State priors
are the products of the three Bernoulli inclusion priors. R and C++ must use
the same nodes, weights, summation order, and log-sum-exp stabilization.

Before locking the deterministic acceptance fixtures, an implementation sanity
check rejected five-node common-scale quadrature because it had no node near a
modest nonzero residual path and numerically aliased that path with the
pleiotropy mixture. The seven-node, parameter-specific rule above is the
numerical freeze used for all acceptance fixtures. No behavioral threshold was
changed.

## Identification statements

### First path

When `vX_role > 0`, `d_k` is independent of `g_k`, and sampling covariance is
known, the population covariance satisfies

```text
Cov(x_k, m_k) = a vX_role.
```

Thus `a` is identified from roles containing RF variation. Weak `vX`, few
variants, selection, or correlated direct molecular effects weaken or violate
this statement.

### Global molecular-outcome path versus residual RF path

After substitution,

```text
y_k = (c + a b) g_k + b d_k + h_k lambda d_k + e_k.
```

Variation in `d_k` that is not collinear with `g_k`, supplied primarily by B
and C variants, separates `b` from `c`. If every `d_k` is zero or proportional
to `g_k`, only `c + a b` is identified.

### Non-aligned shared pleiotropy

For `0 < q < 1` and nonzero `lambda`, the pleiotropy model is a mixture of two
trivariate normal distributions with different molecular-outcome covariance.
The global-path model is a single trivariate normal distribution for a fixed
role. Subject to nondegenerate variances and enough variants, these
distributions are distinct. Separation becomes weak as `q` approaches zero or
one, `lambda` approaches zero, or the number of informative variants falls.

### Exact aligned-pleiotropy equivalence

If `h_k=1` for every variant, a no-mediation model with

```text
b* = 0
lambda* = b
c* = c + a b
```

has exactly the same observation distribution as a mediation model with
`lambda=0`. Therefore exact aligned pleiotropy and a global molecular-outcome
slope are observationally equivalent in these summary statistics. No prior,
additional graph label, or inference algorithm changes this fact.

JG-0.1 consequently does not create a data-derived posterior that claims to
separate these two equivalent explanations. `PP_global_MY` and `PP_two_path`
must carry an explicit exclusion-condition label. An aligned-pleiotropy
sensitivity analysis is required in all later validation.

## Development acceptance tests

The deterministic development fixtures are not a calibration study. They test
whether the implemented model has the behavior implied by its own likelihood.

1. R and C++ state and aggregate posterior probabilities must agree to an
   absolute tolerance of `1e-8`.
2. A clean partial-mediation fixture must have `PP_XM >= 0.80`,
   `PP_global_MY >= 0.80`, and `PP_nonaligned_P <= 0.20`.
3. A non-aligned pleiotropic-mimic fixture must have `PP_XM >= 0.80`,
   `PP_global_MY <= 0.20`, and `PP_nonaligned_P >= 0.80`.
4. A coexistence fixture must have all three marginal mechanism probabilities
   at least `0.70`.
5. A null fixture must have `PP_XM <= 0.20` and `PP_global_MY <= 0.20`.
6. Before examining additional seeds, a development replication sweep is fixed
   at 50 replicates per mechanism, 60 variants per role, and seeds
   `20261001` through `20261050` plus scenario offsets of 0, 100, 200, and 300.
   At least 80% of replicates in each mechanism must satisfy its corresponding
   rules above.

Failure of one of these tests is a model or implementation result. Thresholds
must not be changed merely to obtain a pass.

## Required next extensions

JG-0.1 deliberately does not establish production validity. Subsequent model
versions must add and validate:

- signed-LD covariance across the stacked three-trait vector;
- declared cross-GWAS sampling covariance;
- selection-aware or externally estimated role variances;
- weak-instrument-robust inference;
- directional and sparse idiosyncratic pleiotropy components;
- prior sensitivity and posterior calibration;
- multi-signal regional latent configurations;
- proteome-wide multiple testing and post-selection effect summaries;
- comparison with the frozen factorized method and external competitors.
