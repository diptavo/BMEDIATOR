# Identification and LD Resolution

BMEDIATOR separates three questions:

1. Which structural scenario (M0-M5) best explains the RF, protein, and outcome
   summary statistics?
2. Do the regional protein and outcome associations arise from shared or
   distinct conditional signals?
3. Does the locus contain enough independent information to distinguish a
   mediated path from LD and modeled correlated pleiotropy?

`P_M1` answers the first question within the six-state model.
`regional_PP_shared`, `regional_PP_distinct`, and the `.regional` file address
the second. `P_mediator_ld_resolved` is equal to `P_M1` only when the regional
and instrument requirements for the third question pass; otherwise it is zero.
`P_mediator_identified` is a deprecated compatibility alias.

The probabilities are approximate model weights. The nonregional evidence is
a variational evidence lower bound (ELBO), and the joint regional evidence uses
a Laplace approximation. They are not assumption-free probabilities that a
biological mechanism is true.

## Joint LD-aware structural likelihood

In full mode, the default `--regional-method joint-ld` retains harmonized RF,
protein, and outcome statistics for the same variants in the protein cis
region. This three-trait intersection is distinct from the larger
protein-outcome intersection used for H3/H4, so missing RF rows cannot change
the shared-versus-distinct diagnostic. Signed LD is calculated from the PLINK
reference panel. For each trait, BMEDIATOR iteratively selects conditional lead
variants. The union of RF, protein, and outcome leads defines the regional
components.

For a component LD matrix `R`, marginal effect vector `b_t`, and standard-error
matrix `D_t` for trait `t`, BMEDIATOR forms the approximate joint effects

```text
theta_hat_t = R_lambda^-1 b_t
```

where `R_lambda = (1-lambda)R + lambda I`. The transformed sampling covariance
between traits `s` and `t` is

```text
Cov(theta_hat_s, theta_hat_t)
  = eta_st R_lambda^-1 D_s R_lambda D_t R_lambda^-1,
```

where `eta_st` is a prespecified summary-error correlation. The defaults are
zero, corresponding to nonoverlapping studies or negligible overlap effects.
The three supplied correlations must define a positive-definite matrix.

At component `j`, the latent direct genetic effects are `g_j` on RF, `d_j` on
protein, and `h_j` on outcome. The structural mapping is

```text
theta_Xj = g_j
theta_Mj = beta1 g_j + d_j
theta_Yj = (beta3 + beta1 beta2) g_j + beta2 d_j + h_j.
```

The latent direct effects have zero-mean normal priors with prespecified
variances. Under M5, `d_j` and `h_j` additionally have correlation
`regional_pleiotropy_rho`; under M0-M4 their covariance is zero. Integrating
`g`, `d`, and `h` gives a multivariate normal likelihood for all three traits
and all components. Active beta coefficients are optimized for each scenario,
and a local Hessian supplies the Laplace evidence and approximate coefficient
uncertainty.

The independent genome-wide block uses Set A instruments in the existing CAVI
model. Cis-only Set B and overlap Set C observations are removed from that block
when the joint regional likelihood is available, preventing the same cis data
from being counted twice. Because both blocks use the same beta priors, their
approximate posterior overlap is divided by the prior once when combining
evidence. M1 coefficient estimates combine the two normal posterior
approximations by the same rule.

At least two regional components are required before joint regional evidence is
allowed to affect structural model ranking. A one-component locus is reported
as `UNRESOLVED_SINGLE_COMPONENT`; its calculated regional values are diagnostic
only.

## Shared-versus-distinct diagnostic

The conditional leads are also fine-mapped separately for protein and outcome.
For each signal, BMEDIATOR conditions on the other leads, computes Wakefield
approximate Bayes factors and a credible set, and evaluates every
protein-outcome signal pair under:

- H0: neither trait is associated
- H1: protein only
- H2: outcome only
- H3: both traits, distinct causal variants
- H4: both traits, the same causal variant

The `.regional` file reports all pairs. The protein-level decision prioritizes
a supported shared pair, then a supported distinct pair, then the pair with the
largest H3+H4 evidence. Resolution requires H3+H4 to meet
`--regional-min-both` and H4/(H3+H4), or its complement, to meet
`--regional-min-shared`.

The default association priors are `p_protein = 10^-4`,
`p_outcome = 10^-4`, and `p_shared = 10^-8`. Thus the shared prior equals the
product of the trait-specific priors instead of imposing prior enrichment for
colocalization.

`--regional-method ld-multisignal` retains this conditional H0-H4 procedure but
does not add the joint structural regional likelihood. `--regional-method
single` retains the earlier one-causal-variant regional calculation. Both are
compatibility modes.

## What multiple components resolve

Multiple components provide overidentifying evidence because a single set of
`beta1`, `beta2`, and `beta3` must explain the cross-trait effect pattern across
independent genetic perturbations. For example, under mediation, independent
protein components should induce outcome effects proportional to their protein
effects through a common `beta2`, after accounting for RF and direct outcome
effects. Two distinct variants merely in LD need not satisfy that pattern once
their signed LD is included.

Multiple components do not make mediation assumption-free. A horizontal
pleiotropic mechanism that produces the same proportional effects at every
component can remain observationally equivalent to mediation. M5 tests one
prespecified correlated-pleiotropy alternative, but it cannot enumerate every
biological confounder. External functional, interventional, or replication
evidence is still needed for a causal biological claim.

## Conditional identification assumptions

An LD-resolved mediation report depends on all of the following:

- The reference-panel LD represents the LD in each GWAS sample.
- The conditional signal procedure has not omitted a material signal.
- RF instruments are independent at the configured threshold and sufficiently
  strong; weak-instrument bias is negligible.
- The exclusion restriction holds apart from nuisance pleiotropy represented
  by the six-state model.
- The M5 direct-effect correlation is a scientifically reasonable sensitivity
  model for correlated pleiotropy.
- Alleles and genome builds are aligned across all inputs.
- Supplied overlap correlations adequately represent correlated summary-error
  noise. Their default of zero is inappropriate when material sample overlap is
  expected and cannot otherwise be ignored.
- The normal random-effects priors and Laplace approximation are adequate for
  the number and magnitude of retained components.

These assumptions are not learnable from the same summary statistics alone.

## Output states

- `LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL`: the regional shared-signal
  and instrument eligibility checks pass. The strength of mediation support is
  still given by `P_M1`.
- `LD_DISTINCT_SUPPORTED`: distinct protein and outcome configurations are
  favored; confirmatory mediation selection is disabled.
- `LD_CONFIGURATION_AMBIGUOUS`: H3 versus H4 is not resolved.
- `UNRESOLVED_SINGLE_COMPONENT`: fewer than two regional components are
  available for overidentification.
- `UNRESOLVED_JOINT_MODEL_FAILURE`: the regularized joint likelihood could not
  be evaluated numerically.
- `UNRESOLVED_NO_OUTCOME_SIGNAL` or `UNRESOLVED_NO_PROTEIN_SIGNAL`: a trait has
  no retained conditional regional signal.
- `UNRESOLVED_INSUFFICIENT_*_INSTRUMENTS`: a shared regional signal is supported
  but the RF or cis instrument requirement does not pass.
- `UNRESOLVED_NO_REGIONAL_DATA`: unpruned regional statistics are unavailable.
- `UNRESOLVED_NO_JOINT_REGIONAL_DATA`: fewer than two variants are shared by
  all three GWAS and the reference panel.
- `UNRESOLVED_NO_LD_REFERENCE`: no LD reference was supplied.

Legacy pre-clumped mode cannot provide confirmatory regional resolution.
`--allow-unresolved-selection` restores exploratory selection without this gate
and should not be interpreted as identified mediation.
