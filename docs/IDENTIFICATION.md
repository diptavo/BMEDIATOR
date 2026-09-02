# Identification and LD Resolution

BMEDIATOR separates two questions that must not be conflated:

1. Does the six-state structural model favor M1 over M0-M5?
2. Is the observed protein-outcome association compatible with one shared
   regional causal signal rather than two distinct causal signals in LD?

`P_M1` answers the first question. It is structural model support, not by itself
an identified probability of mediation. `P_mediator_ld_resolved` is nonzero only
when the second question passes, at least two independent observed RF-to-protein
instruments are available, and a cis-only protein instrument exists. It remains
conditional on the assumptions below. `P_mediator_identified` is retained as a
deprecated compatibility alias and must not be read as assumption-free proof.

The `P_M*` values are normalized weights constructed from scenario priors and
variational evidence lower bounds. They are approximate model probabilities,
not exact marginal-likelihood posterior probabilities. Their operating
characteristics must be reported from validation studies separately from the
analytical H3/H4 calculation.

## Regional configuration model

Full mode retains all harmonized cis-region variants shared by the protein GWAS,
outcome GWAS, and PLINK reference. For each variant it computes Wakefield
approximate Bayes factors for the protein and outcome associations. It then
integrates over the standard single-causal-variant configurations:

- H0: neither trait is associated
- H1: protein only
- H2: outcome only
- H3: both traits, distinct causal variants
- H4: both traits, the same causal variant

The output reports `regional_PP_distinct` (H3), `regional_PP_shared` (H4), and
`regional_shared_given_both` = H4 / (H3 + H4). Resolution requires both
H3 + H4 >= `--regional-min-both` and the conditional shared or distinct
probability to exceed `--regional-min-shared`.

The default per-variant priors are `p_protein = 10^-4`,
`p_outcome = 10^-4`, and `p_shared = 10^-8`. Thus
`p_shared = p_protein * p_outcome`: before seeing the regional statistics, the
two causal-status indicators are independent. This removes the 1,000-fold H4
enrichment imposed by the commonly used `p_shared = 10^-5` setting. The
independence prior is conservative for true mediation but analytically
defensible for distinguishing shared signal from distinct variants in high LD;
it was chosen from the prior factorization, not fitted to simulation labels.

This model uses the complete pattern of marginal associations across the
region. Under the single-causal-variant assumption, that pattern contains the
information needed to compare a shared causal variant with distinct variants
in LD; an LD matrix does not enter the approximate-BF formula directly. The
PLINK panel is still required for allele alignment, instrument clumping, proxy
handling, and LD QC.

## Conditional identification assumptions

`LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL` means that the data support
M1 and H4 under all of these assumptions:

- At most one causal protein variant and one causal outcome variant occur in the
  tested region, or the region has already been conditioned into single-signal
  components.
- The retained RF instruments are independent at the configured LD threshold,
  the cis-only signal is sufficiently strong, and weak-instrument bias is
  negligible.
- The exclusion restriction holds: cis instruments affect the outcome through
  the measured protein, apart from modeled sparse pleiotropy.
- Same-variant horizontal pleiotropy is absent or adequately represented by M5
  and the nuisance terms. H4 alone cannot distinguish this from mediation.
- GWAS alleles and genome builds are aligned, and the LD panel matches the study
  ancestry.
- Correlated estimation error from participant overlap between the protein and
  outcome GWAS is negligible. The current regional model does not estimate that
  covariance.

These are scientific assumptions, not quantities that can be learned from the
same summary statistics without additional design information. Violating them
can produce a confident but incorrect result.

## Output states

- `LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL`: shared regional signal,
  at least two RF-to-protein instruments, and a cis-only signal; mediation
  remains conditional on the exclusion assumptions.
- `LD_DISTINCT_SUPPORTED`: distinct protein and outcome causal configurations
  are favored; mediation selection is disabled.
- `LD_CONFIGURATION_AMBIGUOUS`: H3 versus H4 is not resolved.
- `UNRESOLVED_WEAK_REGIONAL_EVIDENCE`: there is insufficient evidence that both
  traits are associated in the region.
- `UNRESOLVED_INSUFFICIENT_*_INSTRUMENTS`: shared signal is supported but the
  RF or cis instrument requirement is not met.
- `UNRESOLVED_NO_REGIONAL_DATA`: unpruned regional statistics are unavailable.
- `UNRESOLVED_NO_LD_REFERENCE`: no LD reference was supplied.

Legacy pre-clumped mode cannot provide confirmatory regional resolution.
`--allow-unresolved-selection` restores the old selection behavior and should
be used only for exploratory compatibility analyses.

## Multiple causal signals

The present regional calculation is not a multi-signal fine-mapping method. A
region with allelic heterogeneity must be conditioned or fine-mapped into
single-signal components before confirmatory interpretation. Native multi-signal
regional inference is planned but should not be claimed for this release.
