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

## How M1 and M5 are separated

M1 and M5 are distinguishable only because their nuisance structures are
restricted differently across the instrument sets. For a protein-specific
Set B instrument `l`, the outcome equation is

```text
Gamma_l = beta2 * alpha_l + phi_l + error_l.
```

Under M1, `beta2` is common across Set B instruments and `phi_l` is a sparse,
mean-zero direct effect independent of `alpha_l`. Under M5, `beta2 = 0`.
Correlated pleiotropy in M5 is represented by correlation between the
RF-to-protein residual `delta_k` and RF-to-outcome residual `psi_k` at
RF-associated Set A/C instruments. It is not allowed to create unrestricted
correlation between protein effects and outcome-direct effects at Set B
instruments. Thus multiple independent Set B signals can support a common
protein-to-outcome slope under M1 while the RF-instrument residual pattern can
support M5.

This restriction is essential. Suppose a proposed M5 model instead allowed a
protein-direct component `d` and an outcome-direct component `h` to have
arbitrary covariance. At a protein-specific component,

```text
M1: M = d,  Y = beta2 * d + h,  Cov(d,h) = 0
M5: M = d,  Y = h,              Cov(d,h) != 0.
```

M1 implies `Cov(M,Y) = beta2 Var(d)`. M5 implies
`Cov(M,Y) = Cov(d,h)`. Choosing the latter covariance to equal the former makes
the models observationally equivalent; their marginal variances can likewise
be matched. A Bayes factor between those formulations would therefore be
determined by prior scale and parameter-count choices, not by identified
mediation evidence. BMEDIATOR does not use that formulation.

The remaining identification claim is explicitly conditional. It requires at
least two independent RF-to-protein observations, at least one protein-specific
cis instrument, and the exclusion/independence assumptions above. A locus with
only a shared variant, or a model that permits outcome-direct effects
proportional to every protein effect, cannot distinguish M1 from M5 using these
summary statistics and must remain unresolved.

## Regional configuration model

Full mode retains all harmonized cis-region variants shared by the protein GWAS,
outcome GWAS, and PLINK reference. The default `ld-multisignal` method extracts
the regional genotypes from the PLINK panel, computes their signed correlation
structure, and iteratively selects conditionally associated signals for each
trait. For each retained signal it conditions on the other lead signals,
computes Wakefield approximate Bayes factors and a posterior credible set, and
then evaluates every protein-outcome signal pair under:

- H0: neither trait is associated
- H1: protein only
- H2: outcome only
- H3: both traits, distinct causal variants
- H4: both traits, the same causal variant

The `.regional` file reports all signal pairs. The `.mediation` regional fields
report the pair selected for the protein-level decision: a supported shared pair
takes priority, followed by the strongest supported distinct pair, followed by
the pair with the largest H3+H4 evidence. Resolution requires H3+H4 >=
`--regional-min-both` and the conditional shared or distinct probability to
exceed `--regional-min-shared`.

The default per-variant priors are `p_protein = 10^-4`,
`p_outcome = 10^-4`, and `p_shared = 10^-8`. Thus
`p_shared = p_protein * p_outcome`: before seeing the regional statistics, the
two causal-status indicators are independent. This removes the 1,000-fold H4
enrichment imposed by the commonly used `p_shared = 10^-5` setting. The
independence prior is conservative for true mediation but analytically
defensible for distinguishing shared signal from distinct variants in high LD;
it was chosen from the prior factorization, not fitted to simulation labels.

The LD matrix now enters the regional calculation directly through conditional
association statistics, signal separation, lead-pair LD, and credible-set-pair
LD. `--regional-method single` reproduces the previous marginal one-causal
calculation for compatibility; it should not be used as the primary analysis
in regions with allelic heterogeneity.

## Conditional identification assumptions

`LD_RESOLVED_SHARED_SIGNAL_ASSUMPTION_CONDITIONAL` means that the data support
M1 and H4 under all of these assumptions:

- The reference-panel LD is a sufficiently accurate estimate of the LD in both
  GWAS samples for conditional signal separation.
- The conditional signal threshold and maximum signal count do not omit a
  material protein or outcome signal.
- The retained RF instruments are independent at the configured LD threshold,
  the cis-only signal is sufficiently strong, and weak-instrument bias is
  negligible.
- The exclusion restriction holds: cis instruments affect the outcome through
  the measured protein, apart from modeled sparse pleiotropy.
- Same-variant horizontal pleiotropy at a protein-specific instrument is absent
  or sparse and independent of instrument strength. M5 represents correlated
  RF-instrument residuals; H4 alone cannot distinguish same-variant horizontal
  pleiotropy from mediation.
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
- `UNRESOLVED_NO_OUTCOME_SIGNAL`: no outcome signal passes the conditional
  regional signal threshold.
- `UNRESOLVED_NO_PROTEIN_SIGNAL`: no protein signal passes the conditional
  regional signal threshold.
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

Native conditional multi-signal inference is the full-mode default. It is a
summary-statistic conditional fine-mapping approximation, not an implementation
of SuSiE. A prior-matched SuSiE/coloc comparison on the bundled UKB-PPP and
deCODE examples agreed on the positive shared PCSK9 result but was more decisive
for distinct IL6R signals. The native method conservatively left IL6R ambiguous.
Scientific validation against broader architectures and ancestry-matched LD
panels remains required.
