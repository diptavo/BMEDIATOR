# Production readiness

## Release decision

`bmediator-joint` is the current implementation of the `JG-0.2.8` joint
Bayesian graph model. The older `bmediator` six-state and factorized modes are
retained for reproducibility and compatibility; they are not substitutes for
the joint model and must not be used to claim calibrated mediation.

Three different claims must remain separate:

1. **Software readiness** means the executable builds reproducibly, rejects
   invalid inputs, fails closed when numerical integration is unreliable, and
   preserves an auditable result for every manifest row.
2. **Model calibration** means a frozen analysis rule meets prespecified error
   and power criteria under the data-generating mechanisms that were tested.
3. **Scientific identification** means the assumptions needed to interpret
   `PP_two_path` as mediation are credible for the analyzed protein. Simulation
   calibration cannot prove those assumptions in real data.

## Implemented safeguards

- One joint trivariate likelihood uses Sets A, B, and C rather than combining
  two independently fitted MR regressions.
- The residual `X -> Y` path is free, so a protein may partially mediate the
  risk-factor effect.
- Signed within-block LD and declared cross-GWAS sampling covariance enter the
  likelihood directly.
- Global protein-outcome, sparse block-level pleiotropy, and directional
  pleiotropy indicators are fitted jointly in 16 factorial states.
- The sparse-state probability is integrated exactly under its Beta prior.
- Uncertain externally derived orientations are integrated within each block.
- At least three independent A-role and three independent B-role blocks are
  required by default.
- Invalid or non-positive-definite LD, excessive cross-block LD, invalid
  sampling correlations, duplicate variants/options, and nonfinite inputs are
  hard errors.
- User options may tighten but cannot relax LD, identification, optimizer, or
  numerical-reportability safeguards.
- Every graph state must optimize successfully. Posterior reporting is
  suppressed when the conservative normalized-posterior quadrature error
  exceeds `0.01` or a relevant state has a successive log-evidence difference
  above one.
- When aggregate posterior uncertainty remains above `0.01`, all state
  evidences are independently recomputed with Smolyak sparse-grid quadrature;
  posterior-influential states are then refined through sparse level 15 without
  changing the acceptance thresholds. Sparse levels, cancellation, and
  tensor-versus-sparse differences are reported.
- Every output records the model version, identification boundary, priors,
  numerical tolerances, role/block counts, state probabilities, and numerical
  diagnostics.
- The manifest runner preserves per-protein failures but fails closed for
  family-wide selection: posterior-FDR selection is unavailable until all
  manifest rows succeed.
- C++ and independent R reference tests, deterministic mechanism fixtures,
  manifest integration tests, and sanitizer checks are included.

## Mandatory analysis contract

A confirmatory analysis is valid only when all of the following were fixed
before outcome results were inspected:

- one genome build and effect-allele orientation across RF, protein, outcome,
  and LD data;
- ancestry-matched signed LD and defensible independent block definitions;
- A/B/C roles selected without reusing the analyzed outcome;
- role variances estimated in an independent, LD-pruned panel;
- orientation signs and probabilities estimated independently of the analyzed
  associations;
- sample overlap correlations declared, including zero only when justified;
- model priors and the protein family defined in advance;
- a complete manifest with no failed or suppressed proteins;
- the exact-aligned-pleiotropy exclusion accepted and reported.

The executable validates the resulting numerical files. It cannot establish
that the external study-design assertions are true.

## Permanent identification boundary

An outcome-direct effect exactly proportional to the protein genetic effect
has the same summary-statistic distribution as a global `M -> Y` effect. No
amount of LD information, sample size, or numerical calibration separates
those mechanisms without additional instruments or external assumptions.
Every result therefore states
`CONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY`.

## Remaining release priorities

### P0: required before routine real-data use

- Complete independent cluster validation of the new sparse-grid fallback and
  its cancellation safeguard. The successive-level posterior diagnostic is
  designed to be conservative but is not a mathematical upper bound on the
  true integration error.
- Freeze `JG-0.2.8` or its successor and pass a new, adequately powered family
  calibration on untouched seeds without relaxing the `0.01` posterior-error
  criterion. Small development simulations do not satisfy this requirement.
- Build and validate a cohort-specific pipeline from raw RF, pQTL, outcome,
  and reference-panel data to the harmonized A/B/C and signed-LD inputs. The
  current preparer validates and joins independently estimated role scales; it
  deliberately does not infer allele harmonization, roles, overlap, or blocks.
- Validate external role-scale and orientation estimation in genuinely
  independent data, including sensitivity to winner's curse and LD/reference
  mismatch.
- Run at least one preregistered real-data analysis and independent replication
  with a complete protein family.
- Add model-averaged posterior estimates and uncertainty for `a`, `b`, the
  residual path, and `a*b`; the current joint output is a mediator-evidence
  classifier, not an effect-estimation report.

### P1: required for a methods-paper performance claim

- Compare the frozen method against SMR/HEIDI, coloc-SuSiE, robust MR methods,
  two-step MR, and multivariable MR on identical inputs and estimands.
- Stress-test prior sensitivity, ancestry/LD mismatch, winner's curse, weak
  instruments, sample overlap, and broader non-aligned pleiotropic mechanisms.
- Demonstrate calibration in data-generating settings not designed around the
  fitted likelihood and report failure rates as well as successful fits.

### P2: packaging

- Decide whether to release the current C++ executables with R orchestration or
  expose the joint core through a formal R package API. The repository is not
  currently an R package.
- Add versioned example joint inputs generated from redistributable data and a
  machine-readable schema.

Passing a simulation gate is necessary evidence for the frozen model. It does
not, by itself, make a real-data mediation claim identified or make the method
publication-ready.
