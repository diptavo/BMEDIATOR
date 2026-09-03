# Changelog

## 1.2.0-dev

- Add an experimental factorized structural engine that estimates RF-to-protein,
  protein-to-outcome, residual RF-to-outcome, and pleiotropy components without
  forcing mutually exclusive M1/M5 labels.
- Add signed-LD errors-in-variables likelihoods, deterministic fixed-prior Bayes
  factors, conservative intersection-union tests, and BY correction across
  proteins for factorized analyses.
- Require separate Set A and Set B evidence and retain regional H3/H4 as an
  independent identification gate; sample-overlap analyses remain unresolved.
- Skip six-state CAVI in factorized mode and add a frozen 2,400-analysis pilot.
- Correct `.regional` row serialization so H3, H4, and downstream columns align
  with the header.

- Withdraw the experimental joint regional structural likelihood because its
  unrestricted M5 protein/outcome direct-effect covariance is nonidentified
  from the M1 protein-to-outcome path. Retain the Set A/B/C structural model and
  separate LD-aware H3/H4 regional diagnostic as the production analysis.
- Add an independent-dataset M1/M5 simulation study. The 36,000-run validation
  failed the production criterion because six-state M5 recall was 32.3% and
  pairwise probabilities were miscalibrated; document the release as research
  evaluation only.

- Freeze scenario and effect priors by default; retain empirical-Bayes and
  data-adaptive local priors as explicit exploratory options.
- Make nuisance pleiotropy terms available in all six states so model evidence
  is compared on a common nuisance structure.
- Add full-mode regional shared-versus-distinct causal-configuration evidence.
- Replace the default one-causal regional calculation with native LD-aware
  conditional multi-signal fine-mapping and per-signal-pair H0-H4 output.
- Add `.regional` output with lead variants, signal-pair posteriors, and
  credible-set cross-LD diagnostics; retain `--regional-method single` for
  compatibility.
- Use the analytical independence prior `p_shared = p_protein * p_outcome` by
  default instead of prespecifying enrichment for shared causal variants.
- Disable confirmatory mediation selection when LD resolution, regional
  association evidence, or multiple-instrument requirements are not met.
- Add explicit conditional-identification states and regression tests covering
  shared and distinct regional signals.
- Document the exclusion, conditional-signal, ancestry, and sample-overlap
  assumptions required for interpretation.

## 1.0.0

- Initial GitHub-ready release.
- C++ command-line implementation of BMEDIATOR.
- Smoke-test dataset and `make test` workflow.
- Supporting scripts for preprocessing, Biowulf runs, aggregation, plotting, and simulation benchmarks.
