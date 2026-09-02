# Changelog

## 1.2.0-dev

- Freeze scenario and effect priors by default; retain empirical-Bayes and
  data-adaptive local priors as explicit exploratory options.
- Make nuisance pleiotropy terms available in all six states so model evidence
  is compared on a common nuisance structure.
- Add full-mode regional shared-versus-distinct causal-configuration evidence.
- Use the analytical independence prior `p_shared = p_protein * p_outcome` by
  default instead of prespecifying enrichment for shared causal variants.
- Disable confirmatory mediation selection when LD resolution, regional
  association evidence, or multiple-instrument requirements are not met.
- Add explicit conditional-identification states and regression tests covering
  shared and distinct regional signals.
- Document the exclusion, single-causal-signal, ancestry, and sample-overlap
  assumptions required for interpretation.

## 1.0.0

- Initial GitHub-ready release.
- C++ command-line implementation of BMEDIATOR.
- Smoke-test dataset and `make test` workflow.
- Supporting scripts for preprocessing, Biowulf runs, aggregation, plotting, and simulation benchmarks.
