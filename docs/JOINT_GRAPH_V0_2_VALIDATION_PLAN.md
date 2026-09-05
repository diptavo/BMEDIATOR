# JG-0.2 frozen development validation plan

## Freeze

This plan is committed before the full `JG-0.2` simulation matrix is run. The
seeds start at `20264001` and have not appeared in `JG-0.1` development. Each
cell has 50 replicates. R generates data and invokes the C++ evaluator with up
to eight local workers.

## Posterior patterns

The following per-replicate patterns are fixed:

| Pattern | Criterion |
| --- | --- |
| `null` | `PP_two_path <= 0.20` and `PP_global_MY <= 0.20` |
| `mediation` | `PP_two_path >= 0.80` and `PP_any_P <= 0.20` |
| `sparse` | `PP_sparse_P >= 0.80` and `PP_global_MY <= 0.20` |
| `directional` | `PP_directional_P >= 0.80` and `PP_global_MY <= 0.20` |
| `two_pleiotropy` | both pleiotropy PPs at least `0.70` and `PP_global_MY <= 0.30` |
| `mediation_sparse` | `PP_two_path >= 0.70` and `PP_sparse_P >= 0.70` |

Every confirmatory cell must match its pattern in at least 80% of replicates.
Diagnostic cells have no pass threshold and cannot be used to tune the frozen
defaults.

## Cells

| Cell | Configuration | Pattern |
| --- | --- | --- |
| `matched_null` | 20 independent LD blocks | null |
| `moderate_mediation` | off-grid `a=b=0.40` | mediation |
| `weak_mediation` | `a=b=0.20`, SE `0.05` | diagnostic |
| `sparse_pleiotropy` | `lambda=0.70`, `q=0.35`, `b=0` | sparse |
| `directional_pleiotropy` | `eta=0.40`, `b=0` | directional |
| `two_pleiotropy` | sparse and directional components, `b=0` | two_pleiotropy |
| `mediation_sparse` | mediation plus sparse pleiotropy | mediation_sparse |
| `declared_overlap_null` | nonzero sampling correlations supplied correctly | null |
| `undeclared_overlap_null` | same errors but correlations reported as zero | diagnostic |
| `ld_null` | within-block AR(1) LD `rho=0.70` | null |
| `ld_mediation` | same LD with moderate mediation | mediation |
| `high_scale_declared_null` | fourfold latent variances supplied correctly | null |
| `high_scale_misspecified_null` | fourfold truth reported at reference scale | diagnostic |
| `low_scale_mediation` | quarter-scale truth supplied correctly | diagnostic |
| `few_blocks_sparse` | six independent pleiotropy blocks | diagnostic |
| `near_aligned_pleiotropy` | `q=0.90`, `b=0` | diagnostic |
| `exact_aligned_pleiotropy` | `q=1`, `b=0` | diagnostic |
| `orientation_70pct` | directional truth with 70% correct orientations | diagnostic |
| `skeptical_prior_mediation` | moderate mediation with lower inclusion priors | mediation |
| `diffuse_prior_null` | null with higher inclusion probabilities and wider priors | null |

The skeptical-prior mediation cell reuses the moderate-mediation datasets, and
the diffuse-prior null cell reuses the matched-null datasets. Only priors
change in these paired comparisons.

The matrix reports optimizer and adaptive-evidence diagnostics for every run.
Any hard numerical or input failure is retained as a failed replicate.

## Interpretation

Passing this matrix would show that the specific `JG-0.1` implementation flaws
have been repaired under matched and selected stress conditions. It would not
establish proteome-wide posterior calibration or justify a production release.
A later frozen, larger, held-out grid must use independent seeds and Biowulf
batch or interactive compute nodes, never the login node.
