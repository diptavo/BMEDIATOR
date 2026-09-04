# Validation Status

## Historical single-signal LD stress test

On 2026-08-17, the compatibility `single` regional method was tested with
`sim/run_regional_ld_stress.py` on Biowulf. Each cell used 500 replicates,
20 cis variants generated from an AR(1) haplotype model, 512 LD-reference
samples, two independent RF instruments, noisy protein/outcome summary
statistics, fixed inference priors, and the analytical regional independence
prior `p_shared = p_protein * p_outcome = 10^-8`.

| Scenario | Realized causal-variant LD | LD-resolved call rate | Distinct-supported rate |
|---|---:|---:|---:|
| Distinct causal variants | 0.19 | 0.0% | 97.6% |
| Distinct causal variants | 0.41 | 0.0% | 97.6% |
| Distinct causal variants | 0.71 | 4.6% | 89.8% |
| Distinct causal variants | 0.87 | 7.4% | 82.8% |
| Shared mediation | 1.00 | 84.2%-86.6% at low/moderate regional LD | 0.0% |
| Shared mediation | 1.00 | 62.6% at target LD 0.90 | 15.6% |
| Shared mediation | 1.00 | 20.2% at target LD 0.98 | 62.8% |

The target-LD parameter controls correlation among neighboring regional
variants. The shared scenario's causal variant is necessarily correlated 1.0
with itself; the decreasing call rate reflects reduced fine-mapping resolution
as nearby variants become nearly interchangeable.

For comparison, the enriched `p_shared = 10^-5` prior produced distinct-causal
false-call rates of 22.0% at realized LD 0.71 and 36.8% at realized LD 0.87.
The independence prior was derived analytically rather than fitted to labels
and substantially removed this prior-driven H4 preference.

## Identification boundary

The same-variant horizontal-pleiotropy arm generated summary statistics that
are observationally equivalent to shared mediation. Its call rates tracked the
shared-mediation arm, as theory requires. This is not a calibration failure
that can be repaired by another data-independent cutoff: discrimination
requires an exclusion assumption, multiple independent cis signals with a
multi-signal model, or external functional or interventional evidence.

An experimental joint regional random-effects model was withdrawn on
2026-09-03 because it assigned M5 unrestricted covariance between
protein-direct and outcome-direct component effects. That covariance duplicates
the observable implication of `beta2` under M1, so the resulting M1-versus-M5
weights were prior-driven and are invalid. In particular, the experimental
PCSK9 results that shifted support from M1 to M5 must not be used. The legacy
compatibility model remains the Set A/B/C CAVI model plus the separate
`ld-multisignal` H3/H4 regional diagnostic; it is not production-calibrated.

## Native multi-signal checks

The default `ld-multisignal` implementation is covered by a deterministic
genotype-based fixture containing two independent shared signals in one region
and a separate distinct-signal region. The test requires recovery of two
signals per trait, shared diagonal signal pairs, distinct off-diagonal pairs,
and rejection of the distinct-causal region.

The legacy implementation has also been compared with a SuSiE/coloc
reference analysis for BMI, CHD, and five proteins measured by each of UKB-PPP
and deCODE, using the same `p1 = 10^-4`, `p2 = 10^-4`, and `p12 = 10^-8`
priors. Both methods supported a shared PCSK9-CHD signal and found no resolved
CHD signal for ANGPTL3 or VEGFA. For IL6R, SuSiE/coloc supported distinct
signals while the native approximation remained conservative and reported an
ambiguous configuration. For NPPB, SuSiE did not retain an outcome credible
set, while the native method retained a weak pair but also reported it as
ambiguous. This is an implementation cross-check, not proof of calibration or
numerical equivalence to SuSiE.

## Dedicated M1-versus-M5 validation

On 2026-09-03, the legacy Set A/B/C structural model was evaluated in
36,000 independent protein simulations. Each protein used a separate synthetic
RF, protein, and outcome dataset; this avoids the invalid practice of mixing
protein-specific RF-outcome effects in one outcome GWAS. The grid contained 50
replicates and 2,000 observations per truth state per cell.

| Cell | Pair AUC | Pair accuracy | M1 sensitivity | M5 specificity | Six-state accuracy | Six-state M1 recall | Six-state M5 recall |
|---|---:|---:|---:|---:|---:|---:|---:|
| Set B = 2 | 0.908 | 0.824 | 0.845 | 0.804 | 0.530 | 0.668 | 0.391 |
| Set B = 4 | 0.976 | 0.901 | 0.956 | 0.847 | 0.626 | 0.852 | 0.400 |
| Set B = 8 | 0.999 | 0.914 | 0.999 | 0.828 | 0.658 | 0.924 | 0.391 |
| Weak second stage, Set B = 4 | 0.797 | 0.721 | 0.688 | 0.755 | 0.263 | 0.359 | 0.167 |
| M5 residual correlation = 0.5 | 0.956 | 0.842 | 0.948 | 0.737 | 0.612 | 0.847 | 0.377 |
| 50% sampling-error correlation | 0.985 | 0.933 | 0.966 | 0.899 | 0.696 | 0.893 | 0.500 |
| Set C-heavy | 0.661 | 0.622 | 0.639 | 0.600 | 0.062 | 0.089 | 0.036 |
| All seven cells | 0.911 | 0.834 | 0.873 | 0.794 | 0.492 | 0.662 | 0.323 |

Pair metrics condition on M1 versus M5 using
`P_M1/(P_M1+P_M5)`. Six-state metrics use the actual M0-M5 argmax. The strong
difference between these columns matters: the model often ranks the true M5
above M1 but assigns still greater support to M2 or M0. More Set B instruments
improve M1/M5 ranking but do not repair M5 identification against all six
states. Pair scores were defined for 94.6% of observations overall but only
62.1% in the Set C-heavy cell; six-state metrics include every observation.

The exactly aligned-pleiotropy boundary generated M5 direct effects
proportional to protein effects and was called M1 in 96.7% of pairwise
comparisons (95% Wilson interval 95.8%-97.4%). This is the expected
nonidentifiability result, not a false-positive rate that can be calibrated
away. The no-Set-B boundary also remained unstable.

The pairwise score was not calibrated as a probability. Its aggregate ten-bin
expected calibration error was 0.089, and several middle and upper bins were
overconfident. Full tables are retained in
`sim/results/m1_m5_identification_20260903/`.

This study therefore fails a production release criterion for M1/M5 inference.
The current model contains useful ranking information in favorable Set B
architectures, but `P_M1` and `P_M5` must not be described as calibrated
probabilities, and M5 is not reliably recovered against the complete six-state
model. Further likelihood development is required before a non-development
release or a strong methodological claim.

## Experimental factorized pilot

The replacement factorized engine was subjected to a frozen 2,400-analysis
pilot on 2026-09-03: 10 replicates, 20 independent datasets per truth state,
four states (M1, M2, M4, M5), and three cells. Unlike the six-state model, the
target was the composite two-stage null rather than forced state
classification.

| Cell | M1 conjunction power at 0.05 | Largest non-M1 rejection at 0.05 | M1 both-BF>10 | Largest non-M1 both-BF>10 |
|---|---:|---:|---:|---:|
| Identified, Set B=4 | 38.5% | 0.5% | 81.5% | 0.5% |
| Balanced Set B pleiotropy | 68.0% | 2.5% | 85.5% | 4.0% |
| Weak second stage | 1.3% | 0.7% | 0.0% | 0.7% |

The frequentist leg tests use a conservative small-instrument t reference;
the fixed-prior Bayes factors retain materially more power in the two
well-identified cells. The full frozen table is in
`sim/results/factorized_pilot_20260903/summary.tsv`.

The complete pilot was independently rebuilt and rerun under GCC on Biowulf
as jobs 28961869 and 28961870. Its summary was byte-identical to the frozen
local table. This establishes implementation reproducibility across the local
Apple Clang and Biowulf GCC builds for that frozen revision.

The pilot Bayes factors above predate the explicit heterogeneity nuisance
model. They remain an audit record for commit `34d74be`, but they are not
performance estimates for the current robust Bayes-factor implementation.

This pilot establishes that separating Set A and Set B repairs the specific
M1/M5 conflation in these independent-instrument simulations. It does not
validate the complete method. Realistic LD, sample overlap, winner's curse,
directional pleiotropy, ancestry mismatch, regional H3/H4 error, and
proteome-wide FDR remain open validation requirements. The factorized mode is
therefore experimental and the development release warning still applies.

## Robust likelihood and batched calibration smoke

A subsequent batched smoke analyzed all proteins in each replicate jointly so
that BY and e-BH operated over the simulated protein family. The initial
fixed-slope Bayes factor showed 16%-23% two-leg BF support in M2/M5 arms with
sparse Set B outliers. This failure motivated a model change: every leg now
contains a heterogeneity SD with the same fixed half-normal prior integrated
under the null and alternative. In a targeted 50-protein rerun of the same
sparse-outlier generator, M2/M5 two-leg BF support fell to 0%; M1 support was
40% under the intentionally broad default heterogeneity prior. These counts
are diagnostic only and are too small for a performance claim.

A historical 1,000-protein least-favorable composite-null smoke used 500 M2 proteins
(XM present, MY null) and 500 M4 proteins (XM null, MY present). It was rerun
after replacing a numerical additive-heterogeneity denominator with the exact
null supremum from the multiplicative-overdispersion e-value model. At nominal
0.05, the null leg rejected in 0.6% and 0.8%, respectively; neither BY nor
e-BH selected a protein. Mean safe e-values for the null legs were 0.210 and
0.206, below the required expectation bound of one. This confirms conservative
behavior in one configuration, not universal calibration.

That null-supremum construction was subsequently retired because it had
negligible power. Its results are retained only as a development audit record;
they do not validate the current fixed Student-t density-ratio e-value.

The durable batched runner and summaries are
`sim/run_factorized_task.py` and
`sim/summarize_factorized_calibration.py`. The smoke configurations include
signed instruments, near-threshold selection, directional and sparse Set B
pleiotropy, declared sample overlap, no-Set-B loci, and the aligned-pleiotropy
boundary.

## Proteome-scale factorized development calibration

The historical Biowulf development grid at
`factorized_calibration_main_20260903_v3` completed 2,600 batched tasks and
1,050,000 protein analyses without failed tasks. It included ten
classification cells and eight 500-protein family-calibration cells, with 100
and 200 replicates per cell, respectively. These data were used for method
development and are not the held-out validation of the final rule.

The robust raw conjunction test had 95.6% power and 4.6% mean FDR in the clean
baseline family, but mean FDR increased to 13.8% with directional pleiotropy,
18.9% with balanced pleiotropy, and 27.6% when mediation was rare. BY and the
safe likelihood e-BH rule made no discoveries in any cell. Applying fixed
analytical p-to-e calibrators to the robust p-values also made no discoveries;
the small-instrument robust p-values were therefore too weak for
proteome-wide e-BH rather than merely paired with a poor likelihood e-value.

A development prior grid converted the robust leg Bayes factors to posterior
probabilities. With `P(XM)=0.50` and `P(MY)=0.25`, mean empirical FDR was at
most 3.5% across all eight cells. Mean power ranged from 47.4% in the baseline
cell to 8.6% under balanced pleiotropy and zero near the selection threshold.
The high quantile of replicate FDR reached 25% in the balanced-pleiotropy cell
because very few traits were selected, so this does not establish uniform
frequentist FDR control. The prior pair has now been frozen for a new-seed
held-out run; it will not be changed in response to that run.

A 1,800-protein local smoke of the strict Gaussian score showed why it is only
a sensitivity analysis. Strict BY achieved 98.6% power with 1.6% mean FDR in
the clean cell, but mean FDR increased to 9.9% with directional pleiotropy and
27.4% with sparse outliers. The fixed strict p-to-e/e-BH rule showed the same
pattern (0%, 5.1%, and 21.9% mean FDR). In contrast, the frozen-prior Bayesian
rule selected 69.9%, 41.8%, and 31.6% of true mediators with no false
selections in this small smoke. The latter is encouraging but too small and
too close to development to support a release claim.

## Factorized multi-signal LD stress

`sim/run_factorized_ld_stress.py` generates genotype-derived PLINK references,
three independent cis blocks, LD-correlated GWAS errors, shared and distinct
causal signals, and an ancestry/reference mismatch. In the initial 12-run
smoke, both matched-reference cells classified shared signals as H4 and
distinct variants in LD as H3. The same-variant pleiotropy arm matched the
shared-mediation arm, as required by the identification theorem. The
completed replicated matrix contained 1,200 analyses. With matched LD,
the native regional model classified 85%-90% of distinct loci as H3 and 80%
of shared-mediation loci as H4; shared loci were incorrectly classified H3 in
1%-8% of runs. Same-variant pleiotropy tracked shared mediation, as required by
the nonidentifiability result. When target LD 0.9 was analyzed with a reference
having LD 0.6, Set B counts and two-leg evidence inflated sharply and false
evidence survived in some distinct-signal runs. An ancestry-matched LD panel
and explicit mismatch sensitivity analysis are therefore release requirements,
not optional refinements.

## Release interpretation

- The historical single-signal gate controls the tested distinct-LD failure well
  through moderate/high realized LD, but the 7.4% extreme-LD rate exceeds a 5%
  target. Those rates must not be attributed to the new multi-signal method.
- The native multi-signal method removes the one-causal-variant assumption but
  still depends on accurate ancestry-matched LD and conditional-signal recovery.
  Near-collinear loci should be reported as sensitivity-limited.
- Structural `P_M1` values are ELBO-based variational approximate model
  weights. The regional stress test does not make them exact posterior
  probabilities.
- These simulations validate implementation behavior against known truth;
  they do not calibrate real-data results from labels.

Accordingly, 1.2.0-dev is suitable only for research evaluation and transparent
GitHub development. It is not a production-calibrated release. The `-dev`
suffix should remain until the frozen factorized posterior rule passes
held-out calibration, LD-mismatch and interval-coverage criteria, followed by
real-data replication checks.

## Frozen Student-t/effect-validity validation

The current Student-t density-ratio e-value and the evidence/effect numerical
validity split were frozen before the new-seed Biowulf run at
`factorized_final_validation_20260903_v4` (seed `20261603`, binary SHA-256
`bd39c66663ca3d7b993ad4964f89344cff765be17d351df0153875c91cbbe193`).
All 3,200 tasks and 1.36 million protein analyses completed. The fixed e-BH
rule made no selections in any 200-replicate proteome family. Thus it was
conservative in this grid but had zero empirical power; the exact analytical
expectation argument, rather than these finite simulations, is the calibration
claim.

The default fixed-prior posterior-FDR rule selected only in the two broad-
instrument families. In the signed family it selected 11.74 proteins per
replicate on average, with mean FDR 0 and mean power 0.236. In the directional
family it selected 13.36, with mean FDR 0.0137, 95th-percentile replicate FDR
0.0914, and mean power 0.266. It selected nothing in the other families. These
are working-model posterior results, not an analytical frequentist guarantee.
In the prespecified sensitivity grid, lowering the directional-pleiotropy prior
from 0.10 to 0.05 increased aggregate power but raised the worst family mean
FDR to 0.0530; increasing it to 0.25 suppressed nearly all selection. The
default was not retuned after observing these results.

The new validity denominator exposed weak effect identification that previous
conditional-coverage summaries obscured. Two-leg evidence was finite for all
broad M1 rows, but stable joint slope/intercept effect estimates were available
for only 41.7%-42.7%. Median slope/intercept collinearity was 0.958-0.963 in the
broad cells and approximately 0.99 in the narrower baseline cells. Among the
subset with stable estimates, broad signed 95% interval coverage was 0.975 for
beta1, 0.991 for beta2, and 1.000 for the indirect effect. Those coverage values
must always be reported with the 41.7% availability denominator.

The companion 1,200-run genotype/LD matrix at
`factorized_ld_final_20260903_v5` retained finite two-leg evidence in every
run, eliminating the former 4%-14% hidden numerical-evidence failure rate.
Stable effect-estimate availability remained only 62%-82%. Matched-reference
distinct loci were classified H3 in 85%-90% of runs; shared mediation and
same-variant pleiotropy continued to produce similar H4 rates, as required by
the identification limit. Under target/reference LD mismatch, Set B counts
inflated from about 3 to 7.2-7.4 and false Bayesian support persisted.

These results support the evidence/effect separation and the H3/H4 diagnostic,
but they do not justify a production release. The analytical e-value is too
low-powered in the tested proteome families, effect estimation is often
unresolved under slope/intercept collinearity, and LD mismatch remains a major
failure mode. Version 1.2.0-dev therefore remains research-only.
