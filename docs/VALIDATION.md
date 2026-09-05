# Validation Status

## Current joint-model status

The current research implementation is the JG-0.2.6 joint model in
`bmediator-joint`. The earlier JG-0.2.4 frozen 50,000-analysis family run passed all
prespecified FDR and power criteria among reportable fits but failed the
numerical-completion gate, with only 38 to 48 of 50 complete families in the
confirmatory scenarios. It is not a production engine. See
[JG-0.2.4 held-out results](JOINT_GRAPH_V0_2_4_HELDOUT_RESULTS.md),
[JG-0.2.6 development results](JOINT_GRAPH_V0_2_6_DEVELOPMENT_RESULTS.md), and
[Production readiness](PRODUCTION_READINESS.md).

The sections below retain historical validation results for the legacy and
factorized implementations. They must not be interpreted as validation of the
joint model.

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

## Analytical-repair LD validation

The frozen 2026-09-04 genotype-based matrix comprised 1,200 analyses: 100
replicates of four truth scenarios under matched high LD, matched moderate LD,
and high-to-moderate reference mismatch. All array tasks and the compute-node
summary completed successfully as Slurm jobs `29061808` and `29061811`.

Under the no-second-stage null, the balanced two-leg nominal rejection rate was
0.04 with matched high LD and 0.06 with matched moderate LD, satisfying the
prespecified upper criterion of 0.075. It was 0.07 under deliberate LD mismatch.
For distinct causal variants in LD, the balanced two-leg test was often positive
(0.87 in matched high LD and 0.46 in matched moderate LD), as expected because
both causal legs can exist without mediation. The regional model classified
these loci as distinct in 0.90 and 0.85 of runs, respectively, and no run
received final factorized mediation support.

Same-variant pleiotropy and shared mediation remained observationally similar:
both frequently had balanced two-leg evidence and shared regional signals. The
final method must therefore continue to label support as exclusion-restriction
conditional rather than claiming that LD alone distinguishes these mechanisms.
The compact frozen table is
`sim/results/analytic_repair_ld_20260904_v1/factorized_ld_stress_summary.tsv`.

## Analytical-repair proteome calibration

The 2026-09-04 confirmatory grid used 200 replicates of 500 jointly analyzed
proteins in each of ten prespecified cells, for 1,000,000 protein analyses.
Arrays `29061505` and `29061533` and compute-node summary `29061536`
completed without failed tasks. The balanced AdaFilter rule retained mean
power of 0.973 in the broad cell and 0.675 in the narrow cell, with mean FDP
of 0.0379 and 0.0399. It nevertheless selected at least one false mediator in
12 of 200 least-favorable all-null families, giving FDR 0.060. This fails the
frozen requirement of no more than 0.05 and prevents designation of AdaFilter
as the confirmatory primary method. The estimate is compatible with Monte
Carlo variation around 0.05, but the endpoint was not passed.

The balanced p-to-e/e-BH candidate was then run from a separately built and
tested source snapshot as arrays `29066362` and `29066364`, with compute-node
summary `29066368`. It made no false selections in the broad, narrow,
least-favorable-null, global-null, heterogeneity, sparse-outlier, dense-
dependence, or directional-pleiotropy cells. It therefore passed the frozen
FDR criteria, but mean power was only 0.0969 in the broad cell and zero in the
narrow cell. It is analytically defensible but not an adequately powered
replacement.

The direct balanced Student-t/e-BH candidate was frozen before its results and
run from a separately built and tested source snapshot as arrays `29075015`
and `29075017`, with compute-node summary `29075025`. All 2,000 array tasks
completed successfully, producing 1,000,000 protein analyses. At 0.05, mean
FDP was 0.0021 in the broad cell and zero in every other identifiable valid or
stress cell; the global and least-favorable all-null cells had no selections.
It therefore passed the frozen FDR checks in this grid. Mean power was 0.6382
in the broad cell and zero in the narrow cell, failing the prespecified 0.80
and 0.50 targets. The rule is materially more powerful than the balanced
p-to-e transform in broad architectures, but it is still not an adequately
powered primary method. At the exactly proportional M5 boundary it selected
in 15 of 200 families; among selected proteins the pooled false fraction was
0.476, consistent with that boundary being nonidentifiable rather than
analytically calibratable.

The arbitrary-dependence balanced BY comparator had broad-cell power 0.792 and
narrow-cell power 0.0038. Directional pleiotropy and exactly proportional M5
effects remained explicit assumption and identification boundaries. Compact
tables are retained under `sim/results/analytic_repair_main_20260904_v1/` and
`sim/results/analytic_repair_p2e_main_20260904_v1/`. The direct Student-t/e-BH
tables are under
`sim/results/analytic_repair_balanced_e_main_20260904_v1/`.

Ordinary BH on the balanced partial-conjunction p-value was subsequently
specified as a new candidate before examining any BH results. Unlike
AdaFilter, it does not filter on the companion leg p-value and does not require
independence between causal legs. Its FDR theorem requires independence or
PRDS across protein-level p-values. The frozen new-seed design and automated
decision criteria are recorded in `sim/configs/analytic_bh_validation.json`,
`sim/ANALYTIC_REPAIR_DECISION_RULES.md`, and
`sim/evaluate_analytic_bh_validation.py`.

That validation was run from frozen commit `dc6eb00` as four 900-task Slurm
arrays (`29082551`, `29082559`, `29082564`, and `29082573`) with compute-node
aggregation in job `29082597`. All 3,600 tasks completed, comprising 1.8
million protein analyses. At `alpha=0.05`, balanced BH passed every frozen
criterion:

| Cell | Mean FDP or null-family FDR | Mean power | Frozen result |
|---|---:|---:|---|
| Broad balanced | 0.0376 | 0.9720 | PASS |
| Narrow balanced | 0.0298 | 0.5945 | PASS |
| Rare mediation | 0.0029 | 0.1039 | PASS; power descriptive |
| Least-favorable null | 0.0090 | n/a | PASS |
| Global null | 0.0000 | n/a | PASS |
| Balanced heterogeneity | 0.0043 | 0.0087 | PASS stress criterion |
| Sparse outliers | 0.0082 | 0.0647 | PASS stress criterion |
| Dense cross-protein dependence | 0.0232 | 0.5574 | PASS descriptive stress criterion |

The least-favorable and global-null one-sided 95% Wilson upper bounds were
0.0154 and 0.0027, below the prespecified 0.065 limit. The dense-dependence
result is empirical stress evidence, not proof of PRDS; generic correlation of
two-sided p-values is insufficient for the ordinary-BH theorem. Directional
pleiotropy produced mean FDP 0.5753, as expected outside the balanced/InSIDE
model. At the exactly proportional M5 boundary, M1 and M5 were selected at
nearly equal rates and the mean FDP was 0.5008, confirming rather than
resolving the stated nonidentifiability. Compact result tables are retained
under `sim/results/analytic_bh_validation_20260905_v1/summary/`.

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

Accordingly, balanced partial-conjunction BH is the first candidate to pass
the prespecified held-out FDR and power endpoints. This supports its use as the
leading confirmatory statistical rule when the balanced/InSIDE leg model and
cross-protein independence/PRDS assumptions are defensible. Version 1.2.0-dev
nevertheless remains suitable only for research evaluation and transparent
GitHub development. The `-dev` suffix should remain until LD-mismatch,
sample-overlap, effect-bias and interval-coverage criteria are completed,
published competitor implementations are benchmarked, and real-data
replication checks succeed. BY and e-BH remain the conservative
arbitrary-dependence sensitivity analyses, although their frozen power was
inadequate.

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
