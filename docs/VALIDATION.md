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
PCSK9 results that shifted support from M1 to M5 must not be used. The
production model remains the Set A/B/C CAVI model plus the separate
`ld-multisignal` H3/H4 regional diagnostic.

## Native multi-signal checks

The default `ld-multisignal` implementation is covered by a deterministic
genotype-based fixture containing two independent shared signals in one region
and a separate distinct-signal region. The test requires recovery of two
signals per trait, shared diagonal signal pairs, distinct off-diagonal pairs,
and rejection of the distinct-causal region.

The production implementation has also been compared with a SuSiE/coloc
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

On 2026-09-03, the production Set A/B/C structural model was evaluated in
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

This pilot establishes that separating Set A and Set B repairs the specific
M1/M5 conflation in these independent-instrument simulations. It does not
validate the complete method. Realistic LD, sample overlap, winner's curse,
directional pleiotropy, ancestry mismatch, regional H3/H4 error, and
proteome-wide FDR remain open validation requirements. The factorized mode is
therefore experimental and the development release warning still applies.

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
suffix should remain until the structural likelihood is revised and passes
M1/M5 discrimination and probability-calibration criteria, followed by
ancestry-mismatch and real-data replication checks.
