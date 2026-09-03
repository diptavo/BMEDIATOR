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

Accordingly, 1.2.0-dev is suitable for research evaluation and transparent
GitHub distribution, but the `-dev` suffix should remain until the new
multi-signal method has a full simulation operating-characteristic study and
sample-overlap, ancestry-mismatch, and real-data replication checks are
complete.
