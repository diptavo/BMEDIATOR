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

## Joint LD-aware implementation checks

The default `joint-ld` implementation is covered by a deterministic
genotype-based fixture containing two independent shared signals in one region
and a separate distinct-signal region. The regression suite requires:

- recovery of the expected conditional protein and outcome signals;
- shared diagonal and distinct off-diagonal H3/H4 signal-pair results;
- finite regional M0-M5 log Bayes factors and a well-conditioned regularized
  component LD matrix;
- rejection of the distinct-causal region by the confirmatory gate;
- explicit `UNRESOLVED_SINGLE_COMPONENT` reporting when only one component is
  retained;
- rejection of a non-positive-definite summary-error correlation matrix;
- invariance of H3/H4 output when a protein-outcome variant is absent from the
  RF summary statistics; and
- preservation of the `ld-multisignal` and `single` compatibility modes.

These are implementation and numerical regression checks. They are not an
operating-characteristic study of the new joint likelihood.

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

## Joint-likelihood example rerun

The BMI-CHD examples were rerun on Biowulf using the same curated inputs and
European GRCh38 LD reference as the earlier `ld-multisignal` analysis. The
shared-versus-distinct classifications were unchanged after keeping the
protein-outcome intersection separate from the three-trait joint intersection.
Structural M1 weights changed because the joint likelihood tests M1 directly
against M5 using the regional cross-trait effect pattern.

| Panel | Protein | Previous P(M1) | Joint-LD P(M1) | Joint-LD P(M5) | Regional state |
|---|---|---:|---:|---:|---|
| UKB-PPP | PCSK9 | 1.000000 | 0.016397 | 0.983476 | shared supported |
| UKB-PPP | IL6R | 0.999995 | 0.155580 | 0.333591 | ambiguous |
| UKB-PPP | NPPB | 0.156211 | 0.002667 | 0.996798 | ambiguous |
| deCODE | PCSK9 | 1.000000 | 0.007852 | 0.991731 | shared supported |
| deCODE | IL6R | 0.999999 | 0.176835 | 0.032336 | ambiguous |
| deCODE | NPPB | 0.012176 | 0.018272 | 0.967386 | ambiguous |

ANGPTL3 and VEGFA had no retained CHD signal in either panel and therefore did
not pass the identification gate. deCODE ANGPTL3 had `P(M1)=0.639449`, but its
confirmatory mediation score remained zero for this reason. No example protein
was selected at the default 5% or 10% mode-specific FDR. These examples show the
intended distinction between colocalization and mediation; they do not by
themselves establish calibration.

## Release interpretation

- The historical single-signal gate controls the tested distinct-LD failure well
  through moderate/high realized LD, but the 7.4% extreme-LD rate exceeds a 5%
  target. Those rates must not be attributed to the new multi-signal method.
- The native conditional multi-signal diagnostic removes the one-causal-variant
  assumption but still depends on accurate ancestry-matched LD and signal
  recovery.
  Near-collinear loci should be reported as sensitivity-limited.
- The joint regional likelihood now models signed LD, RF/protein/outcome effect
  patterns, configurable sampling-error correlation, and M5 correlated
  pleiotropy. Its scenario weights have not yet been calibrated in a broad
  simulation grid or independent empirical replication study.
- Structural `P_M1` values are ELBO-based variational approximate model
  weights. The regional stress test does not make them exact posterior
  probabilities.
- These simulations validate implementation behavior against known truth;
  they do not calibrate real-data results from labels.

Accordingly, 1.2.0-dev is suitable for research evaluation and transparent
GitHub distribution, but the `-dev` suffix should remain until the joint method
has a full simulation operating-characteristic study and sample-overlap,
ancestry-mismatch, prior-sensitivity, and real-data replication checks are
complete.
