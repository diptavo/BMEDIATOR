# JG-0.2.3 final held-out family plan

## Freeze

This plan is frozen before any `JG-0.2.3` family result is generated. The
`JG-0.2.2` family run failed its completion gate and is development evidence.
No seed from that run is reused.

`JG-0.2.3` preserves the same likelihood and priors but evaluates the Beta
prior integral over block-contamination probability `q` exactly. Relevant
continuous-state integrals use successive mode-centered Gauss-Hermite orders
3, 5, 7, 9, 11, and 13. The final successive-order difference must be at most
`0.10` log units. At least three independent A-role and B-role blocks remain
mandatory.

## Independent run and criteria

Seed base `40400000` is disjoint from all previous development and family
runs. The run contains the unchanged 10 scenarios, 50 families per scenario,
and 100 proteins per family, for 50,000 analyses. The primary 5% posterior-FDR
rule and secondary `PP_two_path >= 0.80` rule are unchanged.

Frozen confirmatory criteria are unchanged:

- `baseline`: mean FDR at most 0.05 and mean power at least 0.70;
- `rare`: mean FDR at most 0.05 and mean power at least 0.50;
- `composite_null`: at most 0.05 false discoveries per family and at most 5%
  of families with any discovery;
- `mixed`: mean FDR at most 0.05 and mean power at least 0.60;
- `strong_ld`: mean FDR at most 0.05 and mean power at least 0.65;
- at least 49 of 50 fully completed families in each confirmatory scenario.

The other scenarios and exact-alignment boundary retain their diagnostic
status. Criteria will not be changed after results are observed. All
computation must run on Biowulf batch or allocated interactive nodes.
