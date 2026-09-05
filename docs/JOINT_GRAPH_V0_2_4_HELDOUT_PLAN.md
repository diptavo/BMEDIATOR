# JG-0.2.4 final held-out family plan

## Freeze and numerical rule

This plan is frozen before any `JG-0.2.4` family result is generated. The
`JG-0.2.3` family run failed its completion criterion and is development
evidence. Seed base `50400000` is disjoint from every earlier run.

The structural model, priors, exact Beta integration over `q`, quadrature
orders, and posterior selection rule are unchanged. Only the numerical
reportability rule changes. Treat each state's last successive log-evidence
difference as a symmetric error radius and enumerate all `2^16` corners of
that error box. A posterior is reportable only when the largest resulting
total-variation change in the normalized 16-state posterior is at most `0.01`
and every relevant state's difference is at most `1.0`.

## Run and acceptance

The run contains 10 unchanged scenarios, 50 families per scenario, and 100
proteins per family: 50,000 analyses. The primary 5% posterior-FDR rule and
secondary `PP_two_path >= 0.80` rule are unchanged.

Frozen confirmatory criteria are unchanged:

- `baseline`: mean FDR at most 0.05 and mean power at least 0.70;
- `rare`: mean FDR at most 0.05 and mean power at least 0.50;
- `composite_null`: at most 0.05 false discoveries per family and at most 5%
  of families with any discovery;
- `mixed`: mean FDR at most 0.05 and mean power at least 0.60;
- `strong_ld`: mean FDR at most 0.05 and mean power at least 0.65;
- at least 49 of 50 fully completed families in each confirmatory scenario.

Diagnostic scenarios and the exact-alignment identification boundary retain
their prior status. Criteria will not be changed after results are read. All
Biowulf computation must use batch or allocated interactive nodes.
