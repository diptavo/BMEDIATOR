# JG-0.2.2 held-out family results

## Provenance and decision

The numerical-repair plan and implementation were frozen at `eb1f205`. The
Biowulf analysis ran as build job `29153958`, array job `29153959`, and summary
job `29153961` under
`/data/Dutta_lab/BMEDIATOR_runs/JG_0_2_2_eb1f205`. The source archive SHA-256
was `9b7fb20bd8db679fa00d561d71592244d01388d3f2f77379778e4d30ae2c5edf`.
All 50,000 analyses ran on batch nodes.

The FDR and power criteria passed among reportable fits, but the frozen
completion criterion failed in every confirmatory scenario:

| Scenario | Failed proteins | Mean FDR | Mean power |
| --- | ---: | ---: | ---: |
| baseline | 381/5000 | 0.0000 | 0.806 |
| rare | 371/5000 | 0.0000 | 0.690 |
| composite null | 194/5000 | 0.0000 | NA |
| mixed | 236/5000 | 0.0093 | 0.920 |
| strong LD | 394/5000 | 0.0000 | 0.810 |

Diagnostic failures ranged from 305 to 473 proteins per scenario. These were
controlled posterior suppressions, not process crashes. The result rejects
`JG-0.2.2` as a production engine because 3/5/7/9 tensor quadrature often did
not stabilize the sparse-plus-global-path states.

## Repair consequence

The sparse model mixed over block contamination with an unknown probability
`q` carrying a Beta prior. Numerically integrating transformed `q` was
unnecessary. `JG-0.2.3` computes that integral exactly by dynamic programming
and extends the remaining successive-order check through order 13. No
`JG-0.2.2` result is reused as final validation.
