# JG-0.2.1 frozen patch-validation plan

This patch plan is frozen after the original `JG-0.2` results and before the
full patch run. It uses seeds beginning at `20286001`; none were used for
`JG-0.1` or `JG-0.2` evaluation. There are 50 replicates per cell.

The posterior patterns and 80% minimum correct-pattern rate remain unchanged
from `JOINT_GRAPH_V0_2_VALIDATION_PLAN.md`. The matrix also requires every
confirmatory replicate to complete all 16 state optimizations without Hessian
regularization and with relevant-state adaptive/Laplace discrepancy at most
one log unit.

The original 20 cells are repeated, with three prespecified changes:

1. uncertain orientations are supplied with their independent correctness
   probabilities and exactly integrated;
2. `orientation_70pct` is now a confirmatory directional-pleiotropy cell;
3. the failed 20-block mediation-plus-sparse cell remains diagnostic and is
   paired with a 30-block confirmatory cell using the same latent seed.

Skeptical-prior mediation and diffuse-prior null remain paired to their default
counterparts. Exact and near-aligned pleiotropy, undeclared sample overlap,
misspecified scales, weak mediation, and six-block sparse pleiotropy remain
diagnostic because they violate assumptions or deliberately reduce
information.

Passing this patch matrix would establish repair of the observed development
failures only. Production still requires a substantially larger held-out
calibration and competitor benchmark on batch or interactive compute nodes.
