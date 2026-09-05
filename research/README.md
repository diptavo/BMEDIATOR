# Joint graph research prototype

The files in this directory implement and test the historical JG-0.1
reference model and the current JG-0.2 series. The supported research
executable is `bmediator-joint`; see `docs/JOINT_GRAPH_INPUTS.md`.

Build the standalone C++ evaluator and run the cross-language fixtures with:

```bash
make test-joint-graph
```

Run the prespecified 200-dataset development sweep with:

```bash
Rscript research/evaluate_joint_graph_mechanisms.R
```

Run the frozen step-5 adversarial matrix with:

```bash
Rscript research/evaluate_joint_graph_adversarial.R
```

Its scenarios and interpretation rules are fixed in
`docs/JOINT_GRAPH_V0_1_ADVERSARIAL_PLAN.md`.

The frozen full-run outputs are:

```text
joint_graph_v0_1_adversarial_summary.tsv
joint_graph_v0_1_adversarial_summary_replicates.tsv
```

See `docs/JOINT_GRAPH_V0_1_ADVERSARIAL_RESULTS.md` for interpretation. The
evaluation rejects JG-0.1 as a production engine; the files are retained to
prevent the development failures from being hidden by later model changes.

## JG-0.2 repair

Build and test the adaptive signed-LD implementation with:

```bash
make test-joint-graph-v02
```

Estimate role variances from an independent, LD-pruned estimation panel with:

```bash
Rscript research/estimate_joint_graph_v02_scales.R external.tsv scales.tsv
```

Run the frozen new-seed development matrix with:

```bash
Rscript research/evaluate_joint_graph_v02.R
```

The model specification is in `docs/JOINT_GRAPH_MODEL_V0_2.md`. The final
frozen JG-0.2.4 family run passed its FDR and power criteria but failed its
numerical-completion gate; see
`docs/JOINT_GRAPH_V0_2_4_HELDOUT_RESULTS.md`. JG-0.2.4 is not a production
engine.

`JG-0.2.7` replaces the JG-0.2.6 high-order tensor fallback with deterministic
Smolyak sparse-grid integration for difficult posteriors. It retains the same
joint model and fixed reporting gates. Cluster replay and new-seed development
validation are required before a frozen confirmatory run.

The evaluator accepts a tab-delimited file with columns:

```text
variant role beta_x se_x beta_m se_m beta_y se_y
```

Run the implementations directly with:

```bash
Rscript research/joint_graph_reference.R input.tsv output.r.tsv
build/joint_graph_cli input.tsv output.cpp.tsv
```

`joint_graph_v0_1_development_sweep.tsv` and its replicate file preserve the
development behavior check. These are not publication-calibration results.
