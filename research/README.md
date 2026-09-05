# Joint graph research prototype

The files in this directory implement and test the frozen JG-0.1 reference
model described in `docs/JOINT_GRAPH_MODEL_V0_1.md`.

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
