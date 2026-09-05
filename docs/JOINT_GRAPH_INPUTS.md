# Joint-model input contract

`bmediator-joint` runs the `JG-0.2.8` joint Bayesian model for one protein. It
does not select variants or estimate nuisance scales from the same analyzed
associations. Inputs must be harmonized and prepared before this command is
used.

## Per-protein table

The input TSV has one row per retained variant and exactly these required
columns:

```text
variant ld_block role beta_x se_x beta_m se_m beta_y se_y v_x v_m v_y orientation orientation_probability rho_xm rho_xy rho_my
```

- `beta_x`, `beta_m`, and `beta_y` are effect-allele-aligned associations for
  the risk factor, protein, and outcome.
- `role` is `A`, `B`, or `C`. A is RF-dominant evidence, B is protein-dominant
  evidence, and C is overlapping evidence. At least three independent A-role
  and three independent B-role blocks are required by default.
- `ld_block` identifies approximately independent blocks. Absolute signed LD
  between different blocks must not exceed `0.05` unless a prespecified option
  lowers that threshold.
- `v_x`, `v_m`, and `v_y` are nonnegative role-specific latent variances from
  an independent, LD-pruned estimation panel. They are not estimated from the
  selected protein analysis.
- `orientation` is `-1` or `1`, fixed in independent discovery data for the
  directional-pleiotropy component. `orientation_probability` is its
  externally justified probability of being correct and lies in `[0.5,1]`.
- `rho_xm`, `rho_xy`, and `rho_my` are declared cross-GWAS sampling-error
  correlations. They are constant within one protein and must define a
  positive-definite correlation matrix.

Rows with missing or nonfinite values are rejected. Alleles, genome build,
variant order, LD signs, and effect directions must already agree.

After harmonization, independently estimated role scales can be joined with:

```bash
Rscript research/prepare_joint_graph_input.R \
    protein.harmonized.tsv external_scales.tsv protein.joint.tsv
```

This validates the schema and independent-block requirement. It does not infer
allele alignment, LD blocks, sample overlap, or orientation probabilities.

## LD table

The LD TSV is a square signed-correlation matrix. Its first row and first
column contain variant identifiers; these may be in a different order from the
analysis table because the executable matches them by identifier. Every input
variant must occur exactly once.

## One protein

```bash
./bmediator-joint \
    --input protein.joint.tsv \
    --ld protein.ld.tsv \
    --out protein.result.tsv
```

An optional tab-delimited key/value file can override prespecified priors and
make numerical safeguards stricter:

```bash
./bmediator-joint \
    --input protein.joint.tsv \
    --ld protein.ld.tsv \
    --out protein.result.tsv \
    --options joint_options.tsv
```

Every result records the model version, identification scope, role and block
counts, all priors, numerical settings, posterior state probabilities, and
quadrature diagnostics. A nonzero exit status means no posterior is
reportable. Options cannot raise the cross-block LD or numerical-error limits,
lower the independent-block requirement, or weaken optimizer settings below
the release defaults.

## Protein manifest

The manifest runner accepts a TSV with `protein`, `input`, `ld`, and optional
`options` columns. Relative paths are resolved from the manifest directory.

```bash
Rscript research/run_joint_graph_manifest.R \
    ./bmediator-joint protein_manifest.tsv results/bmediator_joint 8
```

It writes `results/bmediator_joint.joint.tsv` and
`results/bmediator_joint.failures.tsv`. The first file preserves manifest
order and adds `posterior_lfdr`, `posterior_rank`,
`posterior_cumulative_fdr`, `family_complete`, `posterior_fdr_status`, and
`selected_bfdr05`. The runner returns status 2 and leaves
the family-wide rank, cumulative FDR, and `selected_bfdr05` as `NA` when any
manifest row fails. Repair and rerun those rows before interpreting family-wide
selection. The selection rule controls
posterior expected FDR under the fitted model; it does not remove violations
of the model assumptions.

## Interpretation boundary

`PP_two_path` is support for both `X -> M` and global `M -> Y` paths within the
16-state model. It is always conditional on no exactly aligned pleiotropic
effect that is observationally identical to `M -> Y`. Such data are not
identifiable from these summary statistics, regardless of software or sample
size.
