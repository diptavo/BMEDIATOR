#!/usr/bin/env python3
"""Run one M1/M5 task with an independent dataset for every protein."""

from __future__ import annotations

import argparse
import csv
import copy
import shutil
from pathlib import Path

import numpy as np

from benchmark_lib import (
    attach_truth,
    classification_sequence,
    compute_benchmark_metrics,
    generate_proteins,
    load_config,
    read_results,
    read_truth,
    run_bmediator,
    write_dataset,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--cell", required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--structural-method", choices=("legacy-six-state", "factorized"))
    parser.add_argument("--keep-intermediate", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    global_cfg = copy.deepcopy(config["global"])
    benchmark_cfg = config["classification"]
    cell_index, cell = next(
        (index, value)
        for index, value in enumerate(benchmark_cfg["cells"])
        if value["name"] == args.cell
    )
    if args.structural_method:
        binary_options = dict(global_cfg.get("binary_options", {}))
        binary_options["structural_method"] = args.structural_method
        overlap = cell.get("sample_overlap", cell.get("sampling_error_correlation", {}))
        if isinstance(overlap, (int, float)):
            overlap = {"rf_pqtl": overlap, "rf_outcome": overlap, "pqtl_outcome": overlap}
        if isinstance(overlap, dict):
            binary_options["sampling_corr_rf_pqtl"] = overlap.get("rf_pqtl", 0.0)
            binary_options["sampling_corr_rf_outcome"] = overlap.get(
                "rf_outcome", overlap.get("rf_cancer", 0.0)
            )
            binary_options["sampling_corr_pqtl_outcome"] = overlap.get(
                "pqtl_outcome", overlap.get("pqtl_cancer", 0.0)
            )
        global_cfg["binary_options"] = binary_options
    base_seed = int(global_cfg["seed"] if args.seed is None else args.seed)
    rng = np.random.default_rng(base_seed + cell_index * 10000 + args.replicate)
    scenarios = tuple(str(value) for value in benchmark_cfg["scenarios"])
    sequence = classification_sequence(int(benchmark_cfg["proteins_per_scenario"]), scenarios)
    proteins = generate_proteins(rng, "classification", cell, sequence)

    rep_dir = args.outdir / "classification" / args.cell / f"rep_{args.replicate:04d}"
    work_root = rep_dir / "independent_runs"
    rep_dir.mkdir(parents=True, exist_ok=True)
    attached = []

    for index, protein in enumerate(proteins, start=1):
        files = write_dataset(work_root, "protein", f"p{index:05d}", 1, [protein])
        out_prefix = files.truth.parent / "bmediator"
        mediation_path, _ = run_bmediator(args.binary, files, out_prefix, global_cfg)
        result = read_results(mediation_path)
        truth = read_truth(files.truth)
        merged = attach_truth(result, truth)
        if len(merged) != 1:
            raise RuntimeError(f"Expected one result for {protein['protein_id']}, got {len(merged)}")
        attached.extend(merged)

    metrics = compute_benchmark_metrics(attached)
    for row in metrics:
        row["benchmark"] = "classification"
        row["cell"] = args.cell
        row["replicate"] = args.replicate

    task_metrics = rep_dir / "task_metrics.tsv"
    with task_metrics.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(metrics)

    if not args.keep_intermediate:
        shutil.rmtree(work_root)
    print(f"Wrote {task_metrics} ({len(metrics)} independent proteins)")


if __name__ == "__main__":
    main()
