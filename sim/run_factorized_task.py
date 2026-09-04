#!/usr/bin/env python3
"""Run one batched factorized calibration replicate.

Every protein in a replicate is analyzed in one BMEDIATOR invocation so the
reported Benjamini-Yekutieli q-values are calculated over the simulated
proteome, not separately for each protein.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from benchmark_lib import (
    BenchmarkFiles,
    attach_truth,
    benchmark_scenarios,
    calibration_sequence,
    classification_sequence,
    generate_proteins,
    load_config,
    read_results,
    read_truth,
    run_bmediator,
    sampling_error_correlations,
    write_dataset,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--benchmark", choices=("classification", "calibration"), required=True)
    parser.add_argument("--cell", required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--keep-intermediate", action="store_true")
    parser.add_argument("--reuse-output", action="store_true")
    return parser.parse_args()


def factorized_options(global_cfg: dict, cell: dict) -> dict:
    configured = global_cfg.get("binary_options", {})
    if isinstance(configured, list):
        raise ValueError("factorized calibration requires binary_options to be a JSON object")
    options = dict(configured)
    options["structural_method"] = "factorized"
    rf_pqtl, rf_outcome, pqtl_outcome = sampling_error_correlations(cell)
    options["sampling_corr_rf_pqtl"] = rf_pqtl
    options["sampling_corr_rf_outcome"] = rf_outcome
    options["sampling_corr_pqtl_outcome"] = pqtl_outcome
    if bool(cell.get("independent_selection", False)):
        options["factor_independent_selection"] = True
    cell_options = cell.get("binary_options", {})
    if not isinstance(cell_options, dict):
        raise ValueError("cell binary_options must be a JSON object")
    options.update(cell_options)
    return options


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    bench_cfg = config[args.benchmark]
    cell_index, configured_cell = next(
        (index, value)
        for index, value in enumerate(bench_cfg["cells"])
        if value["name"] == args.cell
    )
    cell = copy.deepcopy(configured_cell)
    global_cfg = copy.deepcopy(config["global"])
    configured_options = global_cfg.get("binary_options", {})
    if isinstance(configured_options, dict) and bool(
        configured_options.get("factor_independent_selection", False)
    ):
        cell["independent_selection"] = True
    global_cfg["binary_options"] = factorized_options(global_cfg, cell)

    base_seed = int(global_cfg["seed"] if args.seed is None else args.seed)
    task_seed = base_seed + cell_index * 100_000 + args.replicate
    rng = np.random.default_rng(task_seed)
    scenarios = benchmark_scenarios(config, args.benchmark)
    if args.benchmark == "classification":
        sequence = classification_sequence(int(bench_cfg["proteins_per_scenario"]), scenarios)
    else:
        sequence = calibration_sequence(
            rng,
            int(bench_cfg["proteins_per_replicate"]),
            cell.get("scenario_mix", bench_cfg["scenario_mix"]),
            scenarios,
        )

    rep_dir = args.outdir / args.benchmark / args.cell / f"rep_{args.replicate:04d}"
    if args.reuse_output:
        files = BenchmarkFiles(
            rf=rep_dir / "rf_sumstat.txt",
            pqtl=rep_dir / "pqtl_sumstat.txt",
            cancer=rep_dir / "cancer_sumstat.txt",
            protein_info=rep_dir / "protein_info.txt",
            truth=rep_dir / "truth.tsv",
        )
        proteins = None
    else:
        proteins = generate_proteins(rng, args.benchmark, cell, sequence)
        files = write_dataset(args.outdir, args.benchmark, args.cell, args.replicate, proteins)
    out_prefix = rep_dir / "bmediator"
    mediation_path = out_prefix.with_suffix(".mediation")
    if not args.reuse_output:
        mediation_path, _ = run_bmediator(args.binary, files, out_prefix, global_cfg)
    elif not mediation_path.is_file() or not files.truth.is_file():
        raise RuntimeError("--reuse-output requires existing bmediator.mediation and truth.tsv")
    truth = read_truth(files.truth)
    rows = attach_truth(read_results(mediation_path), truth)
    if len(rows) != len(truth):
        raise RuntimeError(f"expected {len(truth)} result rows, found {len(rows)}")

    for row in rows:
        row["benchmark"] = args.benchmark
        row["cell"] = args.cell
        row["replicate"] = args.replicate
        row["true_mediation"] = int(str(row["true_scenario"]) == "M1")
        row["conjunction_null"] = 1 - row["true_mediation"]

    task_metrics = rep_dir / "task_metrics.tsv"
    with task_metrics.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "benchmark": args.benchmark,
        "cell": args.cell,
        "replicate": args.replicate,
        "seed": task_seed,
        "proteins": len(rows),
        "binary_sha256": sha256_file(args.binary.resolve()),
        "config_sha256": sha256_file(args.config.resolve()),
        "factorized_options": global_cfg["binary_options"],
    }
    (rep_dir / "task_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )

    if not args.keep_intermediate:
        for path in (files.rf, files.pqtl, files.cancer, files.protein_info, files.truth):
            if path.exists():
                path.unlink()
        for path in rep_dir.glob("bmediator.*"):
            path.unlink()
    print(f"Wrote {task_metrics} ({len(rows)} jointly analyzed proteins)")


if __name__ == "__main__":
    main()
