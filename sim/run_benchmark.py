from __future__ import annotations

import argparse
import csv
from pathlib import Path

from benchmark_lib import (
    attach_truth,
    benchmark_scenarios,
    calibration_sequence,
    classification_sequence,
    compute_benchmark_metrics,
    ensure_dir,
    generate_proteins,
    load_config,
    read_results,
    read_truth,
    run_bmediator,
    scenario_alias_map,
    summarize_calibration,
    summarize_classification,
    summarize_fdr_power,
    write_dataset,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run BMEDIATOR simulation benchmarks.")
    parser.add_argument("--config", type=Path, required=True, help="Benchmark config JSON.")
    parser.add_argument("--benchmark", choices=["classification", "calibration", "all"], default="all")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory for simulations.")
    parser.add_argument("--binary", type=Path, default=Path(__file__).resolve().parents[1] / "bmediator")
    parser.add_argument("--seed", type=int, default=None, help="Override global seed.")
    parser.add_argument("--dry-run", action="store_true", help="Generate inputs only, do not run bmediator.")
    return parser.parse_args()


def run_one_benchmark(name: str, cfg: dict, global_cfg: dict, outdir: Path, binary: Path, dry_run: bool, seed: int | None) -> None:
    import numpy as np

    benchmark_dir = ensure_dir(outdir / name)
    summary_dir = ensure_dir(outdir / "summary" / name)
    base_seed = int(global_cfg["seed"] if seed is None else seed)
    scenarios = tuple(str(s) for s in cfg.get("scenarios", benchmark_scenarios({"tmp": cfg}, "tmp")))
    alias_map = {str(k): str(v) for k, v in cfg.get("scenario_aliases", {}).items()}
    summary_scenarios = tuple(dict.fromkeys(alias_map.get(s, s) for s in scenarios))

    all_metrics = []
    all_class_rows = []
    all_confusion_rows = []
    all_calibration_rows = []
    all_bfdr_rows = []
    all_fdr_power_rows = []

    for cell_idx, cell in enumerate(cfg["cells"]):
        cell_name = cell["name"]
        replicates = int(cfg["replicates"])
        for rep in range(1, replicates + 1):
            rng = np.random.default_rng(base_seed + cell_idx * 10000 + rep)
            if name == "classification":
                scenario_sequence = classification_sequence(int(cfg["proteins_per_scenario"]), scenarios)
            else:
                scenario_sequence = calibration_sequence(
                    rng,
                    int(cfg["proteins_per_replicate"]),
                    cfg["scenario_mix"],
                    scenarios,
                )

            proteins = generate_proteins(rng, name, cell, scenario_sequence)
            files = write_dataset(outdir, name, cell_name, rep, proteins)
            out_prefix = files.truth.parent / "bmediator"

            if dry_run:
                continue

            mediation_path, _ = run_bmediator(binary, files, out_prefix, global_cfg)
            truth = read_truth(files.truth)
            results = read_results(mediation_path)
            merged = compute_benchmark_metrics(attach_truth(results, truth), alias_map=alias_map)
            for row in merged:
                row["benchmark"] = name
                row["cell"] = cell_name
                row["replicate"] = rep
            all_metrics.append(merged)
            fdr_power = summarize_fdr_power(merged)
            for row in fdr_power:
                row["benchmark"] = name
                row["cell"] = cell_name
                row["replicate"] = rep
            all_fdr_power_rows.append(fdr_power)

            if name == "classification":
                by_scenario, confusion = summarize_classification(merged, scenarios=summary_scenarios)
                for row in by_scenario:
                    row["cell"] = cell_name
                    row["replicate"] = rep
                for row in confusion:
                    row["cell"] = cell_name
                    row["replicate"] = rep
                all_class_rows.append(by_scenario)
                all_confusion_rows.append(confusion)
            else:
                calibration, bfdr = summarize_calibration(merged)
                for row in calibration:
                    row["cell"] = cell_name
                    row["replicate"] = rep
                for row in bfdr:
                    row["cell"] = cell_name
                    row["replicate"] = rep
                all_calibration_rows.append(calibration)
                all_bfdr_rows.append(bfdr)

    if dry_run:
        return

    if all_metrics:
        rows = [row for chunk in all_metrics for row in chunk]
        with (summary_dir / "protein_level_metrics.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    if all_class_rows:
        rows = [row for chunk in all_class_rows for row in chunk]
        with (summary_dir / "classification_by_scenario.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
    if all_confusion_rows:
        rows = [row for chunk in all_confusion_rows for row in chunk]
        with (summary_dir / "classification_confusion.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
    if all_calibration_rows:
        rows = [row for chunk in all_calibration_rows for row in chunk]
        with (summary_dir / "calibration_bins.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
    if all_bfdr_rows:
        rows = [row for chunk in all_bfdr_rows for row in chunk]
        with (summary_dir / "calibration_bfdr.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
    if all_fdr_power_rows:
        rows = [row for chunk in all_fdr_power_rows for row in chunk]
        with (summary_dir / "fdr_power_by_threshold.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "config.json").write_text(args.config.read_text())
    requested = ["classification", "calibration"] if args.benchmark == "all" else [args.benchmark]
    for benchmark in requested:
        run_one_benchmark(
            benchmark,
            config[benchmark],
            config["global"],
            args.outdir,
            args.binary,
            args.dry_run,
            args.seed,
        )


if __name__ == "__main__":
    main()
