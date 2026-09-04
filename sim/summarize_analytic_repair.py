#!/usr/bin/env python3
"""Summarize replicate-level FDR and power for analytical repair candidates."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


METHODS = {
    "directional_BY": "factor_conjunction_q_BY",
    "balanced_BY": "factor_balanced_conjunction_q_BY",
    "balanced_AdaFilter": "factor_balanced_conjunction_q_AdaFilter",
    "balanced_eBH": "factor_e_q_balanced_EBH",
    "balanced_p2e_eBH": "factor_e_q_p2e_balanced_EBH",
    "directional_eBH": "factor_e_q_EBH",
    "adaptive_eBH": "factor_e_q_adaptive_EBH",
    "bayesian_FDR": "factor_posterior_cum_fdr",
}


def number(value: str) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else math.nan
    except (TypeError, ValueError):
        return math.nan


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def sample_se(values: list[float]) -> float:
    if len(values) < 2:
        return math.nan
    center = mean(values)
    variance = sum((value - center) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(variance / len(values))


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = probability * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thresholds", default="0.01,0.05,0.10")
    args = parser.parse_args()
    thresholds = [float(value) for value in args.thresholds.split(",")]

    grouped: dict[tuple[str, str, float], list[dict[str, float]]] = defaultdict(list)
    by_truth: dict[tuple[str, str, float, str], list[dict[str, float]]] = defaultdict(list)
    paths = sorted(args.input.glob("calibration/*/rep_*/task_metrics.tsv"))
    if not paths:
        raise SystemExit(f"No calibration task metrics under {args.input}")

    for path in paths:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        cell = rows[0]["cell"]
        total_true = sum(row["true_scenario"] == "M1" for row in rows)
        for method, column in METHODS.items():
            for threshold in thresholds:
                selected = [row for row in rows if number(row.get(column, "")) <= threshold]
                true_positive = sum(row["true_scenario"] == "M1" for row in selected)
                false_positive = len(selected) - true_positive
                grouped[(cell, method, threshold)].append({
                    "selected": len(selected),
                    "tp": true_positive,
                    "fp": false_positive,
                    "fdp": false_positive / max(1, len(selected)),
                    "power": true_positive / total_true if total_true else 0.0,
                })
                for scenario in ("M0", "M1", "M2", "M3", "M4", "M5"):
                    scenario_rows = [row for row in rows if row["true_scenario"] == scenario]
                    scenario_selected = sum(
                        number(row.get(column, "")) <= threshold for row in scenario_rows
                    )
                    by_truth[(cell, method, threshold, scenario)].append({
                        "selected": scenario_selected,
                        "total": len(scenario_rows),
                        "rate": scenario_selected / len(scenario_rows) if scenario_rows else 0.0,
                    })

    fields = [
        "cell", "method", "threshold", "replicates", "mean_selected",
        "mean_true_positive", "mean_false_positive", "mean_fdp",
        "se_mean_fdp", "q95_fdp", "fdp_above_target_rate", "pooled_fdp",
        "mean_power", "replicates_with_calls",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (cell, method, threshold), records in sorted(grouped.items()):
            selected = sum(record["selected"] for record in records)
            tp = sum(record["tp"] for record in records)
            fp = sum(record["fp"] for record in records)
            writer.writerow({
                "cell": cell,
                "method": method,
                "threshold": threshold,
                "replicates": len(records),
                "mean_selected": f"{mean([r['selected'] for r in records]):.6g}",
                "mean_true_positive": f"{mean([r['tp'] for r in records]):.6g}",
                "mean_false_positive": f"{mean([r['fp'] for r in records]):.6g}",
                "mean_fdp": f"{mean([r['fdp'] for r in records]):.6g}",
                "se_mean_fdp": f"{sample_se([r['fdp'] for r in records]):.6g}",
                "q95_fdp": f"{quantile([r['fdp'] for r in records], 0.95):.6g}",
                "fdp_above_target_rate": f"{mean([r['fdp'] > threshold for r in records]):.6g}",
                "pooled_fdp": f"{fp / max(1, selected):.6g}",
                "mean_power": f"{mean([r['power'] for r in records]):.6g}",
                "replicates_with_calls": sum(r["selected"] > 0 for r in records),
            })
    truth_output = args.output.with_name(args.output.stem + "_by_truth.tsv")
    truth_fields = [
        "cell", "method", "threshold", "true_scenario", "replicates",
        "total_proteins", "selected_proteins", "mean_replicate_selection_rate",
        "pooled_selection_rate",
    ]
    with truth_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=truth_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for (cell, method, threshold, scenario), records in sorted(by_truth.items()):
            total = sum(record["total"] for record in records)
            selected = sum(record["selected"] for record in records)
            writer.writerow({
                "cell": cell,
                "method": method,
                "threshold": threshold,
                "true_scenario": scenario,
                "replicates": len(records),
                "total_proteins": total,
                "selected_proteins": selected,
                "mean_replicate_selection_rate": f"{mean([r['rate'] for r in records]):.6g}",
                "pooled_selection_rate": f"{selected / total if total else math.nan:.6g}",
            })
    print(f"Wrote {args.output} and {truth_output} from {len(paths)} replicate tasks")


if __name__ == "__main__":
    main()
