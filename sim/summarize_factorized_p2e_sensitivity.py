#!/usr/bin/env python3
"""Evaluate analytically valid p-to-e calibrators on completed factorized runs.

The calculation uses only the valid conjunction p-value. It does not fit a
null distribution or tune a threshold from simulated labels. Simulation truth
is used solely to report operating characteristics after the decision rule is
applied.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from summarize_factorized_calibration import as_float, mean, quantile


CALIBRATORS = {
    "kappa_0.10": (0.10,),
    "kappa_0.25": (0.25,),
    "kappa_0.50": (0.50,),
    "kappa_0.75": (0.75,),
    "mixture_0.10_0.25_0.50_0.75": (0.10, 0.25, 0.50, 0.75),
}


def log_sum_exp(values):
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def log_e_from_p(p_value, kappas):
    if not math.isfinite(p_value) or p_value < 0.0 or p_value > 1.0:
        return math.nan
    p_value = max(p_value, float.fromhex("0x1.0p-1022"))
    terms = [math.log(kappa) + (kappa - 1.0) * math.log(p_value) for kappa in kappas]
    return log_sum_exp(terms) - math.log(len(kappas))


def adjusted_ebh(log_e_values):
    order = sorted(range(len(log_e_values)), key=lambda idx: log_e_values[idx], reverse=True)
    count = len(order)
    adjusted = [math.nan] * count
    running_log_q = math.inf
    for rank in range(count, 0, -1):
        index = order[rank - 1]
        raw_log_q = math.log(count) - math.log(rank) - log_e_values[index]
        running_log_q = min(running_log_q, raw_log_q)
        adjusted[index] = 1.0 if running_log_q >= 0.0 else math.exp(running_log_q)
    return adjusted


def task_metrics(rows, kappas, alpha):
    # BY values are deliberately unavailable when causal-leg GWAS overlap is
    # declared. Use that as an output-level marker that conjunction p-values
    # are not eligible for confirmatory p-to-e conversion in this task.
    if any(row.get("factor_frequentist_status") == "UNRESOLVED_SAMPLE_OVERLAP" for row in rows):
        return None

    tested = []
    for row in rows:
        p_value = as_float(row.get("factor_conjunction_p"))
        if not math.isfinite(p_value):
            continue
        log_e = log_e_from_p(p_value, kappas)
        item = {
            "log_e": log_e,
            "truth": row.get("true_scenario") == "M1",
            "scenario": row.get("true_scenario", ""),
            "boundary": row.get("identification_class") == "nonidentifiable",
        }
        tested.append(item)

    q_values = adjusted_ebh([item["log_e"] for item in tested])
    selected = [item for item, q_value in zip(tested, q_values) if q_value <= alpha]
    identifiable_selected = [item for item in selected if not item["boundary"]]
    true_total = sum(item["truth"] for item in tested)
    true_selected = sum(item["truth"] for item in selected)
    false_selected = len(selected) - true_selected
    identifiable_false = sum(not item["truth"] for item in identifiable_selected)
    boundary = [item for item in tested if item["boundary"]]
    boundary_selected = sum(item["boundary"] for item in selected)

    null_log_e = [
        item["log_e"] for item in tested
        if not item["truth"] and not item["boundary"]
    ]
    maximum = max(null_log_e) if null_log_e else math.nan
    log_mean_null_e = (
        maximum + math.log(mean([math.exp(value - maximum) for value in null_log_e]))
        if null_log_e else math.nan
    )
    return {
        "n_tested": len(tested),
        "n_selected": len(selected),
        "causal_fdr": false_selected / len(selected) if selected else 0.0,
        "identifiable_fdr": (
            identifiable_false / len(identifiable_selected)
            if identifiable_selected else 0.0
        ),
        "power": true_selected / true_total if true_total else math.nan,
        "boundary_selection_rate": boundary_selected / len(boundary) if boundary else math.nan,
        "log_mean_null_e": log_mean_null_e,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()
    if not 0.0 < args.alpha < 1.0:
        raise SystemExit("--alpha must be in (0,1)")

    grouped = defaultdict(list)
    paths = sorted(args.input.glob("*/*/rep_*/task_metrics.tsv"))
    if not paths:
        raise SystemExit("No task outputs under {}".format(args.input))
    for path in paths:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        relative = path.relative_to(args.input)
        benchmark, cell = relative.parts[0], relative.parts[1]
        for name, kappas in CALIBRATORS.items():
            metrics = task_metrics(rows, kappas, args.alpha)
            if metrics is not None:
                grouped[(benchmark, cell, name)].append(metrics)

    fields = [
        "benchmark", "cell", "calibrator", "kappas", "replicates",
        "mean_tested", "mean_selected", "mean_causal_fdr", "causal_fdr_q95",
        "mean_identifiable_fdr", "identifiable_fdr_q95",
        "mean_power", "power_q05", "mean_boundary_selection_rate",
        "max_log_mean_null_e",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for (benchmark, cell, name), records in sorted(grouped.items()):
            writer.writerow({
                "benchmark": benchmark,
                "cell": cell,
                "calibrator": name,
                "kappas": ",".join(str(value) for value in CALIBRATORS[name]),
                "replicates": len(records),
                "mean_tested": mean([record["n_tested"] for record in records]),
                "mean_selected": mean([record["n_selected"] for record in records]),
                "mean_causal_fdr": mean([record["causal_fdr"] for record in records]),
                "causal_fdr_q95": quantile(
                    [record["causal_fdr"] for record in records], 0.95
                ),
                "mean_identifiable_fdr": mean([
                    record["identifiable_fdr"] for record in records
                ]),
                "identifiable_fdr_q95": quantile([
                    record["identifiable_fdr"] for record in records
                ], 0.95),
                "mean_power": mean([record["power"] for record in records]),
                "power_q05": quantile([record["power"] for record in records], 0.05),
                "mean_boundary_selection_rate": mean([
                    record["boundary_selection_rate"] for record in records
                ]),
                "max_log_mean_null_e": max([
                    record["log_mean_null_e"] for record in records
                    if math.isfinite(record["log_mean_null_e"])
                ], default=math.nan),
            })
    print("Wrote p-to-e sensitivity for {} task families to {}".format(len(paths), args.output))


if __name__ == "__main__":
    main()
