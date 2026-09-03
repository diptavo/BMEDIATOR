#!/usr/bin/env python3
"""Summarize factorized two-stage evidence from independent-protein tasks."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    middle = len(finite) // 2
    return finite[middle] if len(finite) % 2 else 0.5 * (finite[middle - 1] + finite[middle])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    paths = sorted(args.input.glob("classification/*/rep_*/task_metrics.tsv"))
    if not paths:
        raise SystemExit(f"No task_metrics.tsv files under {args.input}")
    for path in paths:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                grouped[(row["cell"], row["true_scenario"])].append(row)

    fields = [
        "cell", "true_scenario", "n", "mean_beta1", "mean_beta2",
        "median_p_XM", "median_p_MY", "median_conjunction_p",
        "reject_conjunction_0.01", "reject_conjunction_0.05",
        "reject_conjunction_0.10", "BF_XM_gt_10", "BF_MY_gt_10",
        "both_BF_gt_10",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (cell, scenario), rows in sorted(grouped.items()):
            p_xm = [as_float(row.get("factor_p_XM", "")) for row in rows]
            p_my = [as_float(row.get("factor_p_MY", "")) for row in rows]
            p_conj = [as_float(row.get("factor_conjunction_p", "")) for row in rows]
            log_bf_xm = [as_float(row.get("factor_log_BF_XM", "")) for row in rows]
            log_bf_my = [as_float(row.get("factor_log_BF_MY", "")) for row in rows]
            paired_bf = [
                (as_float(row.get("factor_log_BF_XM", "")),
                 as_float(row.get("factor_log_BF_MY", "")))
                for row in rows
            ]
            valid = [value for value in p_conj if math.isfinite(value)]
            writer.writerow({
                "cell": cell,
                "true_scenario": scenario,
                "n": len(rows),
                "mean_beta1": f"{mean([as_float(r.get('factor_beta1', '')) for r in rows]):.6g}",
                "mean_beta2": f"{mean([as_float(r.get('factor_beta2', '')) for r in rows]):.6g}",
                "median_p_XM": f"{median(p_xm):.6g}",
                "median_p_MY": f"{median(p_my):.6g}",
                "median_conjunction_p": f"{median(p_conj):.6g}",
                "reject_conjunction_0.01": f"{mean([float(p <= 0.01) for p in valid]):.6g}",
                "reject_conjunction_0.05": f"{mean([float(p <= 0.05) for p in valid]):.6g}",
                "reject_conjunction_0.10": f"{mean([float(p <= 0.10) for p in valid]):.6g}",
                "BF_XM_gt_10": f"{mean([float(v > math.log(10.0)) for v in log_bf_xm if math.isfinite(v)]):.6g}",
                "BF_MY_gt_10": f"{mean([float(v > math.log(10.0)) for v in log_bf_my if math.isfinite(v)]):.6g}",
                "both_BF_gt_10": f"{mean([float(a > math.log(10.0) and b > math.log(10.0)) for a, b in paired_bf if math.isfinite(a) and math.isfinite(b)]):.6g}",
            })
    print(f"Wrote {args.output} from {len(paths)} tasks")


if __name__ == "__main__":
    main()
