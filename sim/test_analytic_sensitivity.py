from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


Record = Tuple[float, float, float, int, bool, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate sensitivity-bounded analytical mediation p-values."
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--output-name", default="analytic_sensitivity_test", help="Subdirectory under summary/."
    )
    parser.add_argument(
        "--bias-bounds",
        default="0,0.05,0.10,0.15,0.20,0.25",
        help="Comma-separated absolute protein-to-outcome bias bounds (Delta).",
    )
    parser.add_argument("--thresholds", default="0.01,0.05,0.10")
    return parser.parse_args()


def to_float(value, default: float = float("nan")) -> float:
    try:
        result = float(value)
        if math.isfinite(result):
            return result
    except Exception:
        pass
    return default


def to_int(value, default: int = 0) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def iter_metric_paths(outdir: Path) -> Iterable[Tuple[str, Path]]:
    for benchmark in ("classification", "calibration"):
        path = outdir / "summary" / benchmark / "protein_level_metrics.tsv"
        if path.exists():
            yield benchmark, path


def adjust_pvalues(pvals: Sequence[float], method: str) -> List[float]:
    n = len(pvals)
    if n == 0:
        return []
    harmonic = sum(1.0 / i for i in range(1, n + 1)) if method == "BY" else 1.0
    order = sorted(range(n), key=lambda i: pvals[i])
    adjusted = [1.0] * n
    previous = 1.0
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        value = min(previous, pvals[idx] * n * harmonic / rank)
        adjusted[idx] = min(1.0, value)
        previous = value
    return adjusted


def normal_two_sided_p(z: float) -> float:
    if not math.isfinite(z):
        return 1.0
    return min(1.0, math.erfc(abs(z) / math.sqrt(2.0)))


def interval_null_p(beta: float, se: float, bias_bound: float) -> float:
    """Conservative p-value for H0: |beta| <= bias_bound.

    The doubled least-favourable-tail p-value is conservative at either boundary
    of the composite null and avoids the doubled type-I error at Delta=0.
    """
    if not math.isfinite(beta) or not math.isfinite(se) or se <= 0.0:
        return 1.0
    z_beyond_bound = (abs(beta) - bias_bound) / se
    if z_beyond_bound <= 0.0:
        return 1.0
    return normal_two_sided_p(z_beyond_bound)


def load_records(outdir: Path) -> Dict[Tuple[str, str, str], List[Record]]:
    datasets: Dict[Tuple[str, str, str], List[Record]] = defaultdict(list)
    for benchmark, path in iter_metric_paths(outdir):
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                key = (
                    row.get("benchmark", benchmark) or benchmark,
                    row.get("cell", ""),
                    row.get("replicate", ""),
                )
                p_alpha = to_float(row.get("ivw_rf_to_pp_p"), 1.0)
                beta = to_float(row.get("ivw_pp_to_outcome_beta"))
                se = to_float(row.get("ivw_pp_to_outcome_se"))
                is_true = to_int(row.get("is_true_m1", 0), 0)
                direction_consistent = str(row.get("direction_consistent", "")).upper() == "YES"
                p_m1 = to_float(row.get("P_M1", row.get("p_m1", 0.0)), 0.0)
                datasets[key].append((p_alpha, beta, se, is_true, direction_consistent, p_m1))
    return datasets


def passes_gate(gate: str, record: Record) -> bool:
    direction_consistent = record[4]
    p_m1 = record[5]
    if gate == "standalone":
        return True
    if gate == "direction":
        return direction_consistent
    if gate == "direction_pm1_0.5":
        return direction_consistent and p_m1 >= 0.5
    if gate == "direction_pm1_0.9":
        return direction_consistent and p_m1 >= 0.9
    raise ValueError(f"Unknown gate: {gate}")


def main() -> None:
    args = parse_args()
    bounds = tuple(float(value) for value in args.bias_bounds.split(",") if value)
    thresholds = tuple(float(value) for value in args.thresholds.split(",") if value)
    gates = ("standalone", "direction", "direction_pm1_0.5", "direction_pm1_0.9")
    adjustments = ("BH", "BY")
    datasets = load_records(args.outdir)
    aggregate = defaultdict(
        lambda: {
            "selected": 0,
            "tp": 0,
            "fp": 0,
            "power_sum": 0.0,
            "power_n": 0,
            "reps": 0,
            "reps_with_calls": 0,
        }
    )

    for (benchmark, cell, _replicate), records in datasets.items():
        total_true = sum(record[3] for record in records)
        for bias_bound in bounds:
            # Intersection-union p-value for H0: alpha=0 OR |beta|<=Delta.
            p_med = [
                max(record[0], interval_null_p(record[1], record[2], bias_bound))
                for record in records
            ]
            for adjustment in adjustments:
                q_med = adjust_pvalues(p_med, adjustment)
                for gate in gates:
                    keep = [passes_gate(gate, record) for record in records]
                    for threshold in thresholds:
                        selected = [
                            idx
                            for idx, q_value in enumerate(q_med)
                            if q_value <= threshold and keep[idx]
                        ]
                        tp = sum(records[idx][3] for idx in selected)
                        fp = len(selected) - tp
                        key = (benchmark, cell, bias_bound, adjustment, gate, threshold)
                        values = aggregate[key]
                        values["selected"] += len(selected)
                        values["tp"] += tp
                        values["fp"] += fp
                        values["power_sum"] += (tp / total_true) if total_true else 0.0
                        values["power_n"] += 1 if total_true else 0
                        values["reps"] += 1
                        values["reps_with_calls"] += 1 if selected else 0

    destination = args.outdir / "summary" / args.output_name
    destination.mkdir(parents=True, exist_ok=True)
    summary_path = destination / "analytic_sensitivity_fdr_power.tsv"
    fields = [
        "benchmark",
        "cell",
        "bias_bound",
        "adjustment",
        "gate",
        "threshold",
        "replicates",
        "replicates_with_calls",
        "n_selected",
        "true_positives",
        "false_positives",
        "empirical_fdr",
        "mean_replicate_power",
    ]
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(aggregate):
            values = aggregate[key]
            writer.writerow(
                {
                    "benchmark": key[0],
                    "cell": key[1],
                    "bias_bound": key[2],
                    "adjustment": key[3],
                    "gate": key[4],
                    "threshold": key[5],
                    "replicates": values["reps"],
                    "replicates_with_calls": values["reps_with_calls"],
                    "n_selected": values["selected"],
                    "true_positives": values["tp"],
                    "false_positives": values["fp"],
                    "empirical_fdr": (
                        values["fp"] / values["selected"] if values["selected"] else ""
                    ),
                    "mean_replicate_power": (
                        values["power_sum"] / values["power_n"] if values["power_n"] else ""
                    ),
                }
            )
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
