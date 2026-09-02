from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


Record = Tuple[float, float, int, bool, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate analytical conjunction-null mediation calibration.")
    parser.add_argument("--outdir", type=Path, required=True, help="Simulation output directory.")
    parser.add_argument("--output-name", default="analytic_conjunction_test", help="Subdirectory under summary/.")
    parser.add_argument("--thresholds", default="0.01,0.05,0.10", help="Comma-separated q-value thresholds.")
    return parser.parse_args()


def to_float(value, default: float = 1.0) -> float:
    try:
        x = float(value)
        if math.isfinite(x):
            return x
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


def bh_qvalues(pvals: Sequence[float]) -> List[float]:
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    qvals = [1.0] * n
    previous = 1.0
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        value = min(previous, pvals[idx] * n / rank)
        qvals[idx] = min(1.0, value)
        previous = value
    return qvals


def load_records(outdir: Path) -> Dict[Tuple[str, str, str], List[Record]]:
    datasets: Dict[Tuple[str, str, str], List[Record]] = defaultdict(list)
    for benchmark, path in iter_metric_paths(outdir):
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                key = (row.get("benchmark", benchmark) or benchmark, row.get("cell", ""), row.get("replicate", ""))
                p_alpha = to_float(row.get("ivw_rf_to_pp_p"), 1.0)
                p_beta = to_float(row.get("ivw_pp_to_outcome_p"), 1.0)
                is_true = to_int(row.get("is_true_m1", 0), 0)
                direction_consistent = str(row.get("direction_consistent", "")).upper() == "YES"
                p_m1 = to_float(row.get("P_M1", row.get("p_m1", 0.0)), 0.0)
                datasets[key].append((p_alpha, p_beta, is_true, direction_consistent, p_m1))
    return datasets


def include_variant(name: str, rec: Record) -> bool:
    _, _, _, direction_consistent, p_m1 = rec
    if name == "conjunction":
        return True
    if name == "conjunction_dir_consistent":
        return direction_consistent
    if name == "conjunction_dir_pm1_0.5":
        return direction_consistent and p_m1 >= 0.5
    if name == "conjunction_dir_pm1_0.9":
        return direction_consistent and p_m1 >= 0.9
    raise ValueError(f"Unknown variant: {name}")


def main() -> None:
    args = parse_args()
    thresholds = tuple(float(x) for x in args.thresholds.split(",") if x)
    variants = (
        "conjunction",
        "conjunction_dir_consistent",
        "conjunction_dir_pm1_0.5",
        "conjunction_dir_pm1_0.9",
    )
    datasets = load_records(args.outdir)
    aggregate = defaultdict(lambda: {"selected": 0, "tp": 0, "fp": 0, "power_sum": 0.0, "power_n": 0, "reps": 0, "reps_with_calls": 0})

    for (benchmark, cell, _replicate), records in datasets.items():
        p_med = [max(rec[0], rec[1]) for rec in records]
        q_med = bh_qvalues(p_med)
        total_true = sum(1 for rec in records if rec[2] == 1)
        for variant in variants:
            keep = [include_variant(variant, rec) for rec in records]
            for threshold in thresholds:
                selected = [i for i, q in enumerate(q_med) if q <= threshold and keep[i]]
                tp = sum(1 for i in selected if records[i][2] == 1)
                fp = len(selected) - tp
                key = (benchmark, cell, variant, threshold)
                aggregate[key]["selected"] += len(selected)
                aggregate[key]["tp"] += tp
                aggregate[key]["fp"] += fp
                aggregate[key]["power_sum"] += (tp / total_true) if total_true else 0.0
                aggregate[key]["power_n"] += 1 if total_true else 0
                aggregate[key]["reps"] += 1
                aggregate[key]["reps_with_calls"] += 1 if selected else 0

    outdir = args.outdir / "summary" / args.output_name
    outdir.mkdir(parents=True, exist_ok=True)
    fields = [
        "benchmark",
        "cell",
        "method",
        "threshold",
        "replicates",
        "replicates_with_calls",
        "n_selected",
        "true_positives",
        "false_positives",
        "empirical_fdr",
        "mean_replicate_power",
    ]
    summary_path = outdir / "analytic_conjunction_fdr_power.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(aggregate):
            vals = aggregate[key]
            writer.writerow(
                {
                    "benchmark": key[0],
                    "cell": key[1],
                    "method": key[2],
                    "threshold": key[3],
                    "replicates": vals["reps"],
                    "replicates_with_calls": vals["reps_with_calls"],
                    "n_selected": vals["selected"],
                    "true_positives": vals["tp"],
                    "false_positives": vals["fp"],
                    "empirical_fdr": (vals["fp"] / vals["selected"]) if vals["selected"] else "",
                    "mean_replicate_power": (vals["power_sum"] / vals["power_n"]) if vals["power_n"] else "",
                }
            )
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
