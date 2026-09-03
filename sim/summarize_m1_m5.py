#!/usr/bin/env python3
"""Summarize the dedicated M1-versus-M5 identification study."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def mean(values: list[float]) -> float | None:
    finite = [value for value in values if math.isfinite(value)]
    return statistics.fmean(finite) if finite else None


def rate(flags: list[bool]) -> float | None:
    return statistics.fmean(int(flag) for flag in flags) if flags else None


def wilson(flags: list[bool]) -> tuple[float | None, float | None]:
    if not flags:
        return None, None
    n = len(flags)
    p = sum(int(flag) for flag in flags) / n
    z = 1.959963984540054
    denominator = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denominator
    radius = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denominator
    return center - radius, center + radius


def expected_calibration_error(labels: list[int], scores: list[float], bins: int = 10) -> float | None:
    if not scores:
        return None
    total = len(scores)
    result = 0.0
    for index in range(bins):
        lower = index / bins
        upper = (index + 1) / bins
        members = [
            (label, score) for label, score in zip(labels, scores)
            if lower <= score <= upper and (index == bins - 1 or score < upper)
        ]
        if not members:
            continue
        observed = statistics.fmean(label for label, _ in members)
        predicted = statistics.fmean(score for _, score in members)
        result += len(members) / total * abs(observed - predicted)
    return result


def auc(labels: list[int], scores: list[float]) -> float | None:
    positives = [score for label, score in zip(labels, scores) if label == 1]
    negatives = [score for label, score in zip(labels, scores) if label == 0]
    if not positives or not negatives:
        return None
    wins = 0.0
    for positive in positives:
        for negative in negatives:
            wins += float(positive > negative) + 0.5 * float(positive == negative)
    return wins / (len(positives) * len(negatives))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_rows(outdir: Path) -> list[dict]:
    rows = []
    task_paths = sorted((outdir / "classification").glob("*/rep_*/task_metrics.tsv"))
    aggregate = outdir / "summary" / "classification" / "protein_level_metrics.tsv"
    paths = task_paths if task_paths else ([aggregate] if aggregate.exists() else [])
    for path in paths:
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("true_scenario") not in {"M1", "M5"}:
                    continue
                p1 = as_float(row, "P_M1")
                p5 = as_float(row, "P_M5")
                denominator = p1 + p5
                pair_score = p1 / denominator if denominator > 0.0 else float("nan")
                row["p1"] = p1
                row["p5"] = p5
                row["pair_score"] = pair_score
                row["cell"] = row.get("cell") or path.parents[1].name
                row["identification_class"] = row.get("identification_class") or "identified"
                rows.append(row)
    return rows


def summarize_group(cell: str, rows: list[dict]) -> dict:
    identifiable = all(row["identification_class"] != "nonidentifiable" for row in rows)
    scored = [row for row in rows if math.isfinite(row["pair_score"])]
    m1 = [row for row in scored if row["true_scenario"] == "M1"]
    m5 = [row for row in scored if row["true_scenario"] == "M5"]
    all_m1 = [row for row in rows if row["true_scenario"] == "M1"]
    all_m5 = [row for row in rows if row["true_scenario"] == "M5"]
    labels = [int(row["true_scenario"] == "M1") for row in scored]
    scores = [row["pair_score"] for row in scored]
    pair_correct = [
        (row["pair_score"] >= 0.5) == (row["true_scenario"] == "M1")
        for row in scored
    ]
    brier = mean([(score - label) ** 2 for score, label in zip(scores, labels)])
    log_loss = mean([
        -(label * math.log(max(score, 1e-12))
          + (1 - label) * math.log(max(1.0 - score, 1e-12)))
        for score, label in zip(scores, labels)
    ])
    classes = sorted(set(row["identification_class"] for row in rows))
    m1_calls = [row["pair_score"] >= 0.5 for row in m1]
    m5_calls = [row["pair_score"] < 0.5 for row in m5]
    all_correct = [row.get("pred_scenario") == row["true_scenario"] for row in rows]
    all_m1_correct = [row.get("pred_scenario") == "M1" for row in all_m1]
    all_m5_correct = [row.get("pred_scenario") == "M5" for row in all_m5]
    pair_ci = wilson(pair_correct)
    m1_ci = wilson(m1_calls)
    m5_ci = wilson(m5_calls)
    all_ci = wilson(all_correct)
    all_m1_ci = wilson(all_m1_correct)
    all_m5_ci = wilson(all_m5_correct)
    boundary_calls = [row["pair_score"] >= 0.5 for row in m5]
    boundary_ci = wilson(boundary_calls)
    return {
        "cell": cell,
        "identification_class": ",".join(classes),
        "n": len(rows),
        "n_scored": len(scored),
        "pair_scored_fraction": len(scored) / len(rows) if rows else None,
        "n_m1": len(m1),
        "n_m5": len(m5),
        "pair_auc": auc(labels, scores) if identifiable else None,
        "pair_accuracy": rate(pair_correct) if identifiable else None,
        "pair_accuracy_ci_low": pair_ci[0] if identifiable else None,
        "pair_accuracy_ci_high": pair_ci[1] if identifiable else None,
        "m1_sensitivity": rate(m1_calls) if identifiable else None,
        "m1_sensitivity_ci_low": m1_ci[0] if identifiable else None,
        "m1_sensitivity_ci_high": m1_ci[1] if identifiable else None,
        "m5_specificity": rate(m5_calls) if identifiable else None,
        "m5_specificity_ci_low": m5_ci[0] if identifiable else None,
        "m5_specificity_ci_high": m5_ci[1] if identifiable else None,
        "all_state_accuracy": rate(all_correct) if identifiable else None,
        "all_state_accuracy_ci_low": all_ci[0] if identifiable else None,
        "all_state_accuracy_ci_high": all_ci[1] if identifiable else None,
        "all_state_m1_recall": rate(all_m1_correct) if identifiable else None,
        "all_state_m1_recall_ci_low": all_m1_ci[0] if identifiable else None,
        "all_state_m1_recall_ci_high": all_m1_ci[1] if identifiable else None,
        "all_state_m5_recall": rate(all_m5_correct) if identifiable else None,
        "all_state_m5_recall_ci_low": all_m5_ci[0] if identifiable else None,
        "all_state_m5_recall_ci_high": all_m5_ci[1] if identifiable else None,
        "pair_brier": brier if identifiable else None,
        "pair_log_loss": log_loss if identifiable else None,
        "pair_ece_10": expected_calibration_error(labels, scores) if identifiable else None,
        "mean_pair_score_m1": mean([row["pair_score"] for row in m1]),
        "mean_pair_score_m5": mean([row["pair_score"] for row in m5]),
        "m1_p1_ge_0_8": rate([row["p1"] >= 0.8 for row in m1]),
        "m5_false_p1_ge_0_8": rate([row["p1"] >= 0.8 for row in m5]),
        "m5_p5_ge_0_8": rate([row["p5"] >= 0.8 for row in m5]),
        "boundary_m1_pair_call_rate": rate(boundary_calls) if not identifiable else None,
        "boundary_m1_pair_call_ci_low": boundary_ci[0] if not identifiable else None,
        "boundary_m1_pair_call_ci_high": boundary_ci[1] if not identifiable else None,
    }


def calibration_rows(rows: list[dict]) -> list[dict]:
    usable = [
        row for row in rows
        if row["identification_class"] != "nonidentifiable" and math.isfinite(row["pair_score"])
    ]
    output = []
    for index in range(10):
        lower = index / 10.0
        upper = (index + 1) / 10.0
        subset = [
            row for row in usable
            if lower <= row["pair_score"] <= upper
            if index == 9 or row["pair_score"] < upper
        ]
        if not subset:
            continue
        output.append({
            "bin": index,
            "lower": lower,
            "upper": upper,
            "n": len(subset),
            "mean_pair_score": mean([row["pair_score"] for row in subset]),
            "observed_m1_fraction": rate([row["true_scenario"] == "M1" for row in subset]),
        })
    return output


def main() -> None:
    args = parse_args()
    rows = load_rows(args.outdir)
    if not rows:
        raise SystemExit(f"No M1/M5 task metrics found under {args.outdir}")

    by_cell: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_cell[row["cell"]].append(row)

    summaries = [summarize_group(cell, by_cell[cell]) for cell in sorted(by_cell)]
    identified_rows = [row for row in rows if row["identification_class"] != "nonidentifiable"]
    if identified_rows:
        summaries.append(summarize_group("ALL_IDENTIFIABLE", identified_rows))

    output_dir = args.outdir / "summary" / "m1_m5"
    write_tsv(output_dir / "discrimination.tsv", summaries)
    write_tsv(output_dir / "pair_calibration.tsv", calibration_rows(rows))
    print(f"Wrote {output_dir / 'discrimination.tsv'}")
    print(f"Wrote {output_dir / 'pair_calibration.tsv'}")


if __name__ == "__main__":
    main()
