#!/usr/bin/env python3
"""Evaluate the frozen balanced partial-conjunction BH validation criteria."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


ALPHA = 0.05
MIXED_VALID = ("broad_balanced", "narrow_balanced", "rare_mediation")
PURE_NULL = ("least_favorable_null", "global_null")
STRESS = (
    "balanced_heterogeneity",
    "sparse_outliers",
    "dense_cross_protein_dependence",
)


def wilson_upper(successes: int, total: int, z: float = 1.6448536269514722) -> float:
    if total <= 0:
        return math.nan
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = proportion + z * z / (2.0 * total)
    radius = z * math.sqrt(
        proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)
    )
    return (center + radius) / denominator


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--method", default="balanced_BH")
    parser.add_argument("--threshold", type=float, default=ALPHA)
    parser.add_argument("--require-pass", action="store_true")
    args = parser.parse_args()

    with args.summary.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = {
        row["cell"]: row
        for row in rows
        if row["method"] == args.method
        and math.isclose(float(row["threshold"]), args.threshold)
    }
    required = set(MIXED_VALID + PURE_NULL + STRESS)
    missing = sorted(required - set(selected))
    if missing:
        raise SystemExit(f"Missing {args.method} rows for: {', '.join(missing)}")

    decisions: list[dict[str, object]] = []

    def record(cell: str, criterion: str, observed: float, limit: float, passed: bool) -> None:
        decisions.append({
            "cell": cell,
            "criterion": criterion,
            "observed": f"{observed:.8g}",
            "limit": f"{limit:.8g}",
            "pass": "YES" if passed else "NO",
        })

    for cell in MIXED_VALID:
        row = selected[cell]
        mean_fdp = float(row["mean_fdp"])
        upper = mean_fdp + 1.96 * float(row["se_mean_fdp"])
        record(cell, "mean_FDP_max", mean_fdp, ALPHA, mean_fdp <= ALPHA)
        record(cell, "mean_FDP_normal_upper_max", upper, 0.06, upper <= 0.06)

    for cell in PURE_NULL:
        row = selected[cell]
        replicates = int(row["replicates"])
        families_with_calls = int(row["replicates_with_calls"])
        estimate = families_with_calls / replicates
        upper = wilson_upper(families_with_calls, replicates)
        record(cell, "pure_null_FDR_max", estimate, ALPHA, estimate <= ALPHA)
        record(cell, "pure_null_Wilson_upper_max", upper, 0.065, upper <= 0.065)

    for cell in STRESS:
        mean_fdp = float(selected[cell]["mean_fdp"])
        record(cell, "stress_mean_FDP_max", mean_fdp, ALPHA, mean_fdp <= ALPHA)

    for cell, minimum in (("broad_balanced", 0.80), ("narrow_balanced", 0.50)):
        power = float(selected[cell]["mean_power"])
        record(cell, "mean_power_min", power, minimum, power >= minimum)

    overall = all(row["pass"] == "YES" for row in decisions)
    decisions.append({
        "cell": "ALL",
        "criterion": "overall_frozen_decision",
        "observed": "1" if overall else "0",
        "limit": "1",
        "pass": "YES" if overall else "NO",
    })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("cell", "criterion", "observed", "limit", "pass"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(decisions)
    print(f"Wrote {args.output}: overall={'PASS' if overall else 'FAIL'}")
    if args.require_pass and not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
