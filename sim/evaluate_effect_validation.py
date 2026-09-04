#!/usr/bin/env python3
"""Apply the frozen pass/fail rules for post-selection effect validation."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


VALID_INDEPENDENT = (
    "strong_independent",
    "moderate_independent",
    "few_instruments_independent",
    "balanced_heterogeneity_independent",
)
VALID_BY = (
    "strong_independent",
    "balanced_heterogeneity_independent",
)
DEPENDENT = "dense_cross_protein_dependence"
DECLARED_OVERLAP = "declared_overlap"
EFFECT_LIMITS = {
    "strong_independent": 0.15,
    "moderate_independent": 0.20,
    "few_instruments_independent": 0.25,
    "balanced_heterogeneity_independent": 0.20,
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fcr-summary", type=Path, required=True)
    parser.add_argument("--effect-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-pass", action="store_true")
    args = parser.parse_args()

    fcr_rows = read_rows(args.fcr_summary)
    effect_rows = read_rows(args.effect_summary)
    fcr = {(row["cell"], row["method"]): row for row in fcr_rows}
    effects = {
        (row["cell"], row["true_scenario"], row["subset"],
         row["interval_method"], row["effect"]): row
        for row in effect_rows
    }
    decisions: list[dict[str, str]] = []

    def record(
        cell: str, method: str, criterion: str, observed: float,
        limit: float, relation: str,
    ) -> None:
        if relation == "max":
            passed = math.isfinite(observed) and observed <= limit
        elif relation == "min":
            passed = math.isfinite(observed) and observed >= limit
        elif relation == "equal":
            passed = math.isfinite(observed) and math.isclose(observed, limit)
        else:
            raise ValueError(relation)
        decisions.append({
            "cell": cell,
            "method": method,
            "criterion": criterion,
            "observed": f"{observed:.10g}",
            "relation": relation,
            "limit": f"{limit:.10g}",
            "pass": "YES" if passed else "NO",
        })

    required_fcr = {(cell, "BH") for cell in VALID_INDEPENDENT}
    required_fcr.update((cell, "BY") for cell in VALID_BY)
    required_fcr.add((DEPENDENT, "BY"))
    required_fcr.update((DECLARED_OVERLAP, method) for method in ("BH", "BY"))
    missing_fcr = sorted(required_fcr - set(fcr))
    if missing_fcr:
        raise SystemExit(f"Missing FCR summary rows: {missing_fcr}")

    for cell in VALID_INDEPENDENT:
        methods = ("BH", "BY") if cell in VALID_BY else ("BH",)
        for method in methods:
            row = fcr[(cell, method)]
            mean_fcp = float(row["mean_fcp"])
            upper = mean_fcp + 1.96 * float(row["se_mean_fcp"])
            record(cell, method, "mean_FCP", mean_fcp, 0.05, "max")
            record(cell, method, "mean_FCP_normal_upper", upper, 0.07, "max")
            record(
                cell, method, "selected_interval_complete_rate",
                float(row["interval_complete_rate"]), 1.0, "equal",
            )
            record(
                cell, method, "unbounded_effect_set_rate",
                float(row["unbounded_effect_set_rate"]), 0.01, "max",
            )
            record(
                cell, method, "families_with_calls_fraction",
                int(row["replicates_with_calls"]) / int(row["replicates"]),
                0.90, "min",
            )

    dependent = fcr[(DEPENDENT, "BY")]
    dependent_fcp = float(dependent["mean_fcp"])
    dependent_upper = dependent_fcp + 1.96 * float(dependent["se_mean_fcp"])
    record(DEPENDENT, "BY", "mean_FCP", dependent_fcp, 0.05, "max")
    record(DEPENDENT, "BY", "mean_FCP_normal_upper", dependent_upper, 0.07, "max")
    record(
        DEPENDENT, "BY", "selected_interval_complete_rate",
        float(dependent["interval_complete_rate"]), 1.0, "equal",
    )
    record(
        DEPENDENT, "BY", "unbounded_effect_set_rate",
        float(dependent["unbounded_effect_set_rate"]), 0.01, "max",
    )
    record(
        DEPENDENT, "BY", "families_with_calls_fraction",
        int(dependent["replicates_with_calls"]) / int(dependent["replicates"]),
        0.90, "min",
    )

    for method in ("BH", "BY"):
        row = fcr[(DECLARED_OVERLAP, method)]
        record(
            DECLARED_OVERLAP, method, "selected_intervals_fail_closed",
            float(row["total_finite_fcr_alpha"]), 0.0, "equal",
        )
        record(
            DECLARED_OVERLAP, method, "unresolved_overlap_status_present",
            float(row["total_unresolved_overlap"]), 1.0, "min",
        )

    required_effects = {
        (cell, "M1", subset, interval, effect)
        for cell in VALID_INDEPENDENT
        for subset, interval in (("ALL", "nominal"), ("SELECTED_BH", "nominal"))
        for effect in ("beta1", "beta2", "indirect")
    }
    missing_effects = sorted(required_effects - set(effects))
    if missing_effects:
        raise SystemExit(f"Missing effect summary rows: {missing_effects}")

    for cell in VALID_INDEPENDENT:
        limit = EFFECT_LIMITS[cell]
        for effect in ("beta1", "beta2", "indirect"):
            marginal = effects[(cell, "M1", "ALL", "nominal", effect)]
            selected_bh = effects[(cell, "M1", "SELECTED_BH", "nominal", effect)]
            record(
                cell, f"marginal_{effect}", "estimate_rate",
                float(marginal["estimate_rate"]), 0.99, "min",
            )
            record(
                cell, f"marginal_{effect}", "nominal_coverage",
                float(marginal["coverage"]), 0.925, "min",
            )
            record(
                cell, f"marginal_{effect}", "absolute_standardized_bias",
                float(marginal["standardized_bias"]), limit, "max",
            )
            record(
                cell, f"selected_BH_{effect}", "absolute_standardized_bias",
                float(selected_bh["standardized_bias"]), 0.25, "max",
            )

    overall = all(row["pass"] == "YES" for row in decisions)
    decisions.append({
        "cell": "ALL",
        "method": "ALL",
        "criterion": "overall_frozen_decision",
        "observed": "1" if overall else "0",
        "relation": "equal",
        "limit": "1",
        "pass": "YES" if overall else "NO",
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(decisions[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(decisions)
    print(f"Wrote {args.output}: overall={'PASS' if overall else 'FAIL'}")
    if args.require_pass and not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
