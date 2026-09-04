#!/usr/bin/env python3
"""Summarize marginal and post-selection effect-estimation validation."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from summarize_factorized_calibration import validate_complete_run


ALPHA = 0.05
METHODS = {
    "BH": {
        "q": "factor_balanced_conjunction_q_BH",
        "alpha": "factor_fcr_alpha_BH",
        "status": "factor_fcr_bh_status",
        "expected_status": "PROFILE_FCR_BH_INDEPENDENCE_CONDITIONAL",
        "intervals": {
            "beta1": ("factor_beta1_fcr_ci_lower_BH", "factor_beta1_fcr_ci_upper_BH"),
            "beta2": ("factor_beta2_fcr_ci_lower_BH", "factor_beta2_fcr_ci_upper_BH"),
            "indirect": ("factor_indirect_fcr_ci_lower_BH", "factor_indirect_fcr_ci_upper_BH"),
        },
    },
    "BY": {
        "q": "factor_balanced_conjunction_q_BY",
        "alpha": "factor_fcr_alpha_BY",
        "status": "factor_fcr_by_status",
        "expected_status": "PROFILE_FCR_BY_DEPENDENCE_CONDITIONAL",
        "intervals": {
            "beta1": ("factor_beta1_fcr_ci_lower_BY", "factor_beta1_fcr_ci_upper_BY"),
            "beta2": ("factor_beta2_fcr_ci_lower_BY", "factor_beta2_fcr_ci_upper_BY"),
            "indirect": ("factor_indirect_fcr_ci_lower_BY", "factor_indirect_fcr_ci_upper_BY"),
        },
    },
}
EFFECTS = {
    "beta1": {
        "estimate": "factor_beta1",
        "se": "factor_beta1_se",
        "truth": "true_beta1",
        "nominal": ("factor_beta1_ci_lower", "factor_beta1_ci_upper"),
    },
    "beta2": {
        "estimate": "factor_beta2",
        "se": "factor_beta2_se",
        "truth": "true_beta2",
        "nominal": ("factor_beta2_ci_lower", "factor_beta2_ci_upper"),
    },
    "indirect": {
        "estimate": "factor_indirect",
        "se": "factor_indirect_se",
        "truth": "true_mediated_effect",
        "nominal": ("factor_indirect_ci_lower", "factor_indirect_ci_upper"),
    },
}


def as_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def sample_sd(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if len(finite) < 2:
        return math.nan
    center = mean(finite)
    return math.sqrt(sum((value - center) ** 2 for value in finite) / (len(finite) - 1))


def quantile(values: list[float], probability: float) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    position = probability * (len(finite) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    weight = position - lower
    return finite[lower] * (1.0 - weight) + finite[upper] * weight


def selected(row: dict[str, str], q_column: str) -> bool:
    value = as_float(row.get(q_column))
    return math.isfinite(value) and value <= ALPHA


def interval_covers(
    row: dict[str, str], lower: str, upper: str, truth: str
) -> bool | None:
    values = (as_float(row.get(lower)), as_float(row.get(upper)), as_float(row.get(truth)))
    if math.isnan(values[0]) or math.isnan(values[1]) or not math.isfinite(values[2]):
        return None
    return values[0] <= values[2] <= values[1]


def joint_coverage(
    row: dict[str, str], intervals: dict[str, tuple[str, str]]
) -> bool | None:
    coverage = [
        interval_covers(row, *intervals[name], EFFECTS[name]["truth"])
        for name in EFFECTS
    ]
    if any(value is None for value in coverage):
        return None
    return all(bool(value) for value in coverage)


NOMINAL_INTERVALS = {
    name: value["nominal"] for name, value in EFFECTS.items()
}


def replicate_metrics(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    benchmark = rows[0]["benchmark"]
    cell = rows[0]["cell"]
    replicate = rows[0]["replicate"]
    output: list[dict[str, object]] = []
    tested = [
        row for row in rows
        if math.isfinite(as_float(row.get("factor_balanced_conjunction_p")))
    ]
    for method, spec in METHODS.items():
        selected_rows = [row for row in tested if selected(row, spec["q"])]
        fcr_joint = [joint_coverage(row, spec["intervals"]) for row in selected_rows]
        nominal_joint = [joint_coverage(row, NOMINAL_INTERVALS) for row in selected_rows]
        complete = sum(value is not None for value in fcr_joint)
        false_cover = sum(value is not True for value in fcr_joint)
        nominal_false = sum(value is not True for value in nominal_joint)
        status_counts = Counter(str(row.get(spec["status"], "")) for row in rows)
        n_selected = len(selected_rows)
        output.append({
            "benchmark": benchmark,
            "cell": cell,
            "replicate": replicate,
            "method": method,
            "n_proteins": len(rows),
            "n_tested": len(tested),
            "n_selected": n_selected,
            "n_selected_M1": sum(row.get("true_scenario") == "M1" for row in selected_rows),
            "n_interval_complete": complete,
            "interval_complete_rate": complete / n_selected if n_selected else 1.0,
            "n_false_cover": false_cover,
            "false_coverage_proportion": false_cover / n_selected if n_selected else 0.0,
            "n_nominal_false_cover": nominal_false,
            "nominal_false_coverage_proportion": nominal_false / n_selected if n_selected else 0.0,
            "mean_fcr_alpha": mean([as_float(row.get(spec["alpha"])) for row in selected_rows]),
            "mean_fcr_interval_width_beta1": mean([
                as_float(row.get(spec["intervals"]["beta1"][1])) -
                as_float(row.get(spec["intervals"]["beta1"][0]))
                for row in selected_rows
            ]),
            "mean_fcr_interval_width_beta2": mean([
                as_float(row.get(spec["intervals"]["beta2"][1])) -
                as_float(row.get(spec["intervals"]["beta2"][0]))
                for row in selected_rows
            ]),
            "mean_fcr_interval_width_indirect": mean([
                as_float(row.get(spec["intervals"]["indirect"][1])) -
                as_float(row.get(spec["intervals"]["indirect"][0]))
                for row in selected_rows
            ]),
            "n_expected_status": status_counts[spec["expected_status"]],
            "n_unresolved_overlap": status_counts["UNRESOLVED_SAMPLE_OVERLAP"],
            "n_unresolved_effect": status_counts["UNRESOLVED_FCR_EFFECT_ESTIMATION"],
            "n_unbounded_effect_set": status_counts["UNBOUNDED_FCR_EFFECT_SET"],
            "n_finite_fcr_alpha": sum(
                math.isfinite(as_float(row.get(spec["alpha"]))) for row in rows
            ),
        })
    return output


def summarize_replicates(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["benchmark"]), str(row["cell"]), str(row["method"]))].append(row)
    output: list[dict[str, object]] = []
    for (benchmark, cell, method), group in sorted(grouped.items()):
        fcp = [as_float(row["false_coverage_proportion"]) for row in group]
        nominal_fcp = [as_float(row["nominal_false_coverage_proportion"]) for row in group]
        fcp_sd = sample_sd(fcp)
        total_selected = sum(int(row["n_selected"]) for row in group)
        total_false = sum(int(row["n_false_cover"]) for row in group)
        total_complete = sum(int(row["n_interval_complete"]) for row in group)
        output.append({
            "benchmark": benchmark,
            "cell": cell,
            "method": method,
            "replicates": len(group),
            "replicates_with_calls": sum(int(row["n_selected"]) > 0 for row in group),
            "mean_tested": mean([as_float(row["n_tested"]) for row in group]),
            "mean_selected": mean([as_float(row["n_selected"]) for row in group]),
            "mean_selected_M1": mean([as_float(row["n_selected_M1"]) for row in group]),
            "total_selected": total_selected,
            "interval_complete_rate": total_complete / total_selected if total_selected else 1.0,
            "mean_fcp": mean(fcp),
            "se_mean_fcp": fcp_sd / math.sqrt(len(group)) if math.isfinite(fcp_sd) else math.nan,
            "fcp_q95": quantile(fcp, 0.95),
            "pooled_noncoverage": total_false / total_selected if total_selected else 0.0,
            "mean_nominal_fcp": mean(nominal_fcp),
            "mean_fcr_alpha": mean([as_float(row["mean_fcr_alpha"]) for row in group]),
            "mean_fcr_interval_width_beta1": mean([
                as_float(row["mean_fcr_interval_width_beta1"]) for row in group
            ]),
            "mean_fcr_interval_width_beta2": mean([
                as_float(row["mean_fcr_interval_width_beta2"]) for row in group
            ]),
            "mean_fcr_interval_width_indirect": mean([
                as_float(row["mean_fcr_interval_width_indirect"]) for row in group
            ]),
            "total_expected_status": sum(int(row["n_expected_status"]) for row in group),
            "total_unresolved_overlap": sum(int(row["n_unresolved_overlap"]) for row in group),
            "total_unresolved_effect": sum(int(row["n_unresolved_effect"]) for row in group),
            "total_unbounded_effect_set": sum(
                int(row["n_unbounded_effect_set"]) for row in group
            ),
            "unbounded_effect_set_rate": (
                sum(int(row["n_unbounded_effect_set"]) for row in group) /
                total_selected if total_selected else 0.0
            ),
            "total_finite_fcr_alpha": sum(int(row["n_finite_fcr_alpha"]) for row in group),
        })
    return output


def effect_key_rows(rows: list[dict[str, str]]):
    for row in rows:
        yield "ALL", "nominal", NOMINAL_INTERVALS, row
        for method, spec in METHODS.items():
            if selected(row, spec["q"]):
                yield f"SELECTED_{method}", "nominal", NOMINAL_INTERVALS, row
                yield f"SELECTED_{method}", f"FCR_{method}", spec["intervals"], row


def summarize_effects(paths: list[Path]) -> list[dict[str, object]]:
    accumulators: dict[tuple[str, str, str, str, str, str], dict[str, object]] = {}
    for path in paths:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for subset, interval_method, intervals, row in effect_key_rows(rows):
            for effect, spec in EFFECTS.items():
                key = (
                    row["benchmark"], row["cell"], row["true_scenario"],
                    subset, interval_method, effect,
                )
                acc = accumulators.setdefault(key, {
                    "identification_class": row.get("identification_class", ""),
                    "n": 0,
                    "errors": [],
                    "ses": [],
                    "cover": 0,
                    "interval_n": 0,
                    "widths": [],
                })
                acc["n"] = int(acc["n"]) + 1
                estimate = as_float(row.get(spec["estimate"]))
                truth = as_float(row.get(spec["truth"]))
                se = as_float(row.get(spec["se"]))
                if math.isfinite(estimate) and math.isfinite(truth):
                    acc["errors"].append(estimate - truth)
                if math.isfinite(se) and se > 0.0:
                    acc["ses"].append(se)
                lower, upper = intervals[effect]
                coverage = interval_covers(row, lower, upper, spec["truth"])
                if coverage is not None:
                    acc["interval_n"] = int(acc["interval_n"]) + 1
                    acc["cover"] = int(acc["cover"]) + int(coverage)
                    acc["widths"].append(as_float(row.get(upper)) - as_float(row.get(lower)))

    output: list[dict[str, object]] = []
    for key, acc in sorted(accumulators.items()):
        benchmark, cell, scenario, subset, interval_method, effect = key
        errors = list(acc["errors"])
        bias = mean(errors)
        empirical_sd = sample_sd(errors)
        interval_n = int(acc["interval_n"])
        output.append({
            "benchmark": benchmark,
            "cell": cell,
            "true_scenario": scenario,
            "identification_class": acc["identification_class"],
            "subset": subset,
            "interval_method": interval_method,
            "effect": effect,
            "n": acc["n"],
            "n_estimated": len(errors),
            "estimate_rate": len(errors) / int(acc["n"]) if int(acc["n"]) else math.nan,
            "bias": bias,
            "empirical_sd": empirical_sd,
            "standardized_bias": abs(bias) / empirical_sd
                if math.isfinite(empirical_sd) and empirical_sd > 0.0 else math.nan,
            "rmse": math.sqrt(mean([error * error for error in errors])),
            "mean_se": mean(list(acc["ses"])),
            "n_intervals": interval_n,
            "coverage": int(acc["cover"]) / interval_n if interval_n else math.nan,
            "mean_width": mean(list(acc["widths"])),
        })
    return output


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise SystemExit(f"No rows available for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    if args.require_complete:
        validate_complete_run(args.input, args.config)
    paths = sorted(args.input.glob("*/*/rep_*/task_metrics.tsv"))
    if not paths:
        raise SystemExit(f"No task_metrics.tsv files under {args.input}")

    replicate_rows: list[dict[str, object]] = []
    for path in paths:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        missing = {
            str(spec[key]) for spec in METHODS.values() for key in ("q", "alpha", "status")
            if str(spec[key]) not in rows[0]
        }
        if missing:
            raise SystemExit(f"{path}: missing FCR columns: {', '.join(sorted(missing))}")
        replicate_rows.extend(replicate_metrics(rows))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_table(args.output_dir / "effect_fcr_replicates.tsv", replicate_rows)
    write_table(
        args.output_dir / "effect_fcr_summary.tsv",
        summarize_replicates(replicate_rows),
    )
    write_table(
        args.output_dir / "effect_estimation_summary.tsv",
        summarize_effects(paths),
    )
    print(f"Summarized effect inference from {len(paths)} replicate families")


if __name__ == "__main__":
    main()
