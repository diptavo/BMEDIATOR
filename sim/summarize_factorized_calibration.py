#!/usr/bin/env python3
"""Summarize batched factorized calibration, power, and effect estimation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


def as_float(value: object) -> float:
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
    return finite[middle] if len(finite) % 2 else (finite[middle - 1] + finite[middle]) / 2.0


def quantile(values: list[float], probability: float) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    position = probability * (len(finite) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    weight = position - lower
    return finite[lower] * (1.0 - weight) + finite[upper] * weight


def log_mean_exp(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return math.nan
    maximum = max(finite)
    return maximum + math.log(sum(math.exp(value - maximum) for value in finite)) - math.log(len(finite))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rate(rows: list[dict[str, str]], key: str, threshold: float) -> float:
    values = [as_float(row.get(key)) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    return mean([float(value <= threshold) for value in finite])


def bf_rate(rows: list[dict[str, str]], threshold: float = 10.0) -> float:
    cutoff = math.log(threshold)
    pairs = [
        (as_float(row.get("factor_log_BF_XM")), as_float(row.get("factor_log_BF_MY")))
        for row in rows
    ]
    valid = [(left, right) for left, right in pairs if math.isfinite(left) and math.isfinite(right)]
    return mean([float(left >= cutoff and right >= cutoff) for left, right in valid])


def effect_metrics(rows: list[dict[str, str]], estimate: str, standard_error: str, truth: str) -> tuple[float, float, float]:
    valid = []
    for row in rows:
        est = as_float(row.get(estimate))
        se = as_float(row.get(standard_error))
        true = as_float(row.get(truth))
        if math.isfinite(est) and math.isfinite(se) and math.isfinite(true) and se > 0.0:
            valid.append((est, se, true))
    if not valid:
        return math.nan, math.nan, math.nan
    bias = mean([est - true for est, _, true in valid])
    rmse = math.sqrt(mean([(est - true) ** 2 for est, _, true in valid]))
    coverage = mean([float(abs(est - true) <= 1.959963984540054 * se) for est, se, true in valid])
    return bias, rmse, coverage


def interval_coverage(
    rows: list[dict[str, str]], lower: str, upper: str, truth: str
) -> float:
    covered = []
    for row in rows:
        lo = as_float(row.get(lower))
        hi = as_float(row.get(upper))
        true = as_float(row.get(truth))
        if all(math.isfinite(value) for value in (lo, hi, true)):
            covered.append(float(lo <= true <= hi))
    return mean(covered)


def write_scenario_summary(path: Path, rows: list[dict[str, str]]) -> None:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["benchmark"], row["cell"], row["true_scenario"])].append(row)

    fields = [
        "benchmark", "cell", "true_scenario", "identification_class", "n",
        "finite_two_leg_evidence_rate", "stable_two_leg_effect_rate",
        "median_p_XM", "median_p_MY", "median_conjunction_p",
        "median_p_XM_strict", "median_p_MY_strict", "median_strict_conjunction_p",
        "median_PP_XM", "median_PP_MY", "median_PP_two_stage",
        "median_log_BF_directional_XM", "median_log_BF_directional_MY",
        "median_PP_directional_XM", "median_PP_directional_MY",
        "median_directional_collinearity_XM", "median_directional_collinearity_MY",
        "median_log_e_XM", "median_log_e_MY", "median_log_e_mediation",
        "log_mean_e_XM", "log_mean_e_MY", "log_mean_e_mediation",
        "reject_p_XM_0.05", "reject_p_MY_0.05", "reject_conjunction_0.01",
        "reject_conjunction_0.05", "reject_conjunction_0.10", "selected_BY_0.05",
        "selected_EBH_0.05", "selected_strict_BY_0.05", "selected_p2e_EBH_0.05",
        "selected_posterior_FDR_0.05",
        "both_BF_gt_10", "any_directional_BF_gt_10",
        "beta1_bias", "beta1_rmse", "beta1_coverage95",
        "beta1_ci_coverage95", "beta2_bias", "beta2_rmse", "beta2_coverage95",
        "beta2_ci_coverage95",
        "beta3_bias", "beta3_rmse", "beta3_coverage95",
        "directional_intercept_my_bias", "directional_intercept_my_rmse",
        "directional_intercept_my_coverage95",
        "indirect_bias", "indirect_rmse", "indirect_coverage95",
        "indirect_ci_coverage95",
        "mean_nA", "mean_nB",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (benchmark, cell, scenario), group in sorted(grouped.items()):
            beta1_bias, beta1_rmse, beta1_coverage = effect_metrics(
                group, "factor_beta1", "factor_beta1_se", "true_beta1"
            )
            beta2_bias, beta2_rmse, beta2_coverage = effect_metrics(
                group, "factor_beta2", "factor_beta2_se", "true_beta2"
            )
            beta3_bias, beta3_rmse, beta3_coverage = effect_metrics(
                group, "factor_beta3", "factor_beta3_se", "true_beta3"
            )
            indirect_bias, indirect_rmse, indirect_coverage = effect_metrics(
                group, "factor_indirect", "factor_indirect_se", "true_mediated_effect"
            )
            eta_my_bias, eta_my_rmse, eta_my_coverage = effect_metrics(
                group, "factor_directional_intercept_MY",
                "factor_directional_intercept_MY_se",
                "true_directional_intercept_my"
            )
            writer.writerow({
                "benchmark": benchmark,
                "cell": cell,
                "true_scenario": scenario,
                "identification_class": group[0].get("identification_class", ""),
                "n": len(group),
                "finite_two_leg_evidence_rate": mean([
                    float(all(math.isfinite(as_float(row.get(name))) for name in (
                        "factor_p_XM", "factor_p_MY", "factor_log_BF_XM",
                        "factor_log_BF_MY", "factor_log_e_XM", "factor_log_e_MY",
                        "factor_PP_XM", "factor_PP_MY"
                    ))) for row in group
                ]),
                "stable_two_leg_effect_rate": mean([
                    float(all(math.isfinite(as_float(row.get(name))) for name in (
                        "factor_beta1", "factor_beta1_se", "factor_beta2",
                        "factor_beta2_se", "factor_indirect", "factor_indirect_se"
                    ))) for row in group
                ]),
                "median_p_XM": median([as_float(row.get("factor_p_XM")) for row in group]),
                "median_p_MY": median([as_float(row.get("factor_p_MY")) for row in group]),
                "median_conjunction_p": median([as_float(row.get("factor_conjunction_p")) for row in group]),
                "median_p_XM_strict": median([as_float(row.get("factor_p_XM_strict")) for row in group]),
                "median_p_MY_strict": median([as_float(row.get("factor_p_MY_strict")) for row in group]),
                "median_strict_conjunction_p": median([as_float(row.get("factor_strict_conjunction_p")) for row in group]),
                "median_PP_XM": median([as_float(row.get("factor_PP_XM")) for row in group]),
                "median_PP_MY": median([as_float(row.get("factor_PP_MY")) for row in group]),
                "median_PP_two_stage": median([as_float(row.get("factor_PP_two_stage")) for row in group]),
                "median_log_BF_directional_XM": median([as_float(row.get("factor_log_BF_directional_XM")) for row in group]),
                "median_log_BF_directional_MY": median([as_float(row.get("factor_log_BF_directional_MY")) for row in group]),
                "median_PP_directional_XM": median([as_float(row.get("factor_PP_directional_XM")) for row in group]),
                "median_PP_directional_MY": median([as_float(row.get("factor_PP_directional_MY")) for row in group]),
                "median_directional_collinearity_XM": median([as_float(row.get("factor_directional_collinearity_XM")) for row in group]),
                "median_directional_collinearity_MY": median([as_float(row.get("factor_directional_collinearity_MY")) for row in group]),
                "median_log_e_XM": median([as_float(row.get("factor_log_e_XM")) for row in group]),
                "median_log_e_MY": median([as_float(row.get("factor_log_e_MY")) for row in group]),
                "median_log_e_mediation": median([as_float(row.get("factor_log_e_mediation")) for row in group]),
                "log_mean_e_XM": log_mean_exp([as_float(row.get("factor_log_e_XM")) for row in group]),
                "log_mean_e_MY": log_mean_exp([as_float(row.get("factor_log_e_MY")) for row in group]),
                "log_mean_e_mediation": log_mean_exp([as_float(row.get("factor_log_e_mediation")) for row in group]),
                "reject_p_XM_0.05": rate(group, "factor_p_XM", 0.05),
                "reject_p_MY_0.05": rate(group, "factor_p_MY", 0.05),
                "reject_conjunction_0.01": rate(group, "factor_conjunction_p", 0.01),
                "reject_conjunction_0.05": rate(group, "factor_conjunction_p", 0.05),
                "reject_conjunction_0.10": rate(group, "factor_conjunction_p", 0.10),
                "selected_BY_0.05": rate(group, "factor_conjunction_q_BY", 0.05),
                "selected_EBH_0.05": rate(group, "factor_e_q_EBH", 0.05),
                "selected_strict_BY_0.05": rate(group, "factor_strict_conjunction_q_BY", 0.05),
                "selected_p2e_EBH_0.05": rate(group, "factor_e_q_p2e_EBH", 0.05),
                "selected_posterior_FDR_0.05": rate(group, "factor_posterior_cum_fdr", 0.05),
                "both_BF_gt_10": bf_rate(group),
                "any_directional_BF_gt_10": mean([
                    float(max(as_float(row.get("factor_log_BF_directional_XM")),
                              as_float(row.get("factor_log_BF_directional_MY"))) >= math.log(10.0))
                    for row in group
                    if math.isfinite(as_float(row.get("factor_log_BF_directional_XM"))) and
                       math.isfinite(as_float(row.get("factor_log_BF_directional_MY")))
                ]),
                "beta1_bias": beta1_bias,
                "beta1_rmse": beta1_rmse,
                "beta1_coverage95": beta1_coverage,
                "beta1_ci_coverage95": interval_coverage(
                    group, "factor_beta1_ci_lower", "factor_beta1_ci_upper", "true_beta1"
                ),
                "beta2_bias": beta2_bias,
                "beta2_rmse": beta2_rmse,
                "beta2_coverage95": beta2_coverage,
                "beta2_ci_coverage95": interval_coverage(
                    group, "factor_beta2_ci_lower", "factor_beta2_ci_upper", "true_beta2"
                ),
                "beta3_bias": beta3_bias,
                "beta3_rmse": beta3_rmse,
                "beta3_coverage95": beta3_coverage,
                "directional_intercept_my_bias": eta_my_bias,
                "directional_intercept_my_rmse": eta_my_rmse,
                "directional_intercept_my_coverage95": eta_my_coverage,
                "indirect_bias": indirect_bias,
                "indirect_rmse": indirect_rmse,
                "indirect_coverage95": indirect_coverage,
                "indirect_ci_coverage95": interval_coverage(
                    group, "factor_indirect_ci_lower", "factor_indirect_ci_upper",
                    "true_mediated_effect"
                ),
                "mean_nA": mean([as_float(row.get("factor_nA")) for row in group]),
                "mean_nB": mean([as_float(row.get("factor_nB")) for row in group]),
            })


def calculate_replicate_decisions(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["benchmark"], row["cell"], row["replicate"])].append(row)
    output: list[dict[str, object]] = []
    for (benchmark, cell, replicate), group in sorted(grouped.items()):
        identifiable = [row for row in group if row.get("identification_class") != "nonidentifiable"]
        alternatives = [row for row in identifiable if row["true_scenario"] == "M1"]
        boundary = [row for row in group if row.get("identification_class") == "nonidentifiable"]
        for rule, key in (
            ("raw_0.05", "factor_conjunction_p"),
            ("BY_0.05", "factor_conjunction_q_BY"),
            ("eBH_0.05", "factor_e_q_EBH"),
            ("strict_raw_0.05", "factor_strict_conjunction_p"),
            ("strict_BY_0.05", "factor_strict_conjunction_q_BY"),
            ("p2e_eBH_0.05", "factor_e_q_p2e_EBH"),
            ("posterior_FDR_0.05", "factor_posterior_cum_fdr"),
        ):
            tested = [row for row in identifiable if math.isfinite(as_float(row.get(key)))]
            tested_alternatives = [row for row in tested if row["true_scenario"] == "M1"]
            selected = [row for row in tested if as_float(row.get(key)) <= 0.05]
            true_positive = sum(row["true_scenario"] == "M1" for row in selected)
            false_positive = len(selected) - true_positive
            boundary_valid = [row for row in boundary if math.isfinite(as_float(row.get(key)))]
            posterior_expected_fdp = math.nan
            if rule == "posterior_FDR_0.05" and selected:
                posterior_expected_fdp = mean([
                    as_float(row.get("factor_posterior_local_fdr")) for row in selected
                ])
            empirical_fdr = (
                false_positive / len(selected) if selected
                else (0.0 if tested else math.nan)
            )
            output.append({
                "benchmark": benchmark,
                "cell": cell,
                "replicate": replicate,
                "rule": rule,
                "n_identifiable": len(identifiable),
                "n_tested": len(tested),
                "n_true_mediation": len(alternatives),
                "n_tested_true_mediation": len(tested_alternatives),
                "n_selected": len(selected),
                "true_positive": true_positive,
                "false_positive": false_positive,
                "empirical_fdr": empirical_fdr,
                "posterior_expected_fdp": posterior_expected_fdp,
                "posterior_fdp_gap": (
                    empirical_fdr - posterior_expected_fdp
                    if math.isfinite(empirical_fdr) and
                    math.isfinite(posterior_expected_fdp) else math.nan
                ),
                "power": (
                    true_positive / len(tested_alternatives)
                    if tested_alternatives else math.nan
                ),
                "n_nonidentifiable_boundary": len(boundary),
                "boundary_selection_rate": mean([
                    float(as_float(row.get(key)) <= 0.05) for row in boundary_valid
                ]),
            })
    return output


def write_replicate_decisions(path: Path, rows: list[dict[str, str]]) -> list[dict[str, object]]:
    output = calculate_replicate_decisions(rows)
    fields = list(output[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    return output


def write_decision_summary(path: Path, rows: list[dict[str, object]]) -> None:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["benchmark"]), str(row["cell"]), str(row["rule"]))].append(row)
    fields = [
        "benchmark", "cell", "rule", "replicates", "mean_tested",
        "mean_selected", "mean_fdr", "fdr_q95", "mean_power", "power_q05",
        "mean_boundary_selection_rate", "mean_posterior_expected_fdp",
        "mean_posterior_fdp_gap",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (benchmark, cell, rule), group in sorted(grouped.items()):
            fdr = [as_float(row["empirical_fdr"]) for row in group]
            power = [as_float(row["power"]) for row in group]
            boundary = [as_float(row["boundary_selection_rate"]) for row in group]
            writer.writerow({
                "benchmark": benchmark,
                "cell": cell,
                "rule": rule,
                "replicates": len(group),
                "mean_tested": mean([as_float(row["n_tested"]) for row in group]),
                "mean_selected": mean([as_float(row["n_selected"]) for row in group]),
                "mean_fdr": mean(fdr),
                "fdr_q95": quantile(fdr, 0.95),
                "mean_power": mean(power),
                "power_q05": quantile(power, 0.05),
                "mean_boundary_selection_rate": mean(boundary),
                "mean_posterior_expected_fdp": mean([
                    as_float(row["posterior_expected_fdp"]) for row in group
                ]),
                "mean_posterior_fdp_gap": mean([
                    as_float(row["posterior_fdp_gap"]) for row in group
                ]),
            })


def update_posterior_calibration(
    grouped: dict[tuple[str, str, str, int], dict[str, float]],
    row: dict[str, str],
) -> None:
    posterior = as_float(row.get("factor_PP_two_stage"))
    if not math.isfinite(posterior):
        return
    group = (
        "nonidentifiable"
        if row.get("identification_class") == "nonidentifiable"
        else "identifiable"
    )
    bin_index = min(9, max(0, int(posterior * 10.0)))
    key = (row["benchmark"], row["cell"], group, bin_index)
    if key not in grouped:
        grouped[key] = {"n": 0.0, "sum_posterior": 0.0, "true": 0.0, "brier": 0.0}
    truth = float(row.get("true_scenario") == "M1")
    acc = grouped[key]
    acc["n"] += 1.0
    acc["sum_posterior"] += posterior
    acc["true"] += truth
    acc["brier"] += (posterior - truth) ** 2


def write_posterior_calibration(
    path: Path,
    grouped: dict[tuple[str, str, str, int], dict[str, float]],
) -> None:
    fields = [
        "benchmark", "cell", "identification_group", "bin_lower", "bin_upper",
        "n", "mean_posterior", "observed_mediation_rate", "calibration_gap",
        "mean_brier",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for (benchmark, cell, group, bin_index), acc in sorted(grouped.items()):
            count = acc["n"]
            predicted = acc["sum_posterior"] / count
            observed = acc["true"] / count
            writer.writerow({
                "benchmark": benchmark,
                "cell": cell,
                "identification_group": group,
                "bin_lower": bin_index / 10.0,
                "bin_upper": (bin_index + 1) / 10.0,
                "n": int(count),
                "mean_posterior": predicted,
                "observed_mediation_rate": observed,
                "calibration_gap": observed - predicted,
                "mean_brier": acc["brier"] / count,
            })


def validate_complete_run(input_dir: Path, config_path: Path) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config_sha256 = sha256_file(config_path)
    expected: dict[tuple[str, str, str], int] = {}
    for benchmark in ("classification", "calibration"):
        if benchmark not in config:
            continue
        section = config[benchmark]
        if benchmark == "classification":
            rows_per_task = int(section["proteins_per_scenario"]) * len(section["scenarios"])
        else:
            rows_per_task = int(section["proteins_per_replicate"])
        for cell in section["cells"]:
            replicates = int(cell.get("replicates", section["replicates"]))
            for replicate in range(1, replicates + 1):
                expected[(benchmark, str(cell["name"]), f"rep_{replicate:04d}")] = rows_per_task

    observed: dict[tuple[str, str, str], Path] = {}
    binary_hashes: set[str] = set()
    errors: list[str] = []
    for path in sorted(input_dir.glob("*/*/rep_*/task_metrics.tsv")):
        relative = path.relative_to(input_dir)
        key = (relative.parts[0], relative.parts[1], relative.parts[2])
        if key in observed:
            errors.append(f"duplicate task key {key}: {observed[key]} and {path}")
            continue
        observed[key] = path
        with path.open(newline="") as handle:
            task_rows = list(csv.DictReader(handle, delimiter="\t"))
        expected_rows = expected.get(key)
        if expected_rows is None:
            errors.append(f"unexpected task output: {path}")
        elif len(task_rows) != expected_rows:
            errors.append(f"{path}: expected {expected_rows} rows, found {len(task_rows)}")
        protein_ids = [str(row.get("protein_id", "")) for row in task_rows]
        if len(protein_ids) != len(set(protein_ids)):
            errors.append(f"{path}: duplicate protein_id values")
        metadata_path = path.with_name("task_metadata.json")
        if not metadata_path.is_file():
            errors.append(f"missing task metadata: {metadata_path}")
        else:
            metadata = json.loads(metadata_path.read_text(encoding="ascii"))
            if metadata.get("config_sha256") != config_sha256:
                errors.append(f"{metadata_path}: config hash does not match {config_path}")
            if int(metadata.get("proteins", -1)) != len(task_rows):
                errors.append(f"{metadata_path}: protein count does not match task output")
            identity = (
                metadata.get("benchmark"),
                metadata.get("cell"),
                f"rep_{int(metadata.get('replicate', -1)):04d}",
            )
            if identity != key:
                errors.append(f"{metadata_path}: task identity does not match its path")
            binary_hash = str(metadata.get("binary_sha256", ""))
            if not binary_hash:
                errors.append(f"{metadata_path}: missing binary hash")
            else:
                binary_hashes.add(binary_hash)

    missing = sorted(set(expected) - set(observed))
    if missing:
        preview = ", ".join("/".join(key) for key in missing[:10])
        errors.append(f"missing {len(missing)} task outputs; first: {preview}")
    if len(binary_hashes) > 1:
        errors.append(f"tasks were generated by {len(binary_hashes)} different binaries")
    if errors:
        raise SystemExit("Incomplete or inconsistent calibration run:\n- " + "\n- ".join(errors))


def update_scenario_accumulator(
    grouped: dict[tuple[str, str, str], dict[str, object]],
    row: dict[str, str],
) -> None:
    key = (row["benchmark"], row["cell"], row["true_scenario"])
    if key not in grouped:
        grouped[key] = {
            "identification_class": row.get("identification_class", ""),
            "n": 0,
            "values": defaultdict(list),
            "rates": defaultdict(lambda: [0, 0]),
            "bf": [0, 0],
            "directional_bf": [0, 0],
            "validity": {"evidence": 0, "effect": 0},
            "effects": {
                "beta1": [0, 0.0, 0.0, 0, 0, 0],
                "beta2": [0, 0.0, 0.0, 0, 0, 0],
                "beta3": [0, 0.0, 0.0, 0, 0, 0],
                "directional_intercept_my": [0, 0.0, 0.0, 0, 0, 0],
                "indirect": [0, 0.0, 0.0, 0, 0, 0],
            },
            "counts": {"nA": [0, 0.0], "nB": [0, 0.0]},
        }
    acc = grouped[key]
    acc["n"] = int(acc["n"]) + 1
    acc["validity"]["evidence"] += int(all(
        math.isfinite(as_float(row.get(name))) for name in (
            "factor_p_XM", "factor_p_MY", "factor_log_BF_XM",
            "factor_log_BF_MY", "factor_log_e_XM", "factor_log_e_MY",
            "factor_PP_XM", "factor_PP_MY"
        )
    ))
    acc["validity"]["effect"] += int(all(
        math.isfinite(as_float(row.get(name))) for name in (
            "factor_beta1", "factor_beta1_se", "factor_beta2",
            "factor_beta2_se", "factor_indirect", "factor_indirect_se"
        )
    ))
    values = acc["values"]
    for output_name, input_name in (
        ("p_XM", "factor_p_XM"),
        ("p_MY", "factor_p_MY"),
        ("conjunction_p", "factor_conjunction_p"),
        ("p_XM_strict", "factor_p_XM_strict"),
        ("p_MY_strict", "factor_p_MY_strict"),
        ("strict_conjunction_p", "factor_strict_conjunction_p"),
        ("PP_XM", "factor_PP_XM"),
        ("PP_MY", "factor_PP_MY"),
        ("PP_two_stage", "factor_PP_two_stage"),
        ("log_BF_directional_XM", "factor_log_BF_directional_XM"),
        ("log_BF_directional_MY", "factor_log_BF_directional_MY"),
        ("PP_directional_XM", "factor_PP_directional_XM"),
        ("PP_directional_MY", "factor_PP_directional_MY"),
        ("directional_collinearity_XM", "factor_directional_collinearity_XM"),
        ("directional_collinearity_MY", "factor_directional_collinearity_MY"),
        ("log_e_XM", "factor_log_e_XM"),
        ("log_e_MY", "factor_log_e_MY"),
        ("log_e_mediation", "factor_log_e_mediation"),
    ):
        value = as_float(row.get(input_name))
        if math.isfinite(value):
            values[output_name].append(value)

    rates = acc["rates"]
    for name, input_name, threshold in (
        ("reject_p_XM_0.05", "factor_p_XM", 0.05),
        ("reject_p_MY_0.05", "factor_p_MY", 0.05),
        ("reject_conjunction_0.01", "factor_conjunction_p", 0.01),
        ("reject_conjunction_0.05", "factor_conjunction_p", 0.05),
        ("reject_conjunction_0.10", "factor_conjunction_p", 0.10),
        ("selected_BY_0.05", "factor_conjunction_q_BY", 0.05),
        ("selected_EBH_0.05", "factor_e_q_EBH", 0.05),
        ("selected_strict_BY_0.05", "factor_strict_conjunction_q_BY", 0.05),
        ("selected_p2e_EBH_0.05", "factor_e_q_p2e_EBH", 0.05),
        ("selected_posterior_FDR_0.05", "factor_posterior_cum_fdr", 0.05),
    ):
        value = as_float(row.get(input_name))
        if math.isfinite(value):
            rates[name][1] += 1
            rates[name][0] += int(value <= threshold)

    left = as_float(row.get("factor_log_BF_XM"))
    right = as_float(row.get("factor_log_BF_MY"))
    if math.isfinite(left) and math.isfinite(right):
        acc["bf"][1] += 1
        acc["bf"][0] += int(left >= math.log(10.0) and right >= math.log(10.0))
    directional_xm = as_float(row.get("factor_log_BF_directional_XM"))
    directional_my = as_float(row.get("factor_log_BF_directional_MY"))
    if math.isfinite(directional_xm) and math.isfinite(directional_my):
        acc["directional_bf"][1] += 1
        acc["directional_bf"][0] += int(
            max(directional_xm, directional_my) >= math.log(10.0)
        )

    for label, estimate, standard_error, truth, lower, upper in (
        ("beta1", "factor_beta1", "factor_beta1_se", "true_beta1",
         "factor_beta1_ci_lower", "factor_beta1_ci_upper"),
        ("beta2", "factor_beta2", "factor_beta2_se", "true_beta2",
         "factor_beta2_ci_lower", "factor_beta2_ci_upper"),
        ("beta3", "factor_beta3", "factor_beta3_se", "true_beta3", "", ""),
        ("directional_intercept_my", "factor_directional_intercept_MY",
         "factor_directional_intercept_MY_se", "true_directional_intercept_my",
         "", ""),
        ("indirect", "factor_indirect", "factor_indirect_se", "true_mediated_effect",
         "factor_indirect_ci_lower", "factor_indirect_ci_upper"),
    ):
        est = as_float(row.get(estimate))
        se = as_float(row.get(standard_error))
        true = as_float(row.get(truth))
        if math.isfinite(est) and math.isfinite(se) and math.isfinite(true) and se > 0.0:
            error = est - true
            effect = acc["effects"][label]
            effect[0] += 1
            effect[1] += error
            effect[2] += error * error
            effect[3] += int(abs(error) <= 1.959963984540054 * se)
            lo = as_float(row.get(lower)) if lower else math.nan
            hi = as_float(row.get(upper)) if upper else math.nan
            if math.isfinite(lo) and math.isfinite(hi):
                effect[4] += 1
                effect[5] += int(lo <= true <= hi)

    for label, input_name in (("nA", "factor_nA"), ("nB", "factor_nB")):
        value = as_float(row.get(input_name))
        if math.isfinite(value):
            acc["counts"][label][0] += 1
            acc["counts"][label][1] += value


def write_streaming_scenario_summary(
    path: Path,
    grouped: dict[tuple[str, str, str], dict[str, object]],
) -> None:
    fields = [
        "benchmark", "cell", "true_scenario", "identification_class", "n",
        "finite_two_leg_evidence_rate", "stable_two_leg_effect_rate",
        "median_p_XM", "median_p_MY", "median_conjunction_p",
        "median_p_XM_strict", "median_p_MY_strict", "median_strict_conjunction_p",
        "median_PP_XM", "median_PP_MY", "median_PP_two_stage",
        "median_log_BF_directional_XM", "median_log_BF_directional_MY",
        "median_PP_directional_XM", "median_PP_directional_MY",
        "median_directional_collinearity_XM", "median_directional_collinearity_MY",
        "median_log_e_XM", "median_log_e_MY", "median_log_e_mediation",
        "log_mean_e_XM", "log_mean_e_MY", "log_mean_e_mediation",
        "reject_p_XM_0.05", "reject_p_MY_0.05", "reject_conjunction_0.01",
        "reject_conjunction_0.05", "reject_conjunction_0.10", "selected_BY_0.05",
        "selected_EBH_0.05", "selected_strict_BY_0.05", "selected_p2e_EBH_0.05",
        "selected_posterior_FDR_0.05",
        "both_BF_gt_10", "any_directional_BF_gt_10", "beta1_bias", "beta1_rmse",
        "beta1_coverage95", "beta1_ci_coverage95",
        "beta2_bias", "beta2_rmse", "beta2_coverage95", "beta2_ci_coverage95",
        "beta3_bias", "beta3_rmse", "beta3_coverage95",
        "directional_intercept_my_bias", "directional_intercept_my_rmse",
        "directional_intercept_my_coverage95",
        "indirect_bias", "indirect_rmse", "indirect_coverage95",
        "indirect_ci_coverage95",
        "mean_nA", "mean_nB",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (benchmark, cell, scenario), acc in sorted(grouped.items()):
            values = acc["values"]
            rates = acc["rates"]
            effects = {}
            for label in ("beta1", "beta2", "beta3", "directional_intercept_my", "indirect"):
                count, bias_sum, square_sum, coverage_count, ci_count, ci_covered = \
                    acc["effects"][label]
                effects[label] = (
                    bias_sum / count if count else math.nan,
                    math.sqrt(square_sum / count) if count else math.nan,
                    coverage_count / count if count else math.nan,
                    ci_covered / ci_count if ci_count else math.nan,
                )
            rate_values = {}
            for name in (
                "reject_p_XM_0.05", "reject_p_MY_0.05",
                "reject_conjunction_0.01", "reject_conjunction_0.05",
                "reject_conjunction_0.10", "selected_BY_0.05",
                "selected_EBH_0.05",
                "selected_strict_BY_0.05", "selected_p2e_EBH_0.05",
                "selected_posterior_FDR_0.05",
            ):
                successes, total = rates[name]
                rate_values[name] = successes / total if total else math.nan
            bf_successes, bf_total = acc["bf"]
            directional_successes, directional_total = acc["directional_bf"]
            n_a_count, n_a_sum = acc["counts"]["nA"]
            n_b_count, n_b_sum = acc["counts"]["nB"]
            writer.writerow({
                "benchmark": benchmark,
                "cell": cell,
                "true_scenario": scenario,
                "identification_class": acc["identification_class"],
                "n": acc["n"],
                "finite_two_leg_evidence_rate": acc["validity"]["evidence"] / acc["n"],
                "stable_two_leg_effect_rate": acc["validity"]["effect"] / acc["n"],
                "median_p_XM": median(values["p_XM"]),
                "median_p_MY": median(values["p_MY"]),
                "median_conjunction_p": median(values["conjunction_p"]),
                "median_p_XM_strict": median(values["p_XM_strict"]),
                "median_p_MY_strict": median(values["p_MY_strict"]),
                "median_strict_conjunction_p": median(values["strict_conjunction_p"]),
                "median_PP_XM": median(values["PP_XM"]),
                "median_PP_MY": median(values["PP_MY"]),
                "median_PP_two_stage": median(values["PP_two_stage"]),
                "median_log_BF_directional_XM": median(values["log_BF_directional_XM"]),
                "median_log_BF_directional_MY": median(values["log_BF_directional_MY"]),
                "median_PP_directional_XM": median(values["PP_directional_XM"]),
                "median_PP_directional_MY": median(values["PP_directional_MY"]),
                "median_directional_collinearity_XM": median(values["directional_collinearity_XM"]),
                "median_directional_collinearity_MY": median(values["directional_collinearity_MY"]),
                "median_log_e_XM": median(values["log_e_XM"]),
                "median_log_e_MY": median(values["log_e_MY"]),
                "median_log_e_mediation": median(values["log_e_mediation"]),
                "log_mean_e_XM": log_mean_exp(values["log_e_XM"]),
                "log_mean_e_MY": log_mean_exp(values["log_e_MY"]),
                "log_mean_e_mediation": log_mean_exp(values["log_e_mediation"]),
                **rate_values,
                "both_BF_gt_10": bf_successes / bf_total if bf_total else math.nan,
                "any_directional_BF_gt_10": directional_successes / directional_total
                    if directional_total else math.nan,
                "beta1_bias": effects["beta1"][0],
                "beta1_rmse": effects["beta1"][1],
                "beta1_coverage95": effects["beta1"][2],
                "beta1_ci_coverage95": effects["beta1"][3],
                "beta2_bias": effects["beta2"][0],
                "beta2_rmse": effects["beta2"][1],
                "beta2_coverage95": effects["beta2"][2],
                "beta2_ci_coverage95": effects["beta2"][3],
                "beta3_bias": effects["beta3"][0],
                "beta3_rmse": effects["beta3"][1],
                "beta3_coverage95": effects["beta3"][2],
                "directional_intercept_my_bias": effects["directional_intercept_my"][0],
                "directional_intercept_my_rmse": effects["directional_intercept_my"][1],
                "directional_intercept_my_coverage95": effects["directional_intercept_my"][2],
                "indirect_bias": effects["indirect"][0],
                "indirect_rmse": effects["indirect"][1],
                "indirect_coverage95": effects["indirect"][2],
                "indirect_ci_coverage95": effects["indirect"][3],
                "mean_nA": n_a_sum / n_a_count if n_a_count else math.nan,
                "mean_nB": n_b_sum / n_b_count if n_b_count else math.nan,
            })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    if args.require_complete and args.config is None:
        raise SystemExit("--require-complete requires --config")
    if args.config is not None:
        validate_complete_run(args.input, args.config)

    paths = sorted(args.input.glob("*/*/rep_*/task_metrics.tsv"))
    if not paths:
        raise SystemExit(f"No task_metrics.tsv files under {args.input}")
    scenario_groups: dict[tuple[str, str, str], dict[str, object]] = {}
    replicate_rows: list[dict[str, object]] = []
    posterior_calibration: dict[tuple[str, str, str, int], dict[str, float]] = {}
    total_rows = 0
    for path in paths:
        with path.open(newline="") as handle:
            task_rows = list(csv.DictReader(handle, delimiter="\t"))
        total_rows += len(task_rows)
        for row in task_rows:
            update_scenario_accumulator(scenario_groups, row)
            update_posterior_calibration(posterior_calibration, row)
        replicate_rows.extend(calculate_replicate_decisions(task_rows))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_streaming_scenario_summary(
        args.output_dir / "factorized_scenario_summary.tsv", scenario_groups
    )
    replicate_path = args.output_dir / "factorized_replicate_decisions.tsv"
    with replicate_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(replicate_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(replicate_rows)
    write_decision_summary(args.output_dir / "factorized_decision_summary.tsv", replicate_rows)
    write_posterior_calibration(
        args.output_dir / "factorized_posterior_calibration.tsv",
        posterior_calibration,
    )
    print(f"Summarized {total_rows} proteins from {len(paths)} batched replicates")


if __name__ == "__main__":
    main()
