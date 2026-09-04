#!/usr/bin/env python3
"""Deterministic checks for the experimental analytical calibration tracks."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


G_GRID = (0.25, 1.0, 4.0, 16.0)
KAPPA_GRID = (0.10, 0.25, 0.50, 0.75)


def adaptive_e(z: float) -> float:
    return sum(
        math.exp(-0.5 * math.log1p(g) + 0.5 * g * z * z / (1.0 + g))
        for g in G_GRID
    ) / len(G_GRID)


def p_to_e(p_value: float) -> float:
    return sum(
        kappa * p_value ** (kappa - 1.0) for kappa in KAPPA_GRID
    ) / len(KAPPA_GRID)


def ebh_q(rows: list[dict[str, str]], log_e_column: str) -> dict[str, float]:
    eligible = [row for row in rows if finite(row[log_e_column])]
    eligible.sort(key=lambda row: float(row[log_e_column]), reverse=True)
    raw = [
        len(eligible) / (rank * math.exp(float(row[log_e_column])))
        for rank, row in enumerate(eligible, start=1)
    ]
    adjusted = [1.0] * len(raw)
    running = math.inf
    for index in range(len(raw) - 1, -1, -1):
        running = min(running, raw[index])
        adjusted[index] = min(1.0, running)
    return {row["Protein"]: value for row, value in zip(eligible, adjusted)}


def bh_q(rows: list[dict[str, str]], p_column: str) -> dict[str, float]:
    eligible = [row for row in rows if finite(row[p_column])]
    eligible.sort(key=lambda row: float(row[p_column]))
    raw = [
        float(row[p_column]) * len(eligible) / rank
        for rank, row in enumerate(eligible, start=1)
    ]
    adjusted = [1.0] * len(raw)
    running = 1.0
    for index in range(len(raw) - 1, -1, -1):
        running = min(running, raw[index])
        adjusted[index] = min(1.0, running)
    return {row["Protein"]: value for row, value in zip(eligible, adjusted)}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def close(left: float, right: float, tolerance: float = 2e-6) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def adafilter_q(rows: list[dict[str, str]]) -> dict[str, float]:
    eligible = [
        row
        for row in rows
        if finite(row["factor_p_XM_balanced"])
        and finite(row["factor_p_MY_balanced"])
    ]
    eligible.sort(key=lambda row: float(row["factor_balanced_conjunction_p"]))
    raw = []
    for rank, row in enumerate(eligible, start=1):
        selection = float(row["factor_balanced_conjunction_p"])
        adjustment = sum(
            min(float(other["factor_p_XM_balanced"]),
                float(other["factor_p_MY_balanced"]))
            <= selection
            for other in eligible
        )
        raw.append(selection * adjustment / rank)
    adjusted = [1.0] * len(raw)
    running = 1.0
    for index in range(len(raw) - 1, -1, -1):
        running = min(running, raw[index])
        adjusted[index] = running
    return {row["Protein"]: value for row, value in zip(eligible, adjusted)}


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_analytic_calibration.py MAIN.mediation OVERLAP.mediation")

    # Integrate phi(z) E(z). Every component is a normalized N(0, 1+g)
    # density ratio, so the expectation must be one. The wide deterministic
    # grid also catches sign and normalizing-constant errors in the formula.
    step = 0.002
    total = 0.0
    for index in range(40001):
        z = -40.0 + index * step
        weight = 0.5 if index in (0, 40000) else 1.0
        null_weighted_e = sum(
            math.exp(
                -0.5 * math.log(2.0 * math.pi)
                -0.5 * math.log1p(g)
                -0.5 * z * z / (1.0 + g)
            )
            for g in G_GRID
        ) / len(G_GRID)
        total += weight * null_weighted_e
    total *= step
    if not close(total, 1.0, 2e-5):
        raise SystemExit(f"adaptive e-value null expectation is {total}, expected 1")

    rows = read_rows(Path(sys.argv[1]))
    expected_balanced_bh = bh_q(rows, "factor_balanced_conjunction_p")
    expected_adafilter = adafilter_q(rows)
    expected_balanced_t_ebh = ebh_q(rows, "factor_log_e_mediation_balanced")
    expected_balanced_ebh = ebh_q(rows, "factor_log_e_p2e_balanced_mediation")
    checked = 0
    for row in rows:
        if row["Protein"] not in expected_adafilter:
            continue
        p_xm = float(row["factor_p_XM_balanced"])
        p_my = float(row["factor_p_MY_balanced"])
        conjunction = float(row["factor_balanced_conjunction_p"])
        if not close(conjunction, max(p_xm, p_my)):
            raise SystemExit("balanced conjunction is not max(p_XM, p_MY)")
        reported_bh = float(row["factor_balanced_conjunction_q_BH"])
        if not close(reported_bh, expected_balanced_bh[row["Protein"]]):
            raise SystemExit("balanced BH adjusted value does not match its definition")
        reported = float(row["factor_balanced_conjunction_q_AdaFilter"])
        if not close(reported, expected_adafilter[row["Protein"]]):
            raise SystemExit("AdaFilter adjusted value does not match its definition")
        balanced_log_e = min(
            float(row["factor_log_e_XM_balanced"]),
            float(row["factor_log_e_MY_balanced"]),
        )
        if not close(balanced_log_e, float(row["factor_log_e_mediation_balanced"])):
            raise SystemExit("balanced mediation e-value is not the minimum leg e-value")
        reported_balanced_t_q = float(row["factor_e_q_balanced_EBH"])
        if not close(
            reported_balanced_t_q, expected_balanced_t_ebh[row["Protein"]]
        ):
            raise SystemExit("balanced Student-t e-BH value does not match its definition")
        expected_log_e = math.log(p_to_e(conjunction))
        reported_log_e = float(row["factor_log_e_p2e_balanced_mediation"])
        if not close(reported_log_e, expected_log_e):
            raise SystemExit("balanced p-to-e value does not match its fixed calibrator")
        reported_e_q = float(row["factor_e_q_p2e_balanced_EBH"])
        if not close(reported_e_q, expected_balanced_ebh[row["Protein"]]):
            raise SystemExit("balanced p-to-e e-BH value does not match its definition")
        log_e = min(
            float(row["factor_log_e_XM_adaptive"]),
            float(row["factor_log_e_MY_adaptive"]),
        )
        if not close(log_e, float(row["factor_log_e_mediation_adaptive"])):
            raise SystemExit("mediation e-value is not the minimum leg e-value")
        checked += 1
    if checked == 0:
        raise SystemExit("no eligible analytical calibration rows were checked")

    overlap_rows = read_rows(Path(sys.argv[2]))
    for row in overlap_rows:
        if row["factor_adaptive_ebh_status"] == "UNRESOLVED_SAMPLE_OVERLAP":
            if finite(row["factor_e_q_adaptive_EBH"]):
                raise SystemExit("adaptive e-BH value must fail closed under sample overlap")
        if row["factor_adafilter_status"] == "UNRESOLVED_SAMPLE_OVERLAP":
            if finite(row["factor_balanced_conjunction_q_AdaFilter"]):
                raise SystemExit("AdaFilter value must fail closed under sample overlap")
        if row["factor_balanced_bh_status"] == "UNRESOLVED_SAMPLE_OVERLAP":
            if finite(row["factor_balanced_conjunction_q_BH"]):
                raise SystemExit("balanced BH value must fail closed under sample overlap")
        if row["factor_balanced_p2e_status"] == "UNRESOLVED_SAMPLE_OVERLAP":
            if finite(row["factor_e_q_p2e_balanced_EBH"]):
                raise SystemExit("balanced p-to-e e-BH value must fail closed under sample overlap")
        if row["factor_balanced_ebh_status"] == "UNRESOLVED_SAMPLE_OVERLAP":
            if finite(row["factor_e_q_balanced_EBH"]):
                raise SystemExit("balanced Student-t e-BH value must fail closed under sample overlap")


if __name__ == "__main__":
    main()
