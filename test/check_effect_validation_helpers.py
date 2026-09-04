#!/usr/bin/env python3
"""Deterministic checks for effect-validation bookkeeping."""

from __future__ import annotations

import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "sim"))

from run_factorized_task import factorized_options  # noqa: E402
from summarize_effect_validation import replicate_metrics  # noqa: E402


def base_row(protein: str) -> dict[str, str]:
    row = {
        "benchmark": "calibration",
        "cell": "unit",
        "replicate": "1",
        "protein_id": protein,
        "true_scenario": "M1",
        "true_beta1": "0.3",
        "true_beta2": "0.2",
        "true_mediated_effect": "0.06",
        "factor_balanced_conjunction_p": "0.001",
        "factor_balanced_conjunction_q_BH": "0.01",
        "factor_balanced_conjunction_q_BY": "0.01",
        "factor_fcr_alpha_BH": "0.025",
        "factor_fcr_alpha_BY": "0.005",
        "factor_fcr_bh_status": "PROFILE_FCR_BH_INDEPENDENCE_CONDITIONAL",
        "factor_fcr_by_status": "PROFILE_FCR_BY_DEPENDENCE_CONDITIONAL",
    }
    for suffix in ("BH", "BY"):
        row[f"factor_beta1_fcr_ci_lower_{suffix}"] = "0.1"
        row[f"factor_beta1_fcr_ci_upper_{suffix}"] = "0.5"
        row[f"factor_beta2_fcr_ci_lower_{suffix}"] = "0.1"
        row[f"factor_beta2_fcr_ci_upper_{suffix}"] = "0.3"
        row[f"factor_indirect_fcr_ci_lower_{suffix}"] = "0.01"
        row[f"factor_indirect_fcr_ci_upper_{suffix}"] = "0.15"
    row.update({
        "factor_beta1_ci_lower": "0.2",
        "factor_beta1_ci_upper": "0.4",
        "factor_beta2_ci_lower": "0.1",
        "factor_beta2_ci_upper": "0.3",
        "factor_indirect_ci_lower": "0.02",
        "factor_indirect_ci_upper": "0.12",
    })
    return row


def main() -> None:
    bounded = base_row("bounded")
    unbounded = base_row("unbounded")
    unbounded["factor_fcr_bh_status"] = "UNBOUNDED_FCR_EFFECT_SET"
    unbounded["factor_fcr_by_status"] = "UNBOUNDED_FCR_EFFECT_SET"
    for suffix in ("BH", "BY"):
        for effect in ("beta1", "beta2", "indirect"):
            unbounded[f"factor_{effect}_fcr_ci_lower_{suffix}"] = "-inf"
            unbounded[f"factor_{effect}_fcr_ci_upper_{suffix}"] = "inf"

    missing = base_row("missing")
    missing["factor_beta2_fcr_ci_lower_BH"] = "nan"
    metrics = {row["method"]: row for row in replicate_metrics([bounded, unbounded, missing])}
    bh = metrics["BH"]
    if bh["n_selected"] != 3 or bh["n_interval_complete"] != 2:
        raise SystemExit("unbounded and missing selected intervals were misclassified")
    if not math.isclose(float(bh["false_coverage_proportion"]), 1.0 / 3.0):
        raise SystemExit("missing selected interval was not retained as noncoverage")
    if bh["n_unbounded_effect_set"] != 1:
        raise SystemExit("unbounded confidence set status was not counted")

    global_cfg = {"binary_options": {}}
    cell = {
        "sample_overlap": {"rf_pqtl": 0.3, "pqtl_outcome": 0.4},
        "analysis_sample_overlap": 0.0,
    }
    options = factorized_options(global_cfg, cell)
    if options["sampling_corr_rf_pqtl"] != 0.0 or options["sampling_corr_pqtl_outcome"] != 0.0:
        raise SystemExit("analysis overlap override did not remain separate from DGP overlap")


if __name__ == "__main__":
    main()
