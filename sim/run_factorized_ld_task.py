#!/usr/bin/env python3
"""Run one cell/scenario/replicate of the factorized genotype-based LD stress test."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
from pathlib import Path

import numpy as np

from run_factorized_ld_stress import (
    LD_CELLS,
    SCENARIOS,
    as_float,
    make_inputs,
    run_one,
    sha256_file,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cell", choices=tuple(LD_CELLS), required=True)
    parser.add_argument("--scenario", choices=SCENARIOS, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--reference-samples", type=int, default=1000)
    parser.add_argument("--analysis-samples", type=int, default=4000)
    parser.add_argument("--blocks", type=int, default=3)
    parser.add_argument("--variants-per-block", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--keep-inputs", action="store_true")
    args = parser.parse_args()

    cell_index = tuple(LD_CELLS).index(args.cell)
    scenario_index = SCENARIOS.index(args.scenario)
    seed = args.seed + cell_index * 1_000_000 + scenario_index * 100_000 + args.replicate
    rng = np.random.default_rng(seed)
    target_rho, reference_rho = LD_CELLS[args.cell]
    rep_dir = args.outdir / args.cell / args.scenario / f"rep_{args.replicate:04d}"
    work = rep_dir / "work"
    ld_prefix, realized_ld = make_inputs(
        work,
        rng,
        args.scenario,
        target_rho,
        reference_rho,
        args.analysis_samples,
        args.reference_samples,
        args.blocks,
        args.variants_per_block,
    )
    result = run_one(args.binary.resolve(), work, ld_prefix)
    row = {
        "cell": args.cell,
        "scenario": args.scenario,
        "replicate": args.replicate,
        "seed": seed,
        "target_rho": target_rho,
        "reference_rho": reference_rho,
        "analysis_samples": args.analysis_samples,
        "reference_samples": args.reference_samples,
        "blocks": args.blocks,
        "variants_per_block": args.variants_per_block,
        "binary_sha256": sha256_file(args.binary.resolve()),
        "realized_causal_pair_ld": realized_ld,
        "factor_nA": as_float(result, "factor_nA"),
        "factor_nB": as_float(result, "factor_nB"),
        "factor_cross_set_max_r2": as_float(result, "factor_cross_set_max_r2"),
        "factor_beta1": as_float(result, "factor_beta1"),
        "factor_beta2": as_float(result, "factor_beta2"),
        "factor_conjunction_p": as_float(result, "factor_conjunction_p"),
        "factor_min_log_BF": as_float(result, "factor_min_log_BF"),
        "factor_log_e_mediation": as_float(result, "factor_log_e_mediation"),
        "factor_p_XM_balanced": as_float(result, "factor_p_XM_balanced"),
        "factor_p_MY_balanced": as_float(result, "factor_p_MY_balanced"),
        "factor_balanced_conjunction_p": as_float(
            result, "factor_balanced_conjunction_p"
        ),
        "factor_log_e_p2e_balanced_mediation": as_float(
            result, "factor_log_e_p2e_balanced_mediation"
        ),
        "factor_log_e_mediation_balanced": as_float(
            result, "factor_log_e_mediation_balanced"
        ),
        "factor_e_q_balanced_EBH": as_float(result, "factor_e_q_balanced_EBH"),
        "factor_e_q_p2e_balanced_EBH": as_float(
            result, "factor_e_q_p2e_balanced_EBH"
        ),
        "factor_log_e_mediation_adaptive": as_float(
            result, "factor_log_e_mediation_adaptive"
        ),
        "factor_log_BF_heterogeneity_MY": as_float(
            result, "factor_log_BF_heterogeneity_MY"
        ),
        "regional_PP_shared": as_float(result, "regional_PP_shared"),
        "regional_PP_distinct": as_float(result, "regional_PP_distinct"),
        "regional_independent_shared_signals": as_float(
            result, "regional_independent_shared_signals"
        ),
        "mediation_identifiability": result.get("mediation_identifiability", ""),
        "factor_two_stage_status": result.get("factor_two_stage_status", ""),
        "factor_mediation_status": result.get("factor_mediation_status", ""),
        "factor_frequentist_status": result.get("factor_frequentist_status", ""),
        "factor_ebh_status": result.get("factor_ebh_status", ""),
        "factor_balanced_status": result.get("factor_balanced_status", ""),
        "factor_adafilter_status": result.get("factor_adafilter_status", ""),
        "factor_balanced_p2e_status": result.get("factor_balanced_p2e_status", ""),
        "factor_balanced_ebh_status": result.get("factor_balanced_ebh_status", ""),
        "factor_adaptive_ebh_status": result.get("factor_adaptive_ebh_status", ""),
    }
    rep_dir.mkdir(parents=True, exist_ok=True)
    with (rep_dir / "task_metrics.tsv").open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    if not args.keep_inputs:
        shutil.rmtree(work)
    print(f"Wrote {rep_dir / 'task_metrics.tsv'}")


if __name__ == "__main__":
    main()
