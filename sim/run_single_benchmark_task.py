from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from benchmark_lib import (
    attach_truth,
    benchmark_scenarios,
    calibration_sequence,
    classification_sequence,
    compute_benchmark_metrics,
    generate_proteins,
    load_config,
    read_results,
    read_truth,
    run_bmediator,
    scenario_alias_map,
    write_dataset,
)
from run_competitor_benchmark import load_sumstats, method_score_rows, protein_rows, to_float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one BMEDIATOR simulation task.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--benchmark", choices=["classification", "calibration"], required=True)
    parser.add_argument("--cell", required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--binary", type=Path, default=Path(__file__).resolve().parents[1] / "bmediator")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--keep-instruments", action="store_true", help="Keep bulky .instruments files after task metrics are written.")
    parser.add_argument("--keep-intermediate", action="store_true", help="Keep generated inputs and raw BMEDIATOR outputs after compact metrics are written.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "config.json").write_text(args.config.read_text())
    global_cfg = config["global"]
    bench_cfg = config[args.benchmark]
    cell_cfg = next(cell for cell in bench_cfg["cells"] if cell["name"] == args.cell)
    if any(
        float(cell_cfg.get(key, 0.0)) != 0.0
        for key in ("linked_pqtl_outcome_leakage", "linked_rf_outcome_leakage")
    ):
        raise SystemExit(
            "This cell uses the deprecated pseudo-LD outcome-injection model. "
            "Use sim/run_regional_ld_stress.py, which generates correlated genotypes, "
            "distinct causal variants, regional summary statistics, and a PLINK panel."
        )
    scenarios = benchmark_scenarios(config, args.benchmark)
    alias_map = scenario_alias_map(config, args.benchmark)

    base_seed = int(global_cfg["seed"] if args.seed is None else args.seed)
    cell_index = next(i for i, cell in enumerate(bench_cfg["cells"]) if cell["name"] == args.cell)
    rng = np.random.default_rng(base_seed + cell_index * 10000 + args.replicate)

    if args.benchmark == "classification":
        scenario_sequence = classification_sequence(int(bench_cfg["proteins_per_scenario"]), scenarios)
    else:
        scenario_sequence = calibration_sequence(
            rng,
            int(bench_cfg["proteins_per_replicate"]),
            bench_cfg["scenario_mix"],
            scenarios,
        )

    proteins = generate_proteins(rng, args.benchmark, cell_cfg, scenario_sequence)
    files = write_dataset(args.outdir, args.benchmark, args.cell, args.replicate, proteins)
    out_prefix = files.truth.parent / "bmediator"
    mediation_path, _ = run_bmediator(args.binary, files, out_prefix, global_cfg)

    truth = read_truth(files.truth)
    results = read_results(mediation_path)
    merged = compute_benchmark_metrics(attach_truth(results, truth), alias_map=alias_map)
    for row in merged:
        row["benchmark"] = args.benchmark
        row["cell"] = args.cell
        row["replicate"] = args.replicate
    task_metrics = files.truth.parent / "task_metrics.tsv"
    with task_metrics.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(merged[0].keys()), delimiter="\t")
        writer.writeheader()
        for row in merged:
            writer.writerow(row)

    rf, outcome, pqtl_by_protein, info = load_sumstats(files.truth.parent)
    bmed_p = {}
    for row in results:
        protein = str(row.get("Protein", row.get("protein_id", "")))
        if protein:
            bmed_p[protein] = to_float(row.get("P_M1", row.get("prob_M1")))
    competitor_rows = []
    for truth_row in truth:
        protein = truth_row["protein_id"]
        all_rows, rf_inst, cis_inst = protein_rows(
            protein,
            rf,
            outcome,
            pqtl_by_protein,
            info,
            float(global_cfg["rf_p_threshold"]),
            float(global_cfg["cis_p_threshold"]),
            int(global_cfg.get("cis_window_bp", 1_000_000)),
        )
        rows = method_score_rows(protein, truth_row["true_scenario"], all_rows, rf_inst, cis_inst, bmed_p.get(protein))
        for row in rows:
            row["benchmark"] = args.benchmark
            row["cell"] = args.cell
            row["replicate"] = args.replicate
            competitor_rows.append(row)
    competitor_path = files.truth.parent / "competitor_scores.tsv"
    if competitor_rows:
        with competitor_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(competitor_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            for row in competitor_rows:
                writer.writerow(row)

    instruments = Path(f"{out_prefix}.instruments")
    if instruments.exists() and not args.keep_instruments:
        instruments.unlink()
    if not args.keep_intermediate:
        for path in (
            files.rf,
            files.pqtl,
            files.cancer,
            files.protein_info,
            files.truth,
            Path(f"{out_prefix}.mediation"),
            Path(f"{out_prefix}.hyp"),
        ):
            if path.exists():
                path.unlink()
    print(f"Wrote {task_metrics}")


if __name__ == "__main__":
    main()
