from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

from benchmark_lib import (
    attach_truth,
    benchmark_scenarios,
    compute_benchmark_metrics,
    read_results,
    read_truth,
    scenario_alias_map,
    summarize_calibration,
    summarize_classification,
    summarize_fdr_power,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize BMEDIATOR simulation benchmarks.")
    parser.add_argument("--outdir", type=Path, required=True, help="Benchmark output directory.")
    parser.add_argument("--rebuild", action="store_true", help="Rebuild summary tables from replicate outputs before plotting.")
    return parser.parse_args()


def rebuild_summary_tables(outdir: Path) -> None:
    summary_root = outdir / "summary"
    summary_root.mkdir(parents=True, exist_ok=True)
    config = None
    manifest_path = outdir / "manifests" / "benchmark_manifest.tsv"
    if manifest_path.exists():
        pass

    for benchmark in ("classification", "calibration"):
        benchmark_root = outdir / benchmark
        if not benchmark_root.exists():
            continue
        alias_map = {}
        summary_scenarios = None
        config_path = outdir / "config.json"
        if config_path.exists():
            from benchmark_lib import load_config

            config = load_config(config_path)
            alias_map = scenario_alias_map(config, benchmark)
            scenarios = benchmark_scenarios(config, benchmark)
            summary_scenarios = tuple(dict.fromkeys(alias_map.get(s, s) for s in scenarios))
        task_frames = []
        class_frames = []
        confusion_frames = []
        calibration_frames = []
        bfdr_frames = []
        fdr_power_frames = []

        for rep_dir in benchmark_root.glob("*/*"):
            if not rep_dir.is_dir():
                continue
            task_metrics_path = rep_dir / "task_metrics.tsv"
            if task_metrics_path.exists():
                merged = read_results(task_metrics_path)
            else:
                truth_path = rep_dir / "truth.tsv"
                mediation_path = rep_dir / "bmediator.mediation"
                if not truth_path.exists() or not mediation_path.exists():
                    continue
                truth = read_truth(truth_path)
                results = read_results(mediation_path)
                merged = compute_benchmark_metrics(attach_truth(results, truth), alias_map=alias_map)
                for row in merged:
                    row["benchmark"] = benchmark
                    row["cell"] = rep_dir.parent.name
                    row["replicate"] = int(rep_dir.name.split("_")[-1])
            task_frames.append(merged)
            fdr_power = summarize_fdr_power(merged)
            for row in fdr_power:
                row["benchmark"] = benchmark
                row["cell"] = rep_dir.parent.name
                row["replicate"] = int(rep_dir.name.split("_")[-1])
            fdr_power_frames.append(fdr_power)

            if benchmark == "classification":
                by_scenario, confusion = summarize_classification(merged, scenarios=summary_scenarios)
                for row in by_scenario:
                    row["cell"] = rep_dir.parent.name
                    row["replicate"] = merged[0]["replicate"]
                for row in confusion:
                    row["cell"] = rep_dir.parent.name
                    row["replicate"] = merged[0]["replicate"]
                class_frames.append(by_scenario)
                confusion_frames.append(confusion)
            else:
                calibration, bfdr = summarize_calibration(merged)
                for row in calibration:
                    row["cell"] = rep_dir.parent.name
                    row["replicate"] = merged[0]["replicate"]
                for row in bfdr:
                    row["cell"] = rep_dir.parent.name
                    row["replicate"] = merged[0]["replicate"]
                calibration_frames.append(calibration)
                bfdr_frames.append(bfdr)

        bench_summary = summary_root / benchmark
        bench_summary.mkdir(parents=True, exist_ok=True)
        if task_frames:
            rows = [row for chunk in task_frames for row in chunk]
            with (bench_summary / "protein_level_metrics.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow(row)
        if class_frames:
            rows = [row for chunk in class_frames for row in chunk]
            with (bench_summary / "classification_by_scenario.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow(row)
        if confusion_frames:
            rows = [row for chunk in confusion_frames for row in chunk]
            with (bench_summary / "classification_confusion.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow(row)
        if calibration_frames:
            rows = [row for chunk in calibration_frames for row in chunk]
            with (bench_summary / "calibration_bins.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow(row)
        if bfdr_frames:
            rows = [row for chunk in bfdr_frames for row in chunk]
            with (bench_summary / "calibration_bfdr.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow(row)
        if fdr_power_frames:
            rows = [row for chunk in fdr_power_frames for row in chunk]
            with (bench_summary / "fdr_power_by_threshold.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow(row)


def plot_classification(summary_dir: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    by_scenario = summary_dir / "classification" / "classification_by_scenario.tsv"
    confusion = summary_dir / "classification" / "classification_confusion.tsv"
    if not by_scenario.exists() or not confusion.exists():
        return

    by_rows = list(csv.DictReader(by_scenario.open(), delimiter="\t"))
    agg_map = {}
    for row in by_rows:
        key = (row["cell"], row["true_scenario"])
        agg_map.setdefault(key, []).append(float(row["accuracy"]))
    agg = [{"cell": cell, "true_scenario": scenario, "accuracy": sum(vals) / len(vals)} for (cell, scenario), vals in agg_map.items()]
    cells = list(dict.fromkeys(row["cell"] for row in agg))
    conf_rows = list(csv.DictReader(confusion.open(), delimiter="\t"))
    scenarios = list(dict.fromkeys(
        [row["true_scenario"] for row in conf_rows]
        + [row["pred_scenario"] for row in conf_rows]
        + [row["true_scenario"] for row in by_rows]
    ))

    fig, ax = plt.subplots(figsize=(10, 4.8))
    width = 0.22
    xs = range(len(scenarios))
    for idx, cell in enumerate(cells):
        cur = {row["true_scenario"]: row["accuracy"] for row in agg if row["cell"] == cell}
        ax.bar([x + (idx - (len(cells) - 1) / 2) * width for x in xs], [cur.get(s, 0.0) for s in scenarios], width=width, label=cell)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(scenarios)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Classification accuracy")
    ax.set_title("BMEDIATOR classification benchmark")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(summary_dir / "classification" / "classification_accuracy.png", dpi=200)
    plt.close(fig)

    pivot = np.zeros((len(scenarios), len(scenarios)))
    for row in conf_rows:
        i = scenarios.index(row["true_scenario"])
        j = scenarios.index(row["pred_scenario"])
        pivot[i, j] += float(row["count"])

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(pivot, cmap="Blues")
    ax.set_xticks(range(len(scenarios)))
    ax.set_yticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios)
    ax.set_yticklabels(scenarios)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title("Classification confusion")
    for i in range(len(scenarios)):
        for j in range(len(scenarios)):
            ax.text(j, i, int(pivot[i, j]), ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(summary_dir / "classification" / "classification_confusion.png", dpi=200)
    plt.close(fig)


def plot_calibration(summary_dir: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    bins_path = summary_dir / "calibration" / "calibration_bins.tsv"
    bfdr_path = summary_dir / "calibration" / "calibration_bfdr.tsv"
    if not bins_path.exists() or not bfdr_path.exists():
        return

    bin_rows = list(csv.DictReader(bins_path.open(), delimiter="\t"))
    bin_map = {}
    for row in bin_rows:
        key = (row["cell"], row["bin_mid"])
        bin_map.setdefault(key, {"mean_pred": [], "observed_m1": [], "n": 0})
        bin_map[key]["mean_pred"].append(float(row["mean_pred"]))
        bin_map[key]["observed_m1"].append(float(row["observed_m1"]))
        bin_map[key]["n"] += int(float(row["n"]))
    bins_agg = [
        {
            "cell": cell,
            "bin_mid": float(bin_mid),
            "mean_pred": sum(vals["mean_pred"]) / len(vals["mean_pred"]),
            "observed_m1": sum(vals["observed_m1"]) / len(vals["observed_m1"]),
            "n": vals["n"],
        }
        for (cell, bin_mid), vals in bin_map.items()
    ]
    cells = list(dict.fromkeys(row["cell"] for row in bins_agg))

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    ax.plot([0, 1], [0, 1], linestyle="--", color="#666666", linewidth=1)
    for cell in cells:
        cur = sorted([row for row in bins_agg if row["cell"] == cell], key=lambda row: row["bin_mid"])
        ax.plot([row["mean_pred"] for row in cur], [row["observed_m1"] for row in cur], marker="o", linewidth=1.8, label=cell)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Mean predicted P(M1)")
    ax.set_ylabel("Observed fraction of true M1")
    ax.set_title("Calibration of posterior P(M1)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(summary_dir / "calibration" / "calibration_curve.png", dpi=200)
    plt.close(fig)

    bfdr_rows = list(csv.DictReader(bfdr_path.open(), delimiter="\t"))
    bfdr_map = {}
    for row in bfdr_rows:
        key = (row["cell"], row["rank_bin"])
        bfdr_map.setdefault(key, {"mean_estimated_bfdr": [], "empirical_fdr": []})
        bfdr_map[key]["mean_estimated_bfdr"].append(float(row["mean_estimated_bfdr"]))
        bfdr_map[key]["empirical_fdr"].append(float(row["empirical_fdr"]))
    bfdr_agg = [
        {
            "cell": cell,
            "rank_bin": rank_bin,
            "mean_estimated_bfdr": sum(vals["mean_estimated_bfdr"]) / len(vals["mean_estimated_bfdr"]),
            "empirical_fdr": sum(vals["empirical_fdr"]) / len(vals["empirical_fdr"]),
        }
        for (cell, rank_bin), vals in bfdr_map.items()
    ]
    rank_bins = list(dict.fromkeys(row["rank_bin"] for row in bfdr_agg))

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    x_positions = range(len(rank_bins))
    width = 0.35
    for idx, cell in enumerate(cells):
        cur_map = {row["rank_bin"]: row for row in bfdr_agg if row["cell"] == cell}
        cur = [cur_map[rank_bin] for rank_bin in rank_bins if rank_bin in cur_map]
        xs = [x + (idx - (len(cells) - 1) / 2) * width for x in x_positions]
        ax.bar(xs, [row["empirical_fdr"] for row in cur], width=width, alpha=0.8, label=f"{cell} empirical")
        ax.plot(xs, [row["mean_estimated_bfdr"] for row in cur], color="black", marker="o", linewidth=1)
    ax.set_xticks(list(x_positions))
    ax.set_xticklabels(rank_bins)
    ax.set_ylabel("FDR")
    ax.set_title("Estimated vs empirical Bayesian FDR by rank bin")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(summary_dir / "calibration" / "bfdr_comparison.png", dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    os.environ.setdefault("MPLCONFIGDIR", str(args.outdir / ".mplconfig"))
    if args.rebuild:
        rebuild_summary_tables(args.outdir)
    summary_dir = args.outdir / "summary"
    plot_classification(summary_dir)
    plot_calibration(summary_dir)


if __name__ == "__main__":
    main()
