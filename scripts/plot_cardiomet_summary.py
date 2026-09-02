from __future__ import annotations

import argparse
import csv
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    "/Users/duttad5/RFMediation/.codex_tasks/cardiomet_summary/.mplconfig",
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXPOSURE_ORDER = [
    "BMI_GIANT",
    "BMI_MVP",
    "WHR_GIANT",
    "WHRadjBMI_GIANT",
    "SBP_KEATON",
    "SBP_MVP",
    "DBP_KEATON",
    "DBP_MVP",
    "PP_KEATON",
    "HDL_GLGC",
    "LDL_GLGC",
    "TG_GLGC",
    "TC_GLGC",
]

VASCULAR_OUTCOMES = [
    "CHD_FinnGen",
    "MI_FinnGen",
    "Stroke_FinnGen",
    "HeartFail_FinnGen",
    "PADProxy_FinnGen",
]

KIDNEY_OUTCOMES = [
    "CKD_CKDGen",
    "Microalbuminuria_CKDGen",
    "UACR_CKDGen",
    "RenFail_FinnGen",
    "KidneyStones_FinnGen",
    "eGFRCrea_CKDGen",
    "eGFRCys_CKDGen",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def matrix_from_runs(run_rows: list[dict[str, str]], outcomes: list[str]) -> dict[str, list[int]]:
    matrix = {}
    for exposure in EXPOSURE_ORDER:
        vals = []
        for outcome in outcomes:
            row = next((r for r in run_rows if r["exposure"] == exposure and r["outcome"] == outcome), None)
            vals.append(int(row["fdr5_hits"]) if row else 0)
        matrix[exposure] = vals
    return matrix


def plot_heatmap(matrix: dict[str, list[int]], outcomes: list[str], title: str, out_path: Path) -> None:
    exposures = [e for e in EXPOSURE_ORDER if e in matrix]
    data = [matrix[e] for e in exposures]
    vmax = max(max(row) for row in data) if data else 1
    fig, ax = plt.subplots(figsize=(1.35 * len(outcomes) + 3.5, 0.42 * len(exposures) + 2.2), dpi=180)
    im = ax.imshow(data, cmap="YlOrRd", aspect="auto", vmin=0, vmax=max(vmax, 1))
    ax.set_xticks(range(len(outcomes)))
    ax.set_xticklabels(outcomes, rotation=35, ha="right")
    ax.set_yticks(range(len(exposures)))
    ax.set_yticklabels(exposures)
    ax.set_title(title, fontsize=13, weight="bold")
    for i in range(len(exposures)):
        for j in range(len(outcomes)):
            v = data[i][j]
            ax.text(j, i, str(v), ha="center", va="center", fontsize=7, color="#1f1f1f")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("FDR<5% mediators")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_top_genes(hit_rows: list[dict[str, str]], out_path: Path) -> None:
    counter = Counter(r["gene"] for r in hit_rows if r.get("gene"))
    top = counter.most_common(20)
    genes = [g for g, _ in top][::-1]
    counts = [c for _, c in top][::-1]
    fig, ax = plt.subplots(figsize=(9.5, 6.5), dpi=180)
    ax.barh(genes, counts, color="#2c7fb8")
    for y, c in enumerate(counts):
        ax.text(c + 0.1, y, str(c), va="center", fontsize=8)
    ax.set_xlabel("Number of exposure-outcome runs with FDR<5% mediation")
    ax.set_title("Most recurrent FDR<5% mediator genes", fontsize=13, weight="bold")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-dir", required=True)
    ap.add_argument("--plot-dir", required=True)
    args = ap.parse_args()

    summary_dir = Path(args.summary_dir)
    plot_dir = Path(args.plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    run_rows = read_tsv(summary_dir / "run_summary.tsv")
    hit_rows = read_tsv(summary_dir / "fdr5_hits.tsv")

    vascular_runs = [r for r in run_rows if r["domain"] == "vascular"]
    kidney_runs = [r for r in run_rows if r["domain"] == "kidney"]

    plot_heatmap(
        matrix_from_runs(vascular_runs, VASCULAR_OUTCOMES),
        VASCULAR_OUTCOMES,
        "Vascular outcomes: FDR<5% mediator counts",
        plot_dir / "vascular_fdr5_heatmap.png",
    )
    plot_heatmap(
        matrix_from_runs(kidney_runs, KIDNEY_OUTCOMES),
        KIDNEY_OUTCOMES,
        "Kidney outcomes: FDR<5% mediator counts",
        plot_dir / "kidney_fdr5_heatmap.png",
    )
    plot_top_genes(hit_rows, plot_dir / "top_recurrent_genes_fdr5.png")


if __name__ == "__main__":
    main()
