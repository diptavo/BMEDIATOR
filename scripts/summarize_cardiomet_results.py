from __future__ import annotations

import csv
import math
import os
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    "/Users/duttad5/RFMediation/.codex_tasks/cardiomet_summary/.mplconfig",
)

import matplotlib.pyplot as plt


REMOTE_BASE = (
    "/data/Dutta_lab/BMEDIATOR/analysis/"
    "cardiomet_vascular_kidney_ukbb_combined_tss_pm6_r2_0.1_hg38/results"
)
LOCAL_BASE = Path("/Users/duttad5/RFMediation/.codex_tasks/cardiomet_summary")
SUMMARY_DIR = LOCAL_BASE / "summary"
PLOT_DIR = LOCAL_BASE / "plots"


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


def ensure_dirs() -> None:
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_remote_results() -> list[dict[str, str]]:
    remote_py = r"""
from pathlib import Path
import csv
import json

base = Path(""" + repr(REMOTE_BASE) + r""")
rows = []
for med in sorted(base.rglob("*.mediation")):
    rel = med.relative_to(base)
    parts = rel.parts
    domain = parts[0]
    outcome_dir = parts[1]
    exposure_dir = parts[2]
    stem = med.name.replace("\r", "")
    out_prefix = stem[:-10] if stem.endswith(".mediation") else stem
    outcome_label_map = {
        "chd_finngen": "CHD_FinnGen",
        "mi_finngen": "MI_FinnGen",
        "stroke_finngen": "Stroke_FinnGen",
        "heartfail_finngen": "HeartFail_FinnGen",
        "padproxy_finngen": "PADProxy_FinnGen",
        "ckd_ckdgen": "CKD_CKDGen",
        "microalbuminuria_ckdgen": "Microalbuminuria_CKDGen",
        "uacr_ckdgen": "UACR_CKDGen",
        "renfail_finngen": "RenFail_FinnGen",
        "kidneystones_finngen": "KidneyStones_FinnGen",
        "egfrcrea_ckdgen": "eGFRCrea_CKDGen",
        "egfrcys_ckdgen": "eGFRCys_CKDGen",
    }
    exposure_label_map = {
        "bmi_giant": "BMI_GIANT",
        "bmi_mvp": "BMI_MVP",
        "whr_giant": "WHR_GIANT",
        "whradjbmi_giant": "WHRadjBMI_GIANT",
        "sbp_keaton": "SBP_KEATON",
        "sbp_mvp": "SBP_MVP",
        "dbp_keaton": "DBP_KEATON",
        "dbp_mvp": "DBP_MVP",
        "pp_keaton": "PP_KEATON",
        "hdl_glgc": "HDL_GLGC",
        "ldl_glgc": "LDL_GLGC",
        "tg_glgc": "TG_GLGC",
        "tc_glgc": "TC_GLGC",
    }
    outcome = outcome_label_map.get(outcome_dir, outcome_dir)
    exposure = exposure_label_map.get(exposure_dir, exposure_dir)
    hyp = med.with_name(med.name.replace(".mediation", ".hyp"))
    if not hyp.exists():
        hyp = Path(str(med).replace(".mediation", ".hyp"))

    with med.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        hit_rows = list(reader)

    def yes_count(key):
        return sum(1 for r in hit_rows if r.get(key, "").upper() == "YES")

    def safe_float(v):
        try:
            x = float(v)
            if x != x:
                return None
            return x
        except Exception:
            return None

    pm1_vals = [safe_float(r.get("P_M1")) for r in hit_rows]
    pm1_vals = [x for x in pm1_vals if x is not None]
    best = None
    if hit_rows:
        best = max(
            hit_rows,
            key=lambda r: (
                1 if r.get("selected_fdr5", "").upper() == "YES" else 0,
                safe_float(r.get("P_M1")) or -1.0,
            ),
        )

    hyp_map = {}
    if hyp.exists():
        for line in hyp.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            toks = line.split("\t")
            if len(toks) == 2:
                hyp_map[toks[0]] = toks[1]

    rows.append(
        {
            "domain": domain,
            "exposure": exposure,
            "outcome": outcome,
            "out_prefix": out_prefix,
            "protein_rows": str(len(hit_rows)),
            "fdr10_hits": str(yes_count("selected_fdr10")),
            "fdr5_hits": str(yes_count("selected_fdr5")),
            "strong_hits": str(sum(1 for r in hit_rows if r.get("evidence_tier") == "strong")),
            "suggestive_hits": str(sum(1 for r in hit_rows if r.get("evidence_tier") == "suggestive")),
            "top_gene": best.get("Gene", "") if best else "",
            "top_protein": best.get("Protein", "") if best else "",
            "top_pm1": "" if not best else str(safe_float(best.get("P_M1")) or ""),
            "top_selected_fdr5": "" if not best else best.get("selected_fdr5", ""),
            "max_pm1": "" if not pm1_vals else str(max(pm1_vals)),
            "p1": hyp_map.get("p1", ""),
            "p2": hyp_map.get("p2", ""),
            "p3": hyp_map.get("p3", ""),
        }
    )

print(json.dumps(rows))
"""
    result = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "biowulf", "python3", "-c", remote_py],
        check=True,
        capture_output=True,
        text=True,
    )
    import json

    return json.loads(result.stdout)


def fetch_remote_fdr5_hits() -> list[dict[str, str]]:
    remote_py = r"""
from pathlib import Path
import csv
import json

base = Path(""" + repr(REMOTE_BASE) + r""")
rows = []
outcome_label_map = {
    "chd_finngen": "CHD_FinnGen",
    "mi_finngen": "MI_FinnGen",
    "stroke_finngen": "Stroke_FinnGen",
    "heartfail_finngen": "HeartFail_FinnGen",
    "padproxy_finngen": "PADProxy_FinnGen",
    "ckd_ckdgen": "CKD_CKDGen",
    "microalbuminuria_ckdgen": "Microalbuminuria_CKDGen",
    "uacr_ckdgen": "UACR_CKDGen",
    "renfail_finngen": "RenFail_FinnGen",
    "kidneystones_finngen": "KidneyStones_FinnGen",
    "egfrcrea_ckdgen": "eGFRCrea_CKDGen",
    "egfrcys_ckdgen": "eGFRCys_CKDGen",
}
exposure_label_map = {
    "bmi_giant": "BMI_GIANT",
    "bmi_mvp": "BMI_MVP",
    "whr_giant": "WHR_GIANT",
    "whradjbmi_giant": "WHRadjBMI_GIANT",
    "sbp_keaton": "SBP_KEATON",
    "sbp_mvp": "SBP_MVP",
    "dbp_keaton": "DBP_KEATON",
    "dbp_mvp": "DBP_MVP",
    "pp_keaton": "PP_KEATON",
    "hdl_glgc": "HDL_GLGC",
    "ldl_glgc": "LDL_GLGC",
    "tg_glgc": "TG_GLGC",
    "tc_glgc": "TC_GLGC",
}
for med in sorted(base.rglob("*.mediation")):
    rel = med.relative_to(base)
    domain = rel.parts[0]
    outcome = outcome_label_map.get(rel.parts[1], rel.parts[1])
    exposure = exposure_label_map.get(rel.parts[2], rel.parts[2])
    with med.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r.get("selected_fdr5", "").upper() != "YES":
                continue
            rows.append(
                {
                    "domain": domain,
                    "exposure": exposure,
                    "outcome": outcome,
                    "protein": r.get("Protein", ""),
                    "gene": r.get("Gene", ""),
                    "P_M1": r.get("P_M1", ""),
                    "posterior_local_fdr": r.get("posterior_local_fdr", ""),
                    "evidence_tier": r.get("evidence_tier", ""),
                    "mediated_effect": r.get("mediated_effect", ""),
                }
            )
print(json.dumps(rows))
"""
    result = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "biowulf", "python3", "-c", remote_py],
        check=True,
        capture_output=True,
        text=True,
    )
    import json

    return json.loads(result.stdout)


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize_run_tables(run_rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    by_outcome = []
    by_exposure = []

    outcome_map: dict[tuple[str, str], dict[str, int | float]] = defaultdict(
        lambda: {"runs": 0, "fdr5_hits": 0, "fdr10_hits": 0, "runs_with_fdr5": 0, "max_pm1": 0.0}
    )
    exposure_map: dict[tuple[str, str], dict[str, int | float]] = defaultdict(
        lambda: {"runs": 0, "fdr5_hits": 0, "fdr10_hits": 0, "runs_with_fdr5": 0, "max_pm1": 0.0}
    )

    for row in run_rows:
        domain = row["domain"]
        outcome = row["outcome"]
        exposure = row["exposure"]
        fdr5 = int(row["fdr5_hits"])
        fdr10 = int(row["fdr10_hits"])
        max_pm1 = float(row["max_pm1"]) if row["max_pm1"] else 0.0
        o = outcome_map[(domain, outcome)]
        e = exposure_map[(domain, exposure)]
        for agg in (o, e):
            agg["runs"] += 1
            agg["fdr5_hits"] += fdr5
            agg["fdr10_hits"] += fdr10
            agg["runs_with_fdr5"] += 1 if fdr5 > 0 else 0
            agg["max_pm1"] = max(agg["max_pm1"], max_pm1)

    for (domain, outcome), vals in sorted(outcome_map.items()):
        by_outcome.append(
            {
                "domain": domain,
                "outcome": outcome,
                "runs": str(vals["runs"]),
                "fdr5_hits_total": str(vals["fdr5_hits"]),
                "fdr10_hits_total": str(vals["fdr10_hits"]),
                "runs_with_fdr5": str(vals["runs_with_fdr5"]),
                "max_pm1": f"{vals['max_pm1']:.6f}",
            }
        )
    for (domain, exposure), vals in sorted(exposure_map.items()):
        by_exposure.append(
            {
                "domain": domain,
                "exposure": exposure,
                "runs": str(vals["runs"]),
                "fdr5_hits_total": str(vals["fdr5_hits"]),
                "fdr10_hits_total": str(vals["fdr10_hits"]),
                "runs_with_fdr5": str(vals["runs_with_fdr5"]),
                "max_pm1": f"{vals['max_pm1']:.6f}",
            }
        )
    return by_outcome, by_exposure


def matrix_from_runs(run_rows: list[dict[str, str]], outcomes: list[str]) -> dict[str, list[int]]:
    matrix: dict[str, list[int]] = {}
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
    for i, exposure in enumerate(exposures):
        for j, _ in enumerate(outcomes):
            v = data[i][j]
            ax.text(j, i, str(v), ha="center", va="center", fontsize=7, color="#1f1f1f")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("FDR<5% mediators")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_top_genes(hit_rows: list[dict[str, str]], out_path: Path) -> list[tuple[str, int]]:
    counter = Counter(r["gene"] for r in hit_rows if r["gene"])
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
    return top


def main() -> None:
    ensure_dirs()
    run_rows = fetch_remote_results()
    hit_rows = fetch_remote_fdr5_hits()

    run_rows = sorted(run_rows, key=lambda r: (r["domain"], r["outcome"], r["exposure"]))
    hit_rows = sorted(hit_rows, key=lambda r: (r["domain"], r["outcome"], r["exposure"], r["gene"]))

    write_tsv(
        SUMMARY_DIR / "run_summary.tsv",
        run_rows,
        [
            "domain",
            "exposure",
            "outcome",
            "out_prefix",
            "protein_rows",
            "fdr10_hits",
            "fdr5_hits",
            "strong_hits",
            "suggestive_hits",
            "top_gene",
            "top_protein",
            "top_pm1",
            "top_selected_fdr5",
            "max_pm1",
            "p1",
            "p2",
            "p3",
        ],
    )
    write_tsv(
        SUMMARY_DIR / "fdr5_hits.tsv",
        hit_rows,
        ["domain", "exposure", "outcome", "protein", "gene", "P_M1", "posterior_local_fdr", "evidence_tier", "mediated_effect"],
    )

    by_outcome, by_exposure = summarize_run_tables(run_rows)
    write_tsv(
        SUMMARY_DIR / "outcome_summary.tsv",
        by_outcome,
        ["domain", "outcome", "runs", "fdr5_hits_total", "fdr10_hits_total", "runs_with_fdr5", "max_pm1"],
    )
    write_tsv(
        SUMMARY_DIR / "exposure_summary.tsv",
        by_exposure,
        ["domain", "exposure", "runs", "fdr5_hits_total", "fdr10_hits_total", "runs_with_fdr5", "max_pm1"],
    )

    vascular_runs = [r for r in run_rows if r["domain"] == "vascular"]
    kidney_runs = [r for r in run_rows if r["domain"] == "kidney"]

    plot_heatmap(
        matrix_from_runs(vascular_runs, VASCULAR_OUTCOMES),
        VASCULAR_OUTCOMES,
        "Vascular outcomes: FDR<5% mediator counts",
        PLOT_DIR / "vascular_fdr5_heatmap.png",
    )
    plot_heatmap(
        matrix_from_runs(kidney_runs, KIDNEY_OUTCOMES),
        KIDNEY_OUTCOMES,
        "Kidney outcomes: FDR<5% mediator counts",
        PLOT_DIR / "kidney_fdr5_heatmap.png",
    )
    top_genes = plot_top_genes(hit_rows, PLOT_DIR / "top_recurrent_genes_fdr5.png")

    with (SUMMARY_DIR / "headline_summary.txt").open("w") as handle:
        handle.write(f"run_count\t{len(run_rows)}\n")
        handle.write(f"fdr5_hit_rows\t{len(hit_rows)}\n")
        handle.write(f"vascular_runs\t{len(vascular_runs)}\n")
        handle.write(f"kidney_runs\t{len(kidney_runs)}\n")
        if top_genes:
            handle.write(f"top_gene\t{top_genes[0][0]}\n")
            handle.write(f"top_gene_run_count\t{top_genes[0][1]}\n")


if __name__ == "__main__":
    main()
