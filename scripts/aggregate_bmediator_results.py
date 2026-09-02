from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


OUTCOME_LABEL_MAP = {
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

EXPOSURE_LABEL_MAP = {
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


def safe_float(v: str) -> float | None:
    try:
        x = float(v)
        if x != x:
            return None
        return x
    except Exception:
        return None


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def aggregate(results_dir: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    run_rows: list[dict[str, str]] = []
    hit_rows: list[dict[str, str]] = []

    for med in sorted(results_dir.rglob("*.mediation")):
        rel = med.relative_to(results_dir)
        domain = rel.parts[0]
        outcome = OUTCOME_LABEL_MAP.get(rel.parts[1], rel.parts[1])
        exposure = EXPOSURE_LABEL_MAP.get(rel.parts[2], rel.parts[2])
        stem = med.name.replace("\r", "")
        out_prefix = stem[:-10] if stem.endswith(".mediation") else stem

        hyp = med.with_name(med.name.replace(".mediation", ".hyp"))
        if not hyp.exists():
            hyp = Path(str(med).replace(".mediation", ".hyp"))

        with med.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        pm1_vals = [safe_float(r.get("P_M1", "")) for r in rows]
        pm1_vals = [x for x in pm1_vals if x is not None]

        best = None
        if rows:
            best = max(
                rows,
                key=lambda r: (
                    1 if r.get("selected_fdr5", "").upper() == "YES" else 0,
                    safe_float(r.get("P_M1", "")) or -1.0,
                ),
            )

        hyp_map: dict[str, str] = {}
        if hyp.exists():
            for line in hyp.read_text().splitlines():
                if not line or line.startswith("#"):
                    continue
                toks = line.split("\t")
                if len(toks) == 2:
                    hyp_map[toks[0]] = toks[1]

        run_rows.append(
            {
                "domain": domain,
                "exposure": exposure,
                "outcome": outcome,
                "out_prefix": out_prefix,
                "protein_rows": str(len(rows)),
                "fdr10_hits": str(sum(1 for r in rows if r.get("selected_fdr10", "").upper() == "YES")),
                "fdr5_hits": str(sum(1 for r in rows if r.get("selected_fdr5", "").upper() == "YES")),
                "strong_hits": str(sum(1 for r in rows if r.get("evidence_tier") == "strong")),
                "suggestive_hits": str(sum(1 for r in rows if r.get("evidence_tier") == "suggestive")),
                "top_gene": best.get("Gene", "") if best else "",
                "top_protein": best.get("Protein", "") if best else "",
                "top_pm1": "" if not best else str(safe_float(best.get("P_M1", "")) or ""),
                "top_selected_fdr5": "" if not best else best.get("selected_fdr5", ""),
                "max_pm1": "" if not pm1_vals else str(max(pm1_vals)),
                "p1": hyp_map.get("p1", ""),
                "p2": hyp_map.get("p2", ""),
                "p3": hyp_map.get("p3", ""),
            }
        )

        for r in rows:
            if r.get("selected_fdr5", "").upper() != "YES":
                continue
            hit_rows.append(
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

    return run_rows, hit_rows


def summarize(run_rows: list[dict[str, str]], hit_rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    by_outcome: list[dict[str, str]] = []
    by_exposure: list[dict[str, str]] = []
    gene_rows: list[dict[str, str]] = []

    outcome_map = defaultdict(lambda: {"runs": 0, "fdr5_hits": 0, "fdr10_hits": 0, "runs_with_fdr5": 0, "max_pm1": 0.0})
    exposure_map = defaultdict(lambda: {"runs": 0, "fdr5_hits": 0, "fdr10_hits": 0, "runs_with_fdr5": 0, "max_pm1": 0.0})
    gene_counter = Counter(r["gene"] for r in hit_rows if r["gene"])

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

    for gene, count in gene_counter.most_common():
        gene_rows.append({"gene": gene, "fdr5_run_count": str(count)})

    return by_outcome, by_exposure, gene_rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_rows, hit_rows = aggregate(results_dir)
    by_outcome, by_exposure, gene_rows = summarize(run_rows, hit_rows)

    run_rows = sorted(run_rows, key=lambda r: (r["domain"], r["outcome"], r["exposure"]))
    hit_rows = sorted(hit_rows, key=lambda r: (r["domain"], r["outcome"], r["exposure"], r["gene"]))

    write_tsv(
        output_dir / "run_summary.tsv",
        run_rows,
        ["domain", "exposure", "outcome", "out_prefix", "protein_rows", "fdr10_hits", "fdr5_hits", "strong_hits", "suggestive_hits", "top_gene", "top_protein", "top_pm1", "top_selected_fdr5", "max_pm1", "p1", "p2", "p3"],
    )
    write_tsv(
        output_dir / "fdr5_hits.tsv",
        hit_rows,
        ["domain", "exposure", "outcome", "protein", "gene", "P_M1", "posterior_local_fdr", "evidence_tier", "mediated_effect"],
    )
    write_tsv(
        output_dir / "outcome_summary.tsv",
        by_outcome,
        ["domain", "outcome", "runs", "fdr5_hits_total", "fdr10_hits_total", "runs_with_fdr5", "max_pm1"],
    )
    write_tsv(
        output_dir / "exposure_summary.tsv",
        by_exposure,
        ["domain", "exposure", "runs", "fdr5_hits_total", "fdr10_hits_total", "runs_with_fdr5", "max_pm1"],
    )
    write_tsv(output_dir / "gene_recurrence.tsv", gene_rows, ["gene", "fdr5_run_count"])

    with (output_dir / "headline_summary.txt").open("w") as handle:
        handle.write(f"run_count\t{len(run_rows)}\n")
        handle.write(f"fdr5_hit_rows\t{len(hit_rows)}\n")
        handle.write(f"vascular_runs\t{sum(1 for r in run_rows if r['domain'] == 'vascular')}\n")
        handle.write(f"kidney_runs\t{sum(1 for r in run_rows if r['domain'] == 'kidney')}\n")
        if gene_rows:
            handle.write(f"top_gene\t{gene_rows[0]['gene']}\n")
            handle.write(f"top_gene_run_count\t{gene_rows[0]['fdr5_run_count']}\n")


if __name__ == "__main__":
    main()
