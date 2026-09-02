from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from benchmark_lib import read_results, read_table, read_truth


METHODS = (
    "bmediator",
    "two_step_product",
    "mr_ivw_pqtl_outcome",
    "mr_egger_pqtl_outcome",
    "smr_top_cis",
    "coloc_proxy",
    "mvmr_rf_pqtl",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run competitor-method summaries on BMEDIATOR simulation outputs.")
    parser.add_argument("--outdir", type=Path, required=True, help="Simulation output directory from sim/run_benchmark.py.")
    parser.add_argument("--benchmark", choices=["classification", "calibration", "all"], default="all")
    parser.add_argument("--rf-p-threshold", type=float, default=5e-6)
    parser.add_argument("--cis-p-threshold", type=float, default=5e-6)
    parser.add_argument("--cis-window-bp", type=int, default=1_000_000)
    return parser.parse_args()


def to_float(value, default: float = float("nan")) -> float:
    try:
        return float(value)
    except Exception:
        return default


def two_sided_p_from_z(z: float) -> float:
    if not math.isfinite(z):
        return 1.0
    return max(0.0, min(1.0, math.erfc(abs(z) / math.sqrt(2.0))))


def score_from_p(p_value: float) -> float:
    if not math.isfinite(p_value):
        return 0.0
    return -math.log10(max(p_value, 1e-300))


def weighted_slope_z(x: Sequence[float], y: Sequence[float], se_y: Sequence[float]) -> Tuple[float, float, float]:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    se_arr = np.asarray(se_y, dtype=float)
    keep = np.isfinite(x_arr) & np.isfinite(y_arr) & np.isfinite(se_arr) & (se_arr > 0) & (np.abs(x_arr) > 1e-12)
    x_arr = x_arr[keep]
    y_arr = y_arr[keep]
    se_arr = se_arr[keep]
    if len(x_arr) < 1:
        return float("nan"), float("nan"), 0.0
    w = 1.0 / (se_arr ** 2)
    denom = float(np.sum(w * x_arr * x_arr))
    if denom <= 0:
        return float("nan"), float("nan"), 0.0
    beta = float(np.sum(w * x_arr * y_arr) / denom)
    se = math.sqrt(1.0 / denom)
    return beta, se, beta / se if se > 0 else 0.0


def weighted_regression_z(xs: Sequence[Sequence[float]], y: Sequence[float], se_y: Sequence[float], target_col: int) -> Tuple[float, float, float]:
    x_mat = np.asarray(xs, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    se_arr = np.asarray(se_y, dtype=float)
    if x_mat.ndim == 1:
        x_mat = x_mat.reshape((-1, 1))
    keep = np.isfinite(y_arr) & np.isfinite(se_arr) & (se_arr > 0)
    keep &= np.all(np.isfinite(x_mat), axis=1)
    keep &= np.any(np.abs(x_mat) > 1e-12, axis=1)
    x_mat = x_mat[keep, :]
    y_arr = y_arr[keep]
    se_arr = se_arr[keep]
    if len(y_arr) <= x_mat.shape[1] or target_col >= x_mat.shape[1]:
        return float("nan"), float("nan"), 0.0
    w_sqrt = 1.0 / se_arr
    xw = x_mat * w_sqrt[:, None]
    yw = y_arr * w_sqrt
    try:
        coef, *_ = np.linalg.lstsq(xw, yw, rcond=None)
        xtwx_inv = np.linalg.inv(xw.T @ xw)
    except np.linalg.LinAlgError:
        return float("nan"), float("nan"), 0.0
    beta = float(coef[target_col])
    se = math.sqrt(max(float(xtwx_inv[target_col, target_col]), 0.0))
    return beta, se, beta / se if se > 0 else 0.0


def load_sumstats(rep_dir: Path) -> Tuple[Dict, Dict, Dict[str, List[Dict]], Dict]:
    rf = {row["SNP"]: row for row in read_table(rep_dir / "rf_sumstat.txt", " ")}
    outcome = {row["SNP"]: row for row in read_table(rep_dir / "cancer_sumstat.txt", " ")}
    pqtl_by_protein: Dict[str, List[Dict]] = {}
    for row in read_table(rep_dir / "pqtl_sumstat.txt", " "):
        pqtl_by_protein.setdefault(row["PROTEIN"], []).append(row)
    info = {row["PROTEIN"]: row for row in read_table(rep_dir / "protein_info.txt", " ")}
    return rf, outcome, pqtl_by_protein, info


def protein_rows(
    protein: str,
    rf: Dict,
    outcome: Dict,
    pqtl_by_protein: Dict[str, List[Dict]],
    info: Dict,
    rf_p_threshold: float,
    cis_p_threshold: float,
    cis_window_bp: int,
) -> Tuple[List[Dict], List[Dict], List[Dict]]:
    pinfo = info[protein]
    chrom = str(pinfo["CHR"])
    start = int(float(pinfo["START"])) - cis_window_bp
    end = int(float(pinfo["END"])) + cis_window_bp
    all_rows = []
    rf_inst = []
    cis_inst = []
    for prow in pqtl_by_protein.get(protein, []):
        snp = prow["SNP"]
        if snp not in rf or snp not in outcome:
            continue
        rrow = rf[snp]
        orow = outcome[snp]
        bp = int(float(prow["BP"]))
        is_cis = str(prow["CHR"]) == chrom and start <= bp <= end
        row = {
            "snp": snp,
            "rf_beta": to_float(rrow["BETA"]),
            "rf_se": to_float(rrow["SE"]),
            "rf_p": to_float(rrow["P"]),
            "pqtl_beta": to_float(prow["BETA"]),
            "pqtl_se": to_float(prow["SE"]),
            "pqtl_p": to_float(prow["P"]),
            "outcome_beta": to_float(orow["BETA"]),
            "outcome_se": to_float(orow["SE"]),
            "outcome_p": to_float(orow["P"]),
            "is_cis": is_cis,
        }
        all_rows.append(row)
        if row["rf_p"] <= rf_p_threshold:
            rf_inst.append(row)
        if is_cis and row["pqtl_p"] <= cis_p_threshold:
            cis_inst.append(row)
    return all_rows, rf_inst, cis_inst


def method_score_rows(
    protein: str,
    true_scenario: str,
    all_rows: List[Dict],
    rf_inst: List[Dict],
    cis_inst: List[Dict],
    bmed_p_m1: float | None,
) -> List[Dict]:
    out = []
    is_m1 = int(true_scenario == "M1")

    if bmed_p_m1 is not None and math.isfinite(bmed_p_m1):
        out.append(
            {
                "method": "bmediator",
                "protein_id": protein,
                "true_scenario": true_scenario,
                "is_true_m1": is_m1,
                "score": bmed_p_m1,
                "p_value": 1.0 - bmed_p_m1,
                "n_instruments": len(all_rows),
            }
        )

    beta1, se1, z1 = weighted_slope_z(
        [row["rf_beta"] for row in rf_inst],
        [row["pqtl_beta"] for row in rf_inst],
        [row["pqtl_se"] for row in rf_inst],
    )
    beta2, se2, z2 = weighted_slope_z(
        [row["pqtl_beta"] for row in cis_inst],
        [row["outcome_beta"] for row in cis_inst],
        [row["outcome_se"] for row in cis_inst],
    )
    two_step_z = min(abs(z1), abs(z2)) if math.isfinite(z1) and math.isfinite(z2) else 0.0
    two_step_p = two_sided_p_from_z(two_step_z)
    out.append(
        {
            "method": "two_step_product",
            "protein_id": protein,
            "true_scenario": true_scenario,
            "is_true_m1": is_m1,
            "score": score_from_p(two_step_p),
            "p_value": two_step_p,
            "n_instruments": len(rf_inst) + len(cis_inst),
        }
    )

    ivw_p = two_sided_p_from_z(z2)
    out.append(
        {
            "method": "mr_ivw_pqtl_outcome",
            "protein_id": protein,
            "true_scenario": true_scenario,
            "is_true_m1": is_m1,
            "score": score_from_p(ivw_p),
            "p_value": ivw_p,
            "n_instruments": len(cis_inst),
        }
    )

    egger_x = [[1.0, row["pqtl_beta"]] for row in cis_inst]
    _, _, egger_z = weighted_regression_z(egger_x, [row["outcome_beta"] for row in cis_inst], [row["outcome_se"] for row in cis_inst], 1)
    egger_p = two_sided_p_from_z(egger_z)
    out.append(
        {
            "method": "mr_egger_pqtl_outcome",
            "protein_id": protein,
            "true_scenario": true_scenario,
            "is_true_m1": is_m1,
            "score": score_from_p(egger_p),
            "p_value": egger_p,
            "n_instruments": len(cis_inst),
        }
    )

    if cis_inst:
        top = min(cis_inst, key=lambda row: row["pqtl_p"])
        denom = top["pqtl_beta"]
        if abs(denom) > 1e-12:
            smr_beta = top["outcome_beta"] / denom
            smr_var = (top["outcome_se"] ** 2) / (denom ** 2)
            smr_var += ((top["outcome_beta"] ** 2) * (top["pqtl_se"] ** 2)) / (denom ** 4)
            smr_z = smr_beta / math.sqrt(max(smr_var, 1e-300))
        else:
            smr_z = 0.0
    else:
        smr_z = 0.0
    smr_p = two_sided_p_from_z(smr_z)
    out.append(
        {
            "method": "smr_top_cis",
            "protein_id": protein,
            "true_scenario": true_scenario,
            "is_true_m1": is_m1,
            "score": score_from_p(smr_p),
            "p_value": smr_p,
            "n_instruments": 1 if cis_inst else 0,
        }
    )

    coloc_z = 0.0
    for row in cis_inst:
        z_pqtl = abs(row["pqtl_beta"] / row["pqtl_se"]) if row["pqtl_se"] > 0 else 0.0
        z_outcome = abs(row["outcome_beta"] / row["outcome_se"]) if row["outcome_se"] > 0 else 0.0
        coloc_z = max(coloc_z, min(z_pqtl, z_outcome))
    coloc_p = two_sided_p_from_z(coloc_z)
    out.append(
        {
            "method": "coloc_proxy",
            "protein_id": protein,
            "true_scenario": true_scenario,
            "is_true_m1": is_m1,
            "score": score_from_p(coloc_p),
            "p_value": coloc_p,
            "n_instruments": len(cis_inst),
        }
    )

    mvmr_x = [[row["rf_beta"], row["pqtl_beta"]] for row in all_rows]
    _, _, mvmr_z = weighted_regression_z(mvmr_x, [row["outcome_beta"] for row in all_rows], [row["outcome_se"] for row in all_rows], 1)
    mvmr_p = two_sided_p_from_z(mvmr_z)
    out.append(
        {
            "method": "mvmr_rf_pqtl",
            "protein_id": protein,
            "true_scenario": true_scenario,
            "is_true_m1": is_m1,
            "score": score_from_p(mvmr_p),
            "p_value": mvmr_p,
            "n_instruments": len(all_rows),
        }
    )
    return out


def auroc(rows: List[Dict]) -> float | None:
    positives = sum(row["is_true_m1"] for row in rows)
    negatives = len(rows) - positives
    if positives == 0 or negatives == 0:
        return None
    sorted_rows = sorted(rows, key=lambda row: row["score"])
    rank_sum = 0.0
    for rank, row in enumerate(sorted_rows, start=1):
        if row["is_true_m1"]:
            rank_sum += rank
    return (rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def average_precision(rows: List[Dict]) -> float | None:
    positives = sum(row["is_true_m1"] for row in rows)
    if positives == 0:
        return None
    sorted_rows = sorted(rows, key=lambda row: row["score"], reverse=True)
    tp = 0
    precision_sum = 0.0
    for idx, row in enumerate(sorted_rows, start=1):
        if row["is_true_m1"]:
            tp += 1
            precision_sum += tp / idx
    return precision_sum / positives


def summarize_methods(rows: List[Dict]) -> List[Dict]:
    summary = []
    keys = sorted(set((row["benchmark"], row["cell"], row["method"]) for row in rows))
    for benchmark, cell, method in keys:
        subset = [row for row in rows if row["benchmark"] == benchmark and row["cell"] == cell and row["method"] == method]
        selected = [row for row in subset if row["p_value"] <= 0.05]
        true_m1 = sum(row["is_true_m1"] for row in subset)
        tp = sum(row["is_true_m1"] for row in selected)
        fp = len(selected) - tp
        sorted_subset = sorted(subset, key=lambda row: row["score"], reverse=True)
        top_rows = []
        for top_n in (10, 25, 50):
            cur = sorted_subset[: min(top_n, len(sorted_subset))]
            if cur:
                top_rows.append(
                    {
                        "top_n": top_n,
                        "top_empirical_fdr": 1.0 - float(np.mean([row["is_true_m1"] for row in cur])),
                        "top_power": sum(row["is_true_m1"] for row in cur) / true_m1 if true_m1 else None,
                    }
                )
        base = {
            "benchmark": benchmark,
            "cell": cell,
            "method": method,
            "n": len(subset),
            "n_true_m1": true_m1,
            "auroc": auroc(subset),
            "auprc": average_precision(subset),
            "n_selected_p05": len(selected),
            "empirical_fdr_p05": fp / len(selected) if selected else None,
            "power_p05": tp / true_m1 if true_m1 else None,
        }
        for top in top_rows:
            row = dict(base)
            row.update(top)
            summary.append(row)
    return summary


def write_rows(path: Path, rows: List[Dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_compact_scores(path: Path) -> List[Dict]:
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for row in rows:
        for key in ("is_true_m1", "n_instruments", "replicate"):
            if key in row:
                row[key] = int(float(row[key]))
        for key in ("score", "p_value"):
            if key in row:
                row[key] = to_float(row[key])
    return rows


def main() -> None:
    args = parse_args()
    benchmarks = ("classification", "calibration") if args.benchmark == "all" else (args.benchmark,)
    protein_scores: List[Dict] = []

    for benchmark in benchmarks:
        for rep_dir in sorted((args.outdir / benchmark).glob("*/rep_*")):
            if not rep_dir.is_dir() or not (rep_dir / "truth.tsv").exists():
                compact_scores = rep_dir / "competitor_scores.tsv"
                if compact_scores.exists():
                    protein_scores.extend(read_compact_scores(compact_scores))
                continue
            compact_scores = rep_dir / "competitor_scores.tsv"
            if compact_scores.exists():
                protein_scores.extend(read_compact_scores(compact_scores))
                continue
            rf, outcome, pqtl_by_protein, info = load_sumstats(rep_dir)
            truth_rows = read_truth(rep_dir / "truth.tsv")
            bmed_p: Dict[str, float] = {}
            mediation = rep_dir / "bmediator.mediation"
            if mediation.exists():
                for row in read_results(mediation):
                    protein = str(row.get("Protein", row.get("protein_id", "")))
                    if protein:
                        bmed_p[protein] = to_float(row.get("P_M1", row.get("prob_M1")))

            for truth in truth_rows:
                protein = truth["protein_id"]
                all_rows, rf_inst, cis_inst = protein_rows(
                    protein,
                    rf,
                    outcome,
                    pqtl_by_protein,
                    info,
                    args.rf_p_threshold,
                    args.cis_p_threshold,
                    args.cis_window_bp,
                )
                rows = method_score_rows(protein, truth["true_scenario"], all_rows, rf_inst, cis_inst, bmed_p.get(protein))
                for row in rows:
                    row["benchmark"] = benchmark
                    row["cell"] = rep_dir.parent.name
                    row["replicate"] = int(rep_dir.name.split("_")[-1])
                    protein_scores.append(row)

    summary_dir = args.outdir / "summary" / "competitors"
    write_rows(summary_dir / "protein_level_competitor_scores.tsv", protein_scores)
    write_rows(summary_dir / "method_summary.tsv", summarize_methods(protein_scores))
    print(f"Wrote competitor summaries under {summary_dir}")


if __name__ == "__main__":
    main()
