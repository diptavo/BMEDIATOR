from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


FEATURES = (
    "logit_p_m1",
    "p_m1",
    "p_m2",
    "p_m4",
    "p_m5",
    "m1_margin_max_alt",
    "m1_minus_m2",
    "m1_minus_m5",
    "m1_minus_m4",
    "log1p_nA",
    "log1p_nB",
    "log1p_nC",
    "rf_to_pp_z",
    "pp_to_outcome_z",
    "mediated_z",
    "direction_consistency_prob",
    "rf_to_pp_identifiable",
)

MODEL_NAME = "calibration_model"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fit and apply empirical calibration for BMEDIATOR P(M1).")
    parser.add_argument("--outdir", type=Path, required=True, help="Simulation output directory.")
    parser.add_argument(
        "--train-benchmark",
        default="calibration",
        choices=["classification", "calibration", "all"],
        help="Benchmark rows used to fit the calibration model.",
    )
    parser.add_argument("--holdout-mod", type=int, default=5, help="Replicates divisible by this value are held out.")
    parser.add_argument("--ridge", type=float, default=1.0, help="L2 penalty for logistic calibration.")
    parser.add_argument("--max-iter", type=int, default=50)
    parser.add_argument("--tol", type=float, default=1e-7)
    parser.add_argument("--model-name", default="calibration_model", help="Subdirectory under summary/ for calibration outputs.")
    parser.add_argument(
        "--empirical-min-selected",
        type=int,
        default=50,
        help="Minimum training selections used when estimating worst-cell empirical BFDR.",
    )
    parser.add_argument(
        "--empirical-z",
        type=float,
        default=1.64,
        help="Wilson upper-bound z value for conservative empirical BFDR; use 0 for plug-in FDR.",
    )
    return parser.parse_args()


def to_float(value, default: float = 0.0) -> float:
    try:
        x = float(value)
        if math.isfinite(x):
            return x
    except Exception:
        pass
    return default


def to_int(value, default: int = 0) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def clipped_prob(value: float) -> float:
    return min(1.0 - 1e-8, max(1e-8, value))


def logit(value: float) -> float:
    p = clipped_prob(value)
    return math.log(p / (1.0 - p))


def z_score(beta, se) -> float:
    b = to_float(beta, 0.0)
    s = to_float(se, 0.0)
    if s <= 0:
        return 0.0
    return min(80.0, abs(b / s))


def row_features(row: Dict[str, str]) -> List[float]:
    p1 = to_float(row.get("P_M1", row.get("p_m1", 0.0)), 0.0)
    p2 = to_float(row.get("P_M2", 0.0), 0.0)
    p4 = to_float(row.get("P_M4", 0.0), 0.0)
    p5 = to_float(row.get("P_M5", 0.0), 0.0)
    n_a = to_float(row.get("nA", row.get("nA_true", 0.0)), 0.0)
    n_b = to_float(row.get("nB", row.get("nB_true", 0.0)), 0.0)
    n_c = to_float(row.get("nC", row.get("nC_true", 0.0)), 0.0)
    identifiable = 1.0 if str(row.get("rf_to_pp_identifiable", "")).upper() == "YES" else 0.0
    return [
        logit(p1),
        p1,
        p2,
        p4,
        p5,
        p1 - max(p2, p4, p5),
        p1 - p2,
        p1 - p5,
        p1 - p4,
        math.log1p(max(0.0, n_a)),
        math.log1p(max(0.0, n_b)),
        math.log1p(max(0.0, n_c)),
        z_score(row.get("ivw_rf_to_pp_beta"), row.get("ivw_rf_to_pp_se")),
        z_score(row.get("ivw_pp_to_outcome_beta"), row.get("ivw_pp_to_outcome_se")),
        z_score(row.get("mediated_effect"), row.get("se_mediated")),
        to_float(row.get("direction_consistency_prob", 0.5), 0.5),
        identifiable,
    ]


def iter_metric_paths(outdir: Path) -> Iterable[Tuple[str, Path]]:
    for benchmark in ("classification", "calibration"):
        path = outdir / "summary" / benchmark / "protein_level_metrics.tsv"
        if path.exists():
            yield benchmark, path


def load_rows(outdir: Path) -> Tuple[List[Dict[str, str]], np.ndarray, np.ndarray, np.ndarray]:
    rows: List[Dict[str, str]] = []
    features: List[List[float]] = []
    y: List[int] = []
    holdout: List[bool] = []
    for benchmark, path in iter_metric_paths(outdir):
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                row = dict(row)
                row["benchmark"] = row.get("benchmark", benchmark) or benchmark
                rows.append(row)
                features.append(row_features(row))
                y.append(to_int(row.get("is_true_m1", 0), 0))
                holdout.append(False)
    return rows, np.asarray(features, dtype=float), np.asarray(y, dtype=float), np.asarray(holdout, dtype=bool)


def train_mask_for_rows(rows: Sequence[Dict[str, str]], train_benchmark: str, holdout_mod: int) -> Tuple[np.ndarray, np.ndarray]:
    train = []
    holdout = []
    for row in rows:
        benchmark_ok = train_benchmark == "all" or row.get("benchmark") == train_benchmark
        rep = to_int(row.get("replicate", 0), 0)
        is_holdout = holdout_mod > 0 and rep % holdout_mod == 0
        train.append(bool(benchmark_ok and not is_holdout))
        holdout.append(bool(is_holdout))
    return np.asarray(train, dtype=bool), np.asarray(holdout, dtype=bool)


def standardize(x: np.ndarray, train_mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.mean(x[train_mask], axis=0)
    sd = np.std(x[train_mask], axis=0)
    sd[sd < 1e-8] = 1.0
    z = (x - mean) / sd
    return z, mean, sd


def sigmoid(eta: np.ndarray) -> np.ndarray:
    eta = np.clip(eta, -40.0, 40.0)
    return 1.0 / (1.0 + np.exp(-eta))


def fit_logistic_newton(x: np.ndarray, y: np.ndarray, ridge: float, max_iter: int, tol: float) -> np.ndarray:
    design = np.column_stack([np.ones(x.shape[0]), x])
    beta = np.zeros(design.shape[1], dtype=float)
    penalty = np.eye(design.shape[1], dtype=float) * ridge
    penalty[0, 0] = 0.0
    prev_loss = float("inf")
    for _ in range(max_iter):
        eta = design @ beta
        p = sigmoid(eta)
        w = np.maximum(p * (1.0 - p), 1e-8)
        grad = design.T @ (p - y) + penalty @ beta
        hess = (design.T * w) @ design + penalty
        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(hess, grad, rcond=None)[0]
        candidate = beta - step
        p_cand = sigmoid(design @ candidate)
        loss = -float(np.sum(y * np.log(np.maximum(p_cand, 1e-12)) + (1.0 - y) * np.log(np.maximum(1.0 - p_cand, 1e-12))))
        loss += 0.5 * ridge * float(np.sum(candidate[1:] ** 2))
        if abs(prev_loss - loss) < tol * max(1.0, prev_loss):
            beta = candidate
            break
        if loss > prev_loss:
            beta = beta - 0.5 * step
        else:
            beta = candidate
            prev_loss = loss
    return beta


def predict_logistic(x: np.ndarray, beta: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(x.shape[0]), x])
    return sigmoid(design @ beta)


def auroc(y: np.ndarray, score: np.ndarray) -> float:
    positives = int(np.sum(y == 1))
    negatives = int(np.sum(y == 0))
    if positives == 0 or negatives == 0:
        return float("nan")
    order = np.argsort(score)
    ranks = np.empty(len(score), dtype=float)
    ranks[order] = np.arange(1, len(score) + 1)
    return float((np.sum(ranks[y == 1]) - positives * (positives + 1) / 2.0) / (positives * negatives))


def auprc(y: np.ndarray, score: np.ndarray) -> float:
    positives = int(np.sum(y == 1))
    if positives == 0:
        return float("nan")
    order = np.argsort(-score)
    y_sorted = y[order]
    tp = np.cumsum(y_sorted)
    precision = tp / np.arange(1, len(y_sorted) + 1)
    return float(np.sum(precision[y_sorted == 1]) / positives)


def brier(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean((p - y) ** 2))


def write_model(outdir: Path, beta: np.ndarray, mean: np.ndarray, sd: np.ndarray, args: argparse.Namespace) -> None:
    model_dir = outdir / "summary" / MODEL_NAME
    model_dir.mkdir(parents=True, exist_ok=True)
    with (model_dir / "bmediator_calibration_coefficients.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["term", "coefficient", "feature_mean", "feature_sd"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"term": "intercept", "coefficient": beta[0], "feature_mean": "", "feature_sd": ""})
        for idx, name in enumerate(FEATURES):
            writer.writerow({"term": name, "coefficient": beta[idx + 1], "feature_mean": mean[idx], "feature_sd": sd[idx]})
    with (model_dir / "bmediator_calibration_settings.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["setting", "value"], delimiter="\t")
        writer.writeheader()
        for key in ("train_benchmark", "holdout_mod", "ridge", "max_iter", "tol", "empirical_min_selected", "empirical_z"):
            writer.writerow({"setting": key, "value": getattr(args, key)})


def calibrated_bfdr(prob: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    order = np.argsort(-prob)
    rank = np.empty(len(prob), dtype=int)
    cum = np.empty(len(prob), dtype=float)
    running = np.cumsum(1.0 - prob[order]) / np.arange(1, len(prob) + 1)
    rank[order] = np.arange(1, len(prob) + 1)
    cum[order] = running
    return rank, cum


def wilson_upper(false_count: np.ndarray, n: np.ndarray, z: float) -> np.ndarray:
    n = np.maximum(n.astype(float), 1.0)
    p_hat = false_count.astype(float) / n
    if z <= 0:
        return p_hat
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p_hat + z2 / (2.0 * n)) / denom
    half_width = z * np.sqrt((p_hat * (1.0 - p_hat) / n) + (z2 / (4.0 * n * n))) / denom
    return np.minimum(1.0, center + half_width)


def empirical_worst_cell_bfdr(
    rows: Sequence[Dict[str, str]],
    y: np.ndarray,
    score: np.ndarray,
    train_mask: np.ndarray,
    min_selected: int,
    z: float,
) -> np.ndarray:
    worst = np.ones(len(rows), dtype=float)
    benchmarks = sorted(set(row.get("benchmark", "") for row in rows))
    for benchmark in benchmarks:
        target_idx = np.asarray([i for i, row in enumerate(rows) if row.get("benchmark", "") == benchmark], dtype=int)
        if len(target_idx) == 0:
            continue
        target_worst = np.zeros(len(target_idx), dtype=float)
        cells = sorted(set(row.get("cell", "") for row in rows if row.get("benchmark", "") == benchmark))
        for cell in cells:
            source_idx = np.asarray(
                [
                    i
                    for i, row in enumerate(rows)
                    if train_mask[i] and row.get("benchmark", "") == benchmark and row.get("cell", "") == cell
                ],
                dtype=int,
            )
            if len(source_idx) == 0:
                continue
            order = np.argsort(-score[source_idx])
            source_score = score[source_idx][order]
            source_false = (y[source_idx][order] == 0).astype(float)
            cumulative_false = np.cumsum(source_false)
            selected = np.searchsorted(-source_score, -score[target_idx], side="right")
            selected = np.clip(selected, max(1, min_selected), len(source_idx))
            upper = wilson_upper(cumulative_false[selected - 1], selected.astype(float), z)
            target_worst = np.maximum(target_worst, upper)
        worst[target_idx] = target_worst
    return worst


def empirical_cell_bfdr(
    rows: Sequence[Dict[str, str]],
    y: np.ndarray,
    score: np.ndarray,
    train_mask: np.ndarray,
    min_selected: int,
    z: float,
) -> np.ndarray:
    empirical = np.ones(len(rows), dtype=float)
    benchmarks = sorted(set(row.get("benchmark", "") for row in rows))
    for benchmark in benchmarks:
        cells = sorted(set(row.get("cell", "") for row in rows if row.get("benchmark", "") == benchmark))
        for cell in cells:
            target_idx = np.asarray(
                [i for i, row in enumerate(rows) if row.get("benchmark", "") == benchmark and row.get("cell", "") == cell],
                dtype=int,
            )
            source_idx = np.asarray(
                [
                    i
                    for i, row in enumerate(rows)
                    if train_mask[i] and row.get("benchmark", "") == benchmark and row.get("cell", "") == cell
                ],
                dtype=int,
            )
            if len(target_idx) == 0 or len(source_idx) == 0:
                continue
            order = np.argsort(-score[source_idx])
            source_score = score[source_idx][order]
            source_false = (y[source_idx][order] == 0).astype(float)
            cumulative_false = np.cumsum(source_false)
            selected = np.searchsorted(-source_score, -score[target_idx], side="right")
            selected = np.clip(selected, max(1, min_selected), len(source_idx))
            empirical[target_idx] = wilson_upper(cumulative_false[selected - 1], selected.astype(float), z)
    return empirical


def write_scores(
    outdir: Path,
    rows: Sequence[Dict[str, str]],
    prob: np.ndarray,
    cell_bfdr: np.ndarray,
    worst_cell_bfdr: np.ndarray,
) -> None:
    model_dir = outdir / "summary" / MODEL_NAME
    rank_global, bfdr_global = calibrated_bfdr(prob)
    per_benchmark_rank = np.zeros(len(rows), dtype=int)
    per_benchmark_bfdr = np.zeros(len(rows), dtype=float)
    for benchmark in sorted(set(row.get("benchmark", "") for row in rows)):
        idx = np.asarray([i for i, row in enumerate(rows) if row.get("benchmark", "") == benchmark], dtype=int)
        r, b = calibrated_bfdr(prob[idx])
        per_benchmark_rank[idx] = r
        per_benchmark_bfdr[idx] = b
    fields = [
        "benchmark",
        "cell",
        "replicate",
        "protein_id",
        "true_scenario",
        "is_true_m1",
        "P_M1",
        "P_M2",
        "P_M5",
        "p_m1_calibrated",
        "calibrated_rank_global",
        "calibrated_bfdr_global",
        "calibrated_rank_benchmark",
        "calibrated_bfdr_benchmark",
        "calibrated_bfdr_cell_empirical",
        "calibrated_bfdr_worst_cell",
        "estimated_bfdr_raw",
        "nA",
        "nB",
        "nC",
    ]
    with (model_dir / "bmediator_calibrated_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for i, row in enumerate(rows):
            writer.writerow(
                {
                    "benchmark": row.get("benchmark", ""),
                    "cell": row.get("cell", ""),
                    "replicate": row.get("replicate", ""),
                    "protein_id": row.get("protein_id", row.get("Protein", "")),
                    "true_scenario": row.get("true_scenario", ""),
                    "is_true_m1": row.get("is_true_m1", ""),
                    "P_M1": row.get("P_M1", row.get("p_m1", "")),
                    "P_M2": row.get("P_M2", ""),
                    "P_M5": row.get("P_M5", ""),
                    "p_m1_calibrated": prob[i],
                    "calibrated_rank_global": rank_global[i],
                    "calibrated_bfdr_global": bfdr_global[i],
                    "calibrated_rank_benchmark": per_benchmark_rank[i],
                    "calibrated_bfdr_benchmark": per_benchmark_bfdr[i],
                    "calibrated_bfdr_cell_empirical": cell_bfdr[i],
                    "calibrated_bfdr_worst_cell": worst_cell_bfdr[i],
                    "estimated_bfdr_raw": row.get("estimated_bfdr", ""),
                    "nA": row.get("nA", ""),
                    "nB": row.get("nB", ""),
                    "nC": row.get("nC", ""),
                }
            )


def group_indices(rows: Sequence[Dict[str, str]], masks: Dict[str, np.ndarray]) -> Dict[Tuple[str, str, str], np.ndarray]:
    groups = {}
    for split_name, split_mask in masks.items():
        by_key = defaultdict(list)
        for i, row in enumerate(rows):
            if not split_mask[i]:
                continue
            by_key[(split_name, row.get("benchmark", ""), row.get("cell", ""))].append(i)
        for key, values in by_key.items():
            groups[key] = np.asarray(values, dtype=int)
    return groups


def write_performance(outdir: Path, rows: Sequence[Dict[str, str]], y: np.ndarray, raw: np.ndarray, prob: np.ndarray, train_mask: np.ndarray, holdout_mask: np.ndarray) -> None:
    model_dir = outdir / "summary" / MODEL_NAME
    masks = {
        "train": train_mask,
        "holdout": holdout_mask,
        "all": np.ones(len(rows), dtype=bool),
    }
    groups = group_indices(rows, masks)
    fields = ["split", "benchmark", "cell", "n", "prevalence", "score", "auroc", "auprc", "brier", "mean_prob_true", "mean_prob_false"]
    with (model_dir / "bmediator_calibration_performance.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(groups):
            idx = groups[key]
            for score_name, score in (("raw_P_M1", raw), ("calibrated_P_M1", prob)):
                yy = y[idx]
                ss = score[idx]
                writer.writerow(
                    {
                        "split": key[0],
                        "benchmark": key[1],
                        "cell": key[2],
                        "n": len(idx),
                        "prevalence": float(np.mean(yy)) if len(idx) else float("nan"),
                        "score": score_name,
                        "auroc": auroc(yy, ss),
                        "auprc": auprc(yy, ss),
                        "brier": brier(yy, ss),
                        "mean_prob_true": float(np.mean(ss[yy == 1])) if np.any(yy == 1) else float("nan"),
                        "mean_prob_false": float(np.mean(ss[yy == 0])) if np.any(yy == 0) else float("nan"),
                    }
                )


def write_bins(outdir: Path, rows: Sequence[Dict[str, str]], y: np.ndarray, raw: np.ndarray, prob: np.ndarray, holdout_mask: np.ndarray) -> None:
    model_dir = outdir / "summary" / MODEL_NAME
    fields = ["split", "benchmark", "cell", "score", "bin", "n", "mean_pred", "observed_m1"]
    bins = np.linspace(0.0, 1.0, 11)
    with (model_dir / "bmediator_calibration_bins.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for split_name, split_mask in (("holdout", holdout_mask), ("all", np.ones(len(rows), dtype=bool))):
            for benchmark in sorted(set(row.get("benchmark", "") for row in rows)):
                for cell in sorted(set(row.get("cell", "") for row in rows if row.get("benchmark", "") == benchmark)):
                    idx = np.asarray([i for i, row in enumerate(rows) if split_mask[i] and row.get("benchmark", "") == benchmark and row.get("cell", "") == cell], dtype=int)
                    if len(idx) == 0:
                        continue
                    for score_name, score in (("raw_P_M1", raw), ("calibrated_P_M1", prob)):
                        vals = score[idx]
                        yy = y[idx]
                        for b in range(10):
                            lo = bins[b]
                            hi = bins[b + 1]
                            if b == 9:
                                keep = (vals >= lo) & (vals <= hi)
                            else:
                                keep = (vals >= lo) & (vals < hi)
                            if not np.any(keep):
                                continue
                            writer.writerow(
                                {
                                    "split": split_name,
                                    "benchmark": benchmark,
                                    "cell": cell,
                                    "score": score_name,
                                    "bin": b,
                                    "n": int(np.sum(keep)),
                                    "mean_pred": float(np.mean(vals[keep])),
                                    "observed_m1": float(np.mean(yy[keep])),
                                }
                            )


def write_fdr(
    outdir: Path,
    rows: Sequence[Dict[str, str]],
    y: np.ndarray,
    raw_bfdr: np.ndarray,
    calibrated_bfdr_benchmark: np.ndarray,
    cell_bfdr: np.ndarray,
    worst_cell_bfdr: np.ndarray,
    holdout_mask: np.ndarray,
) -> None:
    model_dir = outdir / "summary" / MODEL_NAME
    fields = ["split", "benchmark", "cell", "threshold", "method", "n_selected", "true_positives", "false_positives", "empirical_fdr", "power"]
    thresholds = (0.01, 0.05, 0.10)
    with (model_dir / "bmediator_calibrated_fdr_power.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for split_name, split_mask in (("holdout", holdout_mask), ("all", np.ones(len(rows), dtype=bool))):
            for benchmark in sorted(set(row.get("benchmark", "") for row in rows)):
                for cell in sorted(set(row.get("cell", "") for row in rows if row.get("benchmark", "") == benchmark)):
                    idx = np.asarray([i for i, row in enumerate(rows) if split_mask[i] and row.get("benchmark", "") == benchmark and row.get("cell", "") == cell], dtype=int)
                    if len(idx) == 0:
                        continue
                    total_true = int(np.sum(y[idx] == 1))
                    for method, bfdr in (
                        ("raw_estimated_bfdr", raw_bfdr),
                        ("calibrated_bfdr", calibrated_bfdr_benchmark),
                        ("calibrated_bfdr_cell_empirical", cell_bfdr),
                        ("calibrated_bfdr_worst_cell", worst_cell_bfdr),
                    ):
                        vals = bfdr[idx]
                        for threshold in thresholds:
                            selected = np.isfinite(vals) & (vals <= threshold)
                            n_selected = int(np.sum(selected))
                            tp = int(np.sum(y[idx][selected] == 1))
                            fp = n_selected - tp
                            writer.writerow(
                                {
                                    "split": split_name,
                                    "benchmark": benchmark,
                                    "cell": cell,
                                    "threshold": threshold,
                                    "method": method,
                                    "n_selected": n_selected,
                                    "true_positives": tp,
                                    "false_positives": fp,
                                    "empirical_fdr": (fp / n_selected) if n_selected else "",
                                    "power": (tp / total_true) if total_true else "",
                                }
                            )


def main() -> None:
    global MODEL_NAME
    args = parse_args()
    MODEL_NAME = args.model_name
    rows, x, y, _ = load_rows(args.outdir)
    if not rows:
        raise SystemExit("No protein_level_metrics.tsv files found.")
    train_mask, holdout_mask = train_mask_for_rows(rows, args.train_benchmark, args.holdout_mod)
    if int(np.sum(train_mask)) == 0:
        raise SystemExit("No training rows selected.")
    x_std, mean, sd = standardize(x, train_mask)
    beta = fit_logistic_newton(x_std[train_mask], y[train_mask], args.ridge, args.max_iter, args.tol)
    prob = predict_logistic(x_std, beta)
    raw = np.asarray([to_float(row.get("P_M1", row.get("p_m1", 0.0)), 0.0) for row in rows], dtype=float)
    raw_bfdr = np.asarray([to_float(row.get("estimated_bfdr", ""), float("nan")) for row in rows], dtype=float)
    _, calibrated_bfdr_benchmark = calibrated_bfdr(prob)
    # Replace with per-benchmark BFDR for reporting.
    calibrated_bfdr_benchmark = np.zeros(len(rows), dtype=float)
    for benchmark in sorted(set(row.get("benchmark", "") for row in rows)):
        idx = np.asarray([i for i, row in enumerate(rows) if row.get("benchmark", "") == benchmark], dtype=int)
        _, b = calibrated_bfdr(prob[idx])
        calibrated_bfdr_benchmark[idx] = b
    cell_bfdr = empirical_cell_bfdr(rows, y, prob, train_mask, args.empirical_min_selected, args.empirical_z)
    worst_cell_bfdr = empirical_worst_cell_bfdr(rows, y, prob, train_mask, args.empirical_min_selected, args.empirical_z)
    write_model(args.outdir, beta, mean, sd, args)
    write_scores(args.outdir, rows, prob, cell_bfdr, worst_cell_bfdr)
    write_performance(args.outdir, rows, y, raw, prob, train_mask, holdout_mask)
    write_bins(args.outdir, rows, y, raw, prob, holdout_mask)
    write_fdr(args.outdir, rows, y, raw_bfdr, calibrated_bfdr_benchmark, cell_bfdr, worst_cell_bfdr, holdout_mask)
    print(f"Wrote calibration outputs under {args.outdir / 'summary' / MODEL_NAME}")


if __name__ == "__main__":
    main()
