#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import numpy as np


def finite_float(value):
    try:
        x = float(value)
    except ValueError:
        return None
    if math.isnan(x) or math.isinf(x):
        return None
    return x


def soft_threshold(x, lam):
    if x > lam:
        return x - lam
    if x < -lam:
        return x + lam
    return 0.0


def read_top_proteins(path, top_n):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            p_m1 = finite_float(row["P_M1"])
            rf_beta = finite_float(row["ivw_rf_to_pp_beta"])
            rf_se = finite_float(row["ivw_rf_to_pp_se"])
            out_beta = finite_float(row["ivw_pp_to_outcome_beta"])
            out_se = finite_float(row["ivw_pp_to_outcome_se"])
            if p_m1 is None or rf_beta is None or rf_se in (None, 0.0) or out_beta is None or out_se in (None, 0.0):
                continue
            row["_P_M1"] = p_m1
            row["_rf_beta"] = rf_beta
            row["_rf_se"] = rf_se
            row["_out_beta"] = out_beta
            row["_out_se"] = out_se
            row["_rf_z"] = rf_beta / rf_se
            row["_out_z"] = out_beta / out_se
            rows.append(row)
    rows.sort(key=lambda r: r["_P_M1"], reverse=True)
    return rows[:top_n]


def load_residual_matrix(path, protein_ids):
    protein_ids = list(protein_ids)
    wanted = set(protein_ids)
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_by_protein = {}
        ordered = []
        for idx, name in enumerate(header[1:], start=1):
            if name in wanted:
                idx_by_protein[name] = idx
                ordered.append(name)
        if not ordered:
            raise RuntimeError("none of the selected proteins were found in the residual matrix")

        data = [[] for _ in ordered]
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= 1:
                continue
            for j, protein in enumerate(ordered):
                idx = idx_by_protein[protein]
                if idx >= len(fields):
                    data[j].append(np.nan)
                    continue
                try:
                    data[j].append(float(fields[idx]))
                except ValueError:
                    data[j].append(np.nan)

    X = np.asarray(data, dtype=float).T
    keep = ~np.any(~np.isfinite(X), axis=1)
    X = X[keep]
    X = X - X.mean(axis=0, keepdims=True)
    scale = X.std(axis=0, ddof=1, keepdims=True)
    scale[scale == 0] = 1.0
    X = X / scale
    R = np.corrcoef(X, rowvar=False)
    return ordered, X.shape[0], R


def weighted_lasso(R, y, weights, lambdas, max_iter=5000, tol=1e-8):
    p = len(y)
    best = None
    path = []
    n_obs = max(p, 2)
    for lam in lambdas:
        beta = np.zeros(p, dtype=float)
        for _ in range(max_iter):
            old = beta.copy()
            for j in range(p):
                rho = y[j] - (R[j, :] @ beta - R[j, j] * beta[j])
                beta[j] = soft_threshold(rho, lam * weights[j]) / max(R[j, j], 1e-8)
            if np.max(np.abs(beta - old)) < tol:
                break
        resid = y - R @ beta
        rss = float(resid @ resid)
        k = int(np.sum(np.abs(beta) > 1e-8))
        ebic = n_obs * math.log(max(rss, 1e-12) / n_obs + 1e-12) + k * math.log(n_obs)
        path.append({"lambda": lam, "beta": beta.copy(), "ebic": ebic, "rss": rss, "k": k})
        if best is None or ebic < best["ebic"]:
            best = {"lambda": lam, "beta": beta.copy(), "ebic": ebic, "rss": rss, "k": k}
    return best, path


def refit_support(R, marginal_beta, support):
    p = len(marginal_beta)
    beta = np.zeros(p, dtype=float)
    if not np.any(support):
        return beta
    idx = np.where(support)[0]
    Rsub = R[np.ix_(idx, idx)] + 1e-4 * np.eye(len(idx))
    bsub = marginal_beta[idx]
    beta[idx] = np.linalg.solve(Rsub, bsub)
    return beta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mediation", required=True)
    ap.add_argument("--residuals", required=True)
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    top_rows = read_top_proteins(args.mediation, args.top_n)
    selected_ids = [row["Protein"] for row in top_rows]

    ordered_ids, n_samples, R = load_residual_matrix(args.residuals, selected_ids)
    row_by_id = {row["Protein"]: row for row in top_rows}
    rows = [row_by_id[p] for p in ordered_ids if p in row_by_id]

    rf_z = np.asarray([abs(row["_rf_z"]) for row in rows], dtype=float)
    out_z = np.asarray([row["_out_z"] for row in rows], dtype=float)
    out_beta = np.asarray([row["_out_beta"] for row in rows], dtype=float)
    rf_beta = np.asarray([row["_rf_beta"] for row in rows], dtype=float)

    y = out_z / max(np.max(np.abs(out_z)), 1.0)
    weights = 1.0 / (rf_z + 1e-3)
    lam_max = float(np.max(np.abs(y) / weights))
    lambdas = np.geomspace(lam_max, max(lam_max * 1e-3, 1e-5), num=60)

    best, path = weighted_lasso(R, y, weights, lambdas)
    if best["k"] == 0:
        nonzero = [item for item in path if item["k"] > 0]
        if nonzero:
            best = min(nonzero, key=lambda item: item["ebic"])
    support = np.abs(best["beta"]) > 1e-8
    joint_beta = refit_support(R, out_beta, support)
    mediated_score = rf_beta * joint_beta

    out_path = Path(args.out_prefix + ".joint.tsv")
    with open(out_path, "w") as out:
        out.write(
            "Protein\tGene\tP_M1\tivw_rf_to_pp_beta\tivw_rf_to_pp_se\tivw_rf_to_pp_p\t"
            "ivw_pp_to_outcome_beta\tivw_pp_to_outcome_se\tivw_pp_to_outcome_p\t"
            "rf_z\tout_z\tlasso_coef\tjoint_pp_to_outcome_beta\tjoint_mediated_score\tselected\n"
        )
        ranked = []
        for i, row in enumerate(rows):
            ranked.append((abs(mediated_score[i]), i))
        for _, i in sorted(ranked, reverse=True):
            row = rows[i]
            out.write(
                f'{row["Protein"]}\t{row["Gene"]}\t{row["_P_M1"]:.6g}\t'
                f'{row["_rf_beta"]:.6g}\t{row["_rf_se"]:.6g}\t{row["ivw_rf_to_pp_p"]}\t'
                f'{row["_out_beta"]:.6g}\t{row["_out_se"]:.6g}\t{row["ivw_pp_to_outcome_p"]}\t'
                f'{row["_rf_z"]:.6g}\t{row["_out_z"]:.6g}\t{best["beta"][i]:.6g}\t'
                f'{joint_beta[i]:.6g}\t{mediated_score[i]:.6g}\t{"YES" if support[i] else "NO"}\n'
            )

    corr_path = Path(args.out_prefix + ".corr.tsv")
    with open(corr_path, "w") as out:
        out.write("Protein\t" + "\t".join(ordered_ids) + "\n")
        for i, protein in enumerate(ordered_ids):
            out.write(protein + "\t" + "\t".join(f"{x:.6g}" for x in R[i]) + "\n")

    print(f"input_top_n\t{args.top_n}")
    print(f"proteins_in_joint_model\t{len(rows)}")
    print(f"aric_samples_used\t{n_samples}")
    print(f"best_lambda\t{best['lambda']:.6g}")
    print(f"selected_proteins\t{int(np.sum(support))}")
    print(f"joint_table\t{out_path}")
    print(f"corr_table\t{corr_path}")


if __name__ == "__main__":
    main()
