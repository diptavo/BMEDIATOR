from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def safe_float(value: str) -> float | None:
    try:
        x = float(value)
    except Exception:
        return None
    if math.isnan(x):
        return None
    return x


def sign_label(value: float | None) -> str:
    if value is None:
        return "NA"
    if value > 0:
        return "POS"
    if value < 0:
        return "NEG"
    return "ZERO"


def agreement_category(
    naive_support: bool,
    bmediator_support: bool,
    sign_agreement: str,
) -> str:
    if naive_support and bmediator_support:
        if sign_agreement == "YES":
            return "concordant_supported"
        if sign_agreement == "NO":
            return "discordant_sign"
        return "supported_sign_na"
    if naive_support and not bmediator_support:
        return "naive_only"
    if (not naive_support) and bmediator_support:
        return "bmediator_only"
    return "concordant_null"


def completed_mediation_files(results_dir: Path) -> list[Path]:
    return sorted(results_dir.rglob("*.mediation"))


def summarize_file(
    med_path: Path,
    outcome: str,
    pm1_threshold: float,
    naive_p_threshold: float,
) -> list[dict[str, str]]:
    exposure = med_path.parent.name
    rows_out: list[dict[str, str]] = []
    with med_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ivw1_beta = safe_float(row.get("ivw_rf_to_pp_beta", ""))
            ivw1_p = safe_float(row.get("ivw_rf_to_pp_p", ""))
            ivw2_beta = safe_float(row.get("ivw_pp_to_outcome_beta", ""))
            ivw2_p = safe_float(row.get("ivw_pp_to_outcome_p", ""))
            pm1 = safe_float(row.get("P_M1", "")) or 0.0
            mediated_effect = safe_float(row.get("mediated_effect", ""))

            naive_indirect = None
            if ivw1_beta is not None and ivw2_beta is not None:
                naive_indirect = ivw1_beta * ivw2_beta

            naive_support = (
                ivw1_p is not None
                and ivw2_p is not None
                and ivw1_p < naive_p_threshold
                and ivw2_p < naive_p_threshold
            )
            bmediator_support = pm1 >= pm1_threshold

            naive_sign = sign_label(naive_indirect)
            bmediator_sign = sign_label(mediated_effect)
            if naive_sign in {"POS", "NEG"} and bmediator_sign in {"POS", "NEG"}:
                sign_agreement = "YES" if naive_sign == bmediator_sign else "NO"
            else:
                sign_agreement = "NA"

            threshold_agreement = (
                "YES" if naive_support == bmediator_support else "NO"
            )

            rows_out.append(
                {
                    "exposure": exposure,
                    "outcome": outcome,
                    "gene": row.get("Gene", ""),
                    "protein": row.get("Protein", ""),
                    "P_M1": f"{pm1:.6f}",
                    "pm1_support": "YES" if bmediator_support else "NO",
                    "selected_fdr10": row.get("selected_fdr10", ""),
                    "selected_fdr5": row.get("selected_fdr5", ""),
                    "ivw_rf_to_pp_beta": row.get("ivw_rf_to_pp_beta", ""),
                    "ivw_rf_to_pp_p": row.get("ivw_rf_to_pp_p", ""),
                    "ivw_pp_to_outcome_beta": row.get("ivw_pp_to_outcome_beta", ""),
                    "ivw_pp_to_outcome_p": row.get("ivw_pp_to_outcome_p", ""),
                    "naive_indirect": "" if naive_indirect is None else f"{naive_indirect:.6f}",
                    "naive_support": "YES" if naive_support else "NO",
                    "naive_sign": naive_sign,
                    "bmediator_sign": bmediator_sign,
                    "sign_agreement": sign_agreement,
                    "threshold_agreement": threshold_agreement,
                    "agreement_flag": agreement_category(
                        naive_support, bmediator_support, sign_agreement
                    ),
                    "mediated_effect": row.get("mediated_effect", ""),
                    "se_mediated": row.get("se_mediated", ""),
                    "evidence_tier": row.get("evidence_tier", ""),
                    "nA": row.get("nA", ""),
                    "nB": row.get("nB", ""),
                    "nC": row.get("nC", ""),
                    "n_rf_to_pp_obs": row.get("n_rf_to_pp_obs", ""),
                }
            )
    return rows_out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--outcome", default="CHD_FinnGen")
    ap.add_argument("--pm1-threshold", type=float, default=0.5)
    ap.add_argument("--naive-p-threshold", type=float, default=0.05)
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows_out: list[dict[str, str]] = []
    for med_path in completed_mediation_files(results_dir):
        rows_out.extend(
            summarize_file(
                med_path,
                args.outcome,
                pm1_threshold=args.pm1_threshold,
                naive_p_threshold=args.naive_p_threshold,
            )
        )

    rows_out.sort(
        key=lambda row: (
            row["exposure"],
            1 if row["agreement_flag"] == "concordant_supported" else 0,
            float(row["P_M1"]),
        ),
        reverse=True,
    )

    out_path = output_dir / "naive_bmediator_agreement.tsv"
    fieldnames = list(rows_out[0].keys()) if rows_out else [
        "exposure",
        "outcome",
        "gene",
        "protein",
        "P_M1",
        "pm1_support",
        "selected_fdr10",
        "selected_fdr5",
        "ivw_rf_to_pp_beta",
        "ivw_rf_to_pp_p",
        "ivw_pp_to_outcome_beta",
        "ivw_pp_to_outcome_p",
        "naive_indirect",
        "naive_support",
        "naive_sign",
        "bmediator_sign",
        "sign_agreement",
        "threshold_agreement",
        "agreement_flag",
        "mediated_effect",
        "se_mediated",
        "evidence_tier",
        "nA",
        "nB",
        "nC",
        "n_rf_to_pp_obs",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    counts: dict[str, int] = {}
    for row in rows_out:
        counts[row["agreement_flag"]] = counts.get(row["agreement_flag"], 0) + 1

    print(out_path)
    print("__FLAG_COUNTS__")
    for key in sorted(counts):
        print(f"{key}\t{counts[key]}")
    print("__TOP_CONCORDANT__")
    top = [
        row for row in rows_out
        if row["agreement_flag"] == "concordant_supported"
    ][:20]
    for row in top:
        print(
            "\t".join(
                row[k]
                for k in [
                    "exposure",
                    "gene",
                    "protein",
                    "P_M1",
                    "naive_indirect",
                    "sign_agreement",
                    "threshold_agreement",
                    "agreement_flag",
                    "selected_fdr10",
                    "selected_fdr5",
                ]
            )
        )


if __name__ == "__main__":
    main()
