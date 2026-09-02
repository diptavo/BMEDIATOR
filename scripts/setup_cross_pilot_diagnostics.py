#!/usr/bin/env python3
import csv
from pathlib import Path


ROOT = Path("/data/Dutta_lab/BMEDIATOR")
INPUT_ROOT = ROOT / "generated" / "bmediator_inputs"
ANALYSIS_ROOT = ROOT / "analysis" / "cardiomet_cross_pilot_fix"

PROTEIN_SUMSTAT = ROOT / "generated" / "UKBB_COMBINED_resource" / "merged" / "ukbb_combined_oid_cis_tss_pm6_r2_0.1.tsv"
PROTEIN_INFO = ROOT / "generated" / "UKBB_COMBINED_resource" / "merged" / "ukbb_combined_oid_cis_tss_pm6_r2_0.1.protein_info.tsv"

EXPOSURES = {
    "BMI_GIANT": (INPUT_ROOT / "exposures" / "bmi_giant_locke_eur_hg38_bmediator.tsv", "bmi_giant"),
    "BMI_MVP": (INPUT_ROOT / "exposures" / "bmi_mvp_gcst90475156_hg38_bmediator.tsv", "bmi_mvp"),
    "WHRadjBMI_GIANT": (INPUT_ROOT / "exposures" / "whradjbmi_giant_2015_eur_hg38_bmediator.tsv", "whradjbmi_giant"),
    "SBP_KEATON": (INPUT_ROOT / "exposures" / "sbp_keaton_gcst90310294_hg38_bmediator.tsv", "sbp_keaton"),
    "SBP_MVP": (INPUT_ROOT / "exposures" / "sbp_mvp_gcst90476403_hg38_bmediator.tsv", "sbp_mvp"),
    "DBP_KEATON": (INPUT_ROOT / "exposures" / "dbp_keaton_gcst90310295_hg38_bmediator.tsv", "dbp_keaton"),
    "DBP_MVP": (INPUT_ROOT / "exposures" / "dbp_mvp_gcst90475252_hg38_bmediator.tsv", "dbp_mvp"),
    "LDL_GLGC": (INPUT_ROOT / "exposures" / "ldl_glgc_no_ukb_eur_hg38_bmediator.tsv", "ldl_glgc"),
    "TG_GLGC": (INPUT_ROOT / "exposures" / "tg_glgc_no_ukb_eur_hg38_bmediator.tsv", "tg_glgc"),
    "TC_GLGC": (INPUT_ROOT / "exposures" / "tc_glgc_no_ukb_eur_hg38_bmediator.tsv", "tc_glgc"),
}

VASCULAR_OUTCOMES = {
    "CHD_FinnGen": (INPUT_ROOT / "outcomes" / "vascular" / "chd_finngen_r12_hg38_bmediator.tsv", "chd_finngen"),
    "MI_FinnGen": (INPUT_ROOT / "outcomes" / "vascular" / "mi_finngen_r12_strict_hg38_bmediator.tsv", "mi_finngen"),
    "Stroke_FinnGen": (INPUT_ROOT / "outcomes" / "vascular" / "stroke_finngen_r12_hg38_bmediator.tsv", "stroke_finngen"),
    "HeartFail_FinnGen": (INPUT_ROOT / "outcomes" / "vascular" / "heartfail_finngen_r12_hg38_bmediator.tsv", "heartfail_finngen"),
}

KIDNEY_OUTCOMES = {
    "CKD_CKDGen": (INPUT_ROOT / "outcomes" / "kidney" / "ckd_ckdgen_2019_hg38_bmediator.tsv", "ckd_ckdgen"),
    "RenFail_FinnGen": (INPUT_ROOT / "outcomes" / "kidney" / "renfail_finngen_r12_hg38_bmediator.tsv", "renfail_finngen"),
    "UACR_CKDGen": (INPUT_ROOT / "outcomes" / "kidney" / "uacr_ckdgen_2019_hg38_bmediator.tsv", "uacr_ckdgen"),
    "Microalbuminuria_CKDGen": (INPUT_ROOT / "outcomes" / "kidney" / "microalbuminuria_ckdgen_2019_hg38_bmediator.tsv", "microalbuminuria_ckdgen"),
    "KidneyStones_FinnGen": (INPUT_ROOT / "outcomes" / "kidney" / "kidneystones_finngen_r12_hg38_bmediator.tsv", "kidneystones_finngen"),
    "eGFRcrea_CKDGen": (INPUT_ROOT / "outcomes" / "kidney" / "egfrcrea_ckdgen_2017_hg38_bmediator.tsv", "egfrcrea_ckdgen"),
    "eGFRcys_CKDGen": (INPUT_ROOT / "outcomes" / "kidney" / "egfrcys_ckdgen_2017_hg38_bmediator.tsv", "egfrcys_ckdgen"),
}

PILOTS = [
    ("vascular", "BMI_GIANT", "CHD_FinnGen"),
    ("vascular", "LDL_GLGC", "CHD_FinnGen"),
    ("vascular", "WHRadjBMI_GIANT", "MI_FinnGen"),
    ("vascular", "TG_GLGC", "Stroke_FinnGen"),
    ("vascular", "SBP_KEATON", "HeartFail_FinnGen"),
    ("kidney", "BMI_GIANT", "CKD_CKDGen"),
    ("kidney", "SBP_MVP", "RenFail_FinnGen"),
    ("kidney", "DBP_KEATON", "UACR_CKDGen"),
    ("kidney", "TG_GLGC", "Microalbuminuria_CKDGen"),
    ("kidney", "TC_GLGC", "KidneyStones_FinnGen"),
    ("kidney", "BMI_MVP", "eGFRcrea_CKDGen"),
    ("kidney", "DBP_MVP", "eGFRcys_CKDGen"),
]


def outcome_lookup(domain, outcome_name):
    if domain == "vascular":
        return VASCULAR_OUTCOMES[outcome_name]
    return KIDNEY_OUTCOMES[outcome_name]


def build_rows():
    rows = []
    for domain, exposure_name, outcome_name in PILOTS:
        rf_sumstat, exp_slug = EXPOSURES[exposure_name]
        outcome_sumstat, out_slug = outcome_lookup(domain, outcome_name)
        out_prefix = ANALYSIS_ROOT / "results" / domain / out_slug / exp_slug / f"{exp_slug}_to_{out_slug}_pilot"
        rows.append(
            {
                "domain": domain,
                "exposure": exposure_name,
                "outcome": outcome_name,
                "rf_sumstat": str(rf_sumstat),
                "outcome_sumstat": str(outcome_sumstat),
                "protein_sumstat": str(PROTEIN_SUMSTAT),
                "protein_info": str(PROTEIN_INFO),
                "out_prefix": str(out_prefix),
            }
        )
    return rows


def main():
    rows = build_rows()
    manifest_dir = ANALYSIS_ROOT / "manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = manifest_dir / "cross_pilot_manifest.tsv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "domain",
                "exposure",
                "outcome",
                "rf_sumstat",
                "outcome_sumstat",
                "protein_sumstat",
                "protein_info",
                "out_prefix",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    summary_path = manifest_dir / "cross_pilot_manifest_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["n_pilots", len(rows)])
        writer.writerow(["n_vascular", sum(1 for row in rows if row["domain"] == "vascular")])
        writer.writerow(["n_kidney", sum(1 for row in rows if row["domain"] == "kidney")])


if __name__ == "__main__":
    main()
