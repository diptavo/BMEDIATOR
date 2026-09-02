#!/usr/bin/env python3
import csv
from pathlib import Path


ROOT = Path("/data/Dutta_lab/BMEDIATOR")
INPUT_ROOT = ROOT / "generated" / "bmediator_inputs"
ANALYSIS_ROOT = ROOT / "analysis" / "cardiomet_vascular_kidney_ukbb_combined_tss_pm6_r2_0.1_hg38"

PROTEIN_SUMSTAT = ROOT / "generated" / "UKBB_COMBINED_resource" / "merged" / "ukbb_combined_oid_cis_tss_pm6_r2_0.1.tsv"
PROTEIN_INFO = ROOT / "generated" / "UKBB_COMBINED_resource" / "merged" / "ukbb_combined_oid_cis_tss_pm6_r2_0.1.protein_info.tsv"

EXPOSURES = [
    ("BMI_GIANT", INPUT_ROOT / "exposures" / "bmi_giant_locke_eur_hg38_bmediator.tsv", "bmi_giant"),
    ("BMI_MVP", INPUT_ROOT / "exposures" / "bmi_mvp_gcst90475156_hg38_bmediator.tsv", "bmi_mvp"),
    ("WHR_GIANT", INPUT_ROOT / "exposures" / "whr_giant_2015_eur_hg38_bmediator.tsv", "whr_giant"),
    ("WHRadjBMI_GIANT", INPUT_ROOT / "exposures" / "whradjbmi_giant_2015_eur_hg38_bmediator.tsv", "whradjbmi_giant"),
    ("SBP_KEATON", INPUT_ROOT / "exposures" / "sbp_keaton_gcst90310294_hg38_bmediator.tsv", "sbp_keaton"),
    ("DBP_KEATON", INPUT_ROOT / "exposures" / "dbp_keaton_gcst90310295_hg38_bmediator.tsv", "dbp_keaton"),
    ("PP_KEATON", INPUT_ROOT / "exposures" / "pp_keaton_gcst90310296_hg38_bmediator.tsv", "pp_keaton"),
    ("SBP_MVP", INPUT_ROOT / "exposures" / "sbp_mvp_gcst90476403_hg38_bmediator.tsv", "sbp_mvp"),
    ("DBP_MVP", INPUT_ROOT / "exposures" / "dbp_mvp_gcst90475252_hg38_bmediator.tsv", "dbp_mvp"),
    ("HDL_GLGC", INPUT_ROOT / "exposures" / "hdl_glgc_no_ukb_eur_hg38_bmediator.tsv", "hdl_glgc"),
    ("LDL_GLGC", INPUT_ROOT / "exposures" / "ldl_glgc_no_ukb_eur_hg38_bmediator.tsv", "ldl_glgc"),
    ("TG_GLGC", INPUT_ROOT / "exposures" / "tg_glgc_no_ukb_eur_hg38_bmediator.tsv", "tg_glgc"),
    ("TC_GLGC", INPUT_ROOT / "exposures" / "tc_glgc_no_ukb_eur_hg38_bmediator.tsv", "tc_glgc"),
]

VASCULAR_OUTCOMES = [
    ("CHD_FinnGen", INPUT_ROOT / "outcomes" / "vascular" / "chd_finngen_r12_hg38_bmediator.tsv", "chd_finngen"),
    ("MI_FinnGen", INPUT_ROOT / "outcomes" / "vascular" / "mi_finngen_r12_strict_hg38_bmediator.tsv", "mi_finngen"),
    ("Stroke_FinnGen", INPUT_ROOT / "outcomes" / "vascular" / "stroke_finngen_r12_hg38_bmediator.tsv", "stroke_finngen"),
    ("HeartFail_FinnGen", INPUT_ROOT / "outcomes" / "vascular" / "heartfail_finngen_r12_hg38_bmediator.tsv", "heartfail_finngen"),
    ("PADProxy_FinnGen", INPUT_ROOT / "outcomes" / "vascular" / "padproxy_finngen_r12_othper_hg38_bmediator.tsv", "padproxy_finngen"),
]

KIDNEY_OUTCOMES = [
    ("CKD_CKDGen", INPUT_ROOT / "outcomes" / "kidney" / "ckd_ckdgen_2019_hg38_bmediator.tsv", "ckd_ckdgen"),
    ("Microalbuminuria_CKDGen", INPUT_ROOT / "outcomes" / "kidney" / "microalbuminuria_ckdgen_2019_hg38_bmediator.tsv", "microalbuminuria_ckdgen"),
    ("UACR_CKDGen", INPUT_ROOT / "outcomes" / "kidney" / "uacr_ckdgen_2019_hg38_bmediator.tsv", "uacr_ckdgen"),
    ("RenFail_FinnGen", INPUT_ROOT / "outcomes" / "kidney" / "renfail_finngen_r12_hg38_bmediator.tsv", "renfail_finngen"),
    ("KidneyStones_FinnGen", INPUT_ROOT / "outcomes" / "kidney" / "kidneystones_finngen_r12_hg38_bmediator.tsv", "kidneystones_finngen"),
]

PILOT_ROWS = [
    ("vascular", "LDL_GLGC", "HeartFail_FinnGen"),
    ("vascular", "TC_GLGC", "Stroke_FinnGen"),
    ("kidney", "SBP_MVP", "CKD_CKDGen"),
    ("kidney", "DBP_KEATON", "UACR_CKDGen"),
]


def build_lookup(items):
    return {name: (str(path), slug) for name, path, slug in items}


def make_rows(domain, outcomes, exposure_lookup):
    rows = []
    for exp_name, rf_sumstat, exp_slug in EXPOSURES:
        for out_name, out_path, out_slug in outcomes:
            out_prefix = ANALYSIS_ROOT / "results" / domain / out_slug / exp_slug / f"{exp_slug}_to_{out_slug}"
            rows.append(
                {
                    "domain": domain,
                    "exposure": exp_name,
                    "outcome": out_name,
                    "rf_sumstat": str(rf_sumstat),
                    "outcome_sumstat": str(out_path),
                    "protein_sumstat": str(PROTEIN_SUMSTAT),
                    "protein_info": str(PROTEIN_INFO),
                    "out_prefix": str(out_prefix),
                }
            )
    return rows


def write_manifest(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
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


def main():
    exposure_lookup = build_lookup(EXPOSURES)
    vascular_lookup = build_lookup(VASCULAR_OUTCOMES)
    kidney_lookup = build_lookup(KIDNEY_OUTCOMES)

    vascular_rows = make_rows("vascular", VASCULAR_OUTCOMES, exposure_lookup)
    kidney_rows = make_rows("kidney", KIDNEY_OUTCOMES, exposure_lookup)

    pilot_rows = []
    for domain, exposure_name, outcome_name in PILOT_ROWS:
        rf_sumstat, exp_slug = exposure_lookup[exposure_name]
        if domain == "vascular":
            outcome_sumstat, out_slug = vascular_lookup[outcome_name]
        else:
            outcome_sumstat, out_slug = kidney_lookup[outcome_name]
        out_prefix = ANALYSIS_ROOT / "pilot" / domain / out_slug / exp_slug / f"{exp_slug}_to_{out_slug}_pilot"
        pilot_rows.append(
            {
                "domain": domain,
                "exposure": exposure_name,
                "outcome": outcome_name,
                "rf_sumstat": rf_sumstat,
                "outcome_sumstat": outcome_sumstat,
                "protein_sumstat": str(PROTEIN_SUMSTAT),
                "protein_info": str(PROTEIN_INFO),
                "out_prefix": str(out_prefix),
            }
        )

    write_manifest(ANALYSIS_ROOT / "manifests" / "pilot_manifest.tsv", pilot_rows)
    write_manifest(ANALYSIS_ROOT / "manifests" / "vascular_manifest.tsv", vascular_rows)
    write_manifest(ANALYSIS_ROOT / "manifests" / "kidney_manifest.tsv", kidney_rows)

    summary_path = ANALYSIS_ROOT / "manifests" / "manifest_summary.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["manifest", "rows"])
        writer.writerow(["pilot", len(pilot_rows)])
        writer.writerow(["vascular", len(vascular_rows)])
        writer.writerow(["kidney", len(kidney_rows)])
        writer.writerow(["total_full", len(vascular_rows) + len(kidney_rows)])


if __name__ == "__main__":
    main()
