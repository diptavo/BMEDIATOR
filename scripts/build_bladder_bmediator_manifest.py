#!/usr/bin/env python3
import csv
from pathlib import Path


ROOT = Path("/data/Dutta_lab/BMEDIATOR")
INPUT_ROOT = ROOT / "generated" / "bmediator_inputs" / "bladder_bmi_smoking"
ANALYSIS_ROOT = ROOT / "analysis" / "bladder_bmi_smoking_ukbb_combined_b37"

PROTEIN_SUMSTAT = ROOT / "generated" / "UKBB_COMBINED_resource" / "merged" / "ukbb_combined_oid_cis_tss_pm6_r2_0.1.tsv"
PROTEIN_INFO = ROOT / "generated" / "UKBB_COMBINED_resource" / "merged" / "ukbb_combined_oid_cis_tss_pm6_r2_0.1.protein_info.tsv"

EXPOSURES = [
    ("BMI_IEU_B_40", "bmi_ieu_b_40", INPUT_ROOT / "exposures" / "bmi_ieu_b_40_b37_bmediator.tsv"),
    ("Smoking_IEU_B_4877", "smoking_ieu_b_4877", INPUT_ROOT / "exposures" / "smoking_ieu_b_4877_b37_bmediator.tsv"),
    ("Smoking_IEU_B_142", "smoking_ieu_b_142", INPUT_ROOT / "exposures" / "smoking_ieu_b_142_b37_bmediator.tsv"),
]

OUTCOMES = [
    ("Bladder_EUR", "bladder_eur", INPUT_ROOT / "outcomes" / "bladder_eur_b37_bmediator.tsv"),
    ("Bladder_MultiAncestry", "bladder_multiancestry", INPUT_ROOT / "outcomes" / "bladder_multiancestry_b37_bmediator.tsv"),
]


def main():
    manifest = ANALYSIS_ROOT / "manifests" / "bladder_bmi_smoking_manifest.tsv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="") as handle:
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
            lineterminator="\n",
        )
        writer.writeheader()
        for exposure, exposure_slug, rf_sumstat in EXPOSURES:
            for outcome, outcome_slug, outcome_sumstat in OUTCOMES:
                out_prefix = ANALYSIS_ROOT / "results" / outcome_slug / exposure_slug / f"{exposure_slug}_to_{outcome_slug}"
                writer.writerow(
                    {
                        "domain": "bladder",
                        "exposure": exposure,
                        "outcome": outcome,
                        "rf_sumstat": str(rf_sumstat),
                        "outcome_sumstat": str(outcome_sumstat),
                        "protein_sumstat": str(PROTEIN_SUMSTAT),
                        "protein_info": str(PROTEIN_INFO),
                        "out_prefix": str(out_prefix),
                    }
                )
    summary = ANALYSIS_ROOT / "manifests" / "manifest_summary.tsv"
    with summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["manifest", "rows"])
        writer.writerow([str(manifest), len(EXPOSURES) * len(OUTCOMES)])
    print(manifest)


if __name__ == "__main__":
    main()
