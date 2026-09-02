#!/usr/bin/env python3
import csv
from pathlib import Path


ANALYSIS_ROOT = Path("/data/Dutta_lab/BMEDIATOR/analysis/cardiomet_cross_pilot_fix")
MANIFEST = ANALYSIS_ROOT / "manifests" / "cross_pilot_manifest.tsv"
SUMMARY_DIR = ANALYSIS_ROOT / "summary"


def read_manifest():
    with MANIFEST.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_run(row):
    mediation_path = Path(row["out_prefix"] + ".mediation")
    hyp_path = Path(row["out_prefix"] + ".hyp")

    summary = {
        "domain": row["domain"],
        "exposure": row["exposure"],
        "outcome": row["outcome"],
        "out_prefix": row["out_prefix"],
        "mediation_exists": "NO",
        "hyp_exists": "NO",
        "proteins_total": 0,
        "proteins_analyzable": 0,
        "proteins_nB_positive": 0,
        "proteins_nC_positive": 0,
        "proteins_rf_obs_positive": 0,
        "proteins_converged_yes": 0,
        "pm1_gt_0.001": 0,
        "pm1_gt_0.01": 0,
        "pm1_gt_0.05": 0,
        "max_pm1": 0.0,
        "top_protein": "",
        "top_gene": "",
        "top_pm1": 0.0,
        "p1": "",
        "p2": "",
        "p3": "",
    }

    if mediation_path.exists():
        summary["mediation_exists"] = "YES"
        with mediation_path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for rec in reader:
                summary["proteins_total"] += 1
                nB = int(rec["nB"])
                nC = int(rec["nC"])
                n_rf = int(rec["n_rf_to_pp_obs"])
                pm1 = float(rec["P_M1"])
                if nB > 0:
                    summary["proteins_nB_positive"] += 1
                if nC > 0:
                    summary["proteins_nC_positive"] += 1
                if n_rf > 0:
                    summary["proteins_rf_obs_positive"] += 1
                if rec["rf_to_pp_identifiable"] == "YES":
                    summary["proteins_analyzable"] += 1
                if rec["converged"] == "YES":
                    summary["proteins_converged_yes"] += 1
                if pm1 > 0.001:
                    summary["pm1_gt_0.001"] += 1
                if pm1 > 0.01:
                    summary["pm1_gt_0.01"] += 1
                if pm1 > 0.05:
                    summary["pm1_gt_0.05"] += 1
                if pm1 > summary["max_pm1"]:
                    summary["max_pm1"] = pm1
                    summary["top_pm1"] = pm1
                    summary["top_protein"] = rec["Protein"]
                    summary["top_gene"] = rec["Gene"]

    if hyp_path.exists():
        summary["hyp_exists"] = "YES"
        with hyp_path.open() as handle:
            next(handle)
            values = next(handle).rstrip("\n").split("\t")
        header = ["p0", "p1", "p2", "p3", "pi1", "pi2_cis", "pi3"]
        hyp = dict(zip(header, values))
        summary["p1"] = hyp.get("p1", "")
        summary["p2"] = hyp.get("p2", "")
        summary["p3"] = hyp.get("p3", "")

    return summary


def main():
    rows = read_manifest()
    summaries = [summarize_run(row) for row in rows]
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    summary_path = SUMMARY_DIR / "cross_pilot_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    headline_path = SUMMARY_DIR / "cross_pilot_headline.txt"
    total = len(summaries)
    done = sum(1 for row in summaries if row["mediation_exists"] == "YES")
    best = max(summaries, key=lambda row: row["max_pm1"])
    with headline_path.open("w") as handle:
        handle.write(f"pilots_total\t{total}\n")
        handle.write(f"pilots_completed\t{done}\n")
        handle.write(
            "best_signal\t{exposure}->{outcome}\t{gene}\t{pm1:.6f}\n".format(
                exposure=best["exposure"],
                outcome=best["outcome"],
                gene=best["top_gene"],
                pm1=best["max_pm1"],
            )
        )


if __name__ == "__main__":
    main()
