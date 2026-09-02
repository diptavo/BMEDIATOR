#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def slug(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lifestyle-manifest", required=True)
    ap.add_argument("--exposure-dir", required=True)
    ap.add_argument("--outcome-eur", required=True)
    ap.add_argument("--outcome-multi", required=True)
    ap.add_argument("--protein-sumstat", required=True)
    ap.add_argument("--protein-info", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out-manifest", required=True)
    args = ap.parse_args()

    exposure_dir = Path(args.exposure_dir)
    results_dir = Path(args.results_dir)
    out_manifest = Path(args.out_manifest)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)

    outcomes = [
        ("bladder_eur", Path(args.outcome_eur)),
        ("bladder_multiancestry", Path(args.outcome_multi)),
    ]

    with Path(args.lifestyle_manifest).open() as handle, out_manifest.open("w", newline="") as out:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = ["domain", "exposure", "outcome", "rf_sumstat", "outcome_sumstat", "protein_sumstat", "protein_info", "out_prefix"]
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()

        for row in reader:
            exposure_id = row["id"]
            category = slug(row["category"])
            trait = slug(row["trait"])
            exposure_name = f"{exposure_id}_{trait}"
            rf_sumstat = exposure_dir / f"{exposure_name}_b37_bmediator.tsv"

            for outcome_name, outcome_sumstat in outcomes:
                out_prefix = results_dir / category / outcome_name / exposure_name / f"{exposure_name}_to_{outcome_name}"
                writer.writerow(
                    {
                        "domain": category,
                        "exposure": exposure_name,
                        "outcome": outcome_name,
                        "rf_sumstat": str(rf_sumstat),
                        "outcome_sumstat": str(outcome_sumstat),
                        "protein_sumstat": args.protein_sumstat,
                        "protein_info": args.protein_info,
                        "out_prefix": str(out_prefix),
                    }
                )


if __name__ == "__main__":
    main()
