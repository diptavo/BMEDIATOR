#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


ALLELES = {"A", "C", "G", "T"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def slug(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def clean_float(value: str) -> float | None:
    try:
        out = float(value)
    except Exception:
        return None
    if not math.isfinite(out):
        return None
    return out


def load_bim(path: Path) -> dict[str, tuple[str, str, str, str]]:
    mapping: dict[str, tuple[str, str, str, str]] = {}
    with path.open() as handle:
        for line in handle:
            chrom, snp, _cm, bp, a1, a2 = line.rstrip("\n").split()[:6]
            mapping[snp] = (chrom, bp, a1.upper(), a2.upper())
    return mapping


def alleles_match(a1: str, a2: str, ref1: str, ref2: str) -> bool:
    query = {a1, a2}
    ref = {ref1, ref2}
    if query == ref:
        return True
    comp = {a1.translate(COMPLEMENT), a2.translate(COMPLEMENT)}
    return comp == ref


def convert_one(input_ma: Path, output_tsv: Path, bim: dict[str, tuple[str, str, str, str]]) -> dict[str, int]:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    stats = {
        "input_rows": 0,
        "written_rows": 0,
        "missing_bim": 0,
        "bad_allele": 0,
        "allele_mismatch": 0,
        "bad_numeric": 0,
    }

    with input_ma.open() as inp, output_tsv.open("w", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "CHR", "BP"])

        for row in reader:
            stats["input_rows"] += 1
            snp = row.get("SNP", "").strip()
            if snp not in bim:
                stats["missing_bim"] += 1
                continue

            a1 = row.get("A1", "").upper()
            a2 = row.get("A2", "").upper()
            if len(a1) != 1 or len(a2) != 1 or a1 not in ALLELES or a2 not in ALLELES:
                stats["bad_allele"] += 1
                continue

            chrom, bp, ref1, ref2 = bim[snp]
            if not alleles_match(a1, a2, ref1, ref2):
                stats["allele_mismatch"] += 1
                continue

            freq = clean_float(row.get("freq", ""))
            beta = clean_float(row.get("b", ""))
            se = clean_float(row.get("se", ""))
            pval = clean_float(row.get("p", ""))
            if freq is None or beta is None or se is None or pval is None:
                stats["bad_numeric"] += 1
                continue
            if not (0.0 < freq < 1.0 and se > 0.0 and 0.0 < pval <= 1.0):
                stats["bad_numeric"] += 1
                continue

            writer.writerow([snp, a1, a2, f"{freq:.12g}", f"{beta:.12g}", f"{se:.12g}", f"{pval:.12g}", chrom, bp])
            stats["written_rows"] += 1

    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--gsmr-dir", required=True)
    ap.add_argument("--bim", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    manifest = Path(args.manifest)
    gsmr_dir = Path(args.gsmr_dir)
    out_dir = Path(args.out_dir)
    summary_path = Path(args.summary)

    bim = load_bim(Path(args.bim))
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    with manifest.open() as handle, summary_path.open("w", newline="") as out:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = ["id", "category", "trait", "input_ma", "output_tsv", "input_rows", "written_rows", "missing_bim", "bad_allele", "allele_mismatch", "bad_numeric"]
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            exposure_id = row["id"]
            trait = slug(row["trait"])
            input_ma = gsmr_dir / f"{exposure_id}_{trait}.exposure.gsmr.dedup.ma"
            output_tsv = out_dir / f"{exposure_id}_{trait}_b37_bmediator.tsv"
            stats = convert_one(input_ma, output_tsv, bim)
            writer.writerow(
                {
                    "id": exposure_id,
                    "category": row["category"],
                    "trait": trait,
                    "input_ma": str(input_ma),
                    "output_tsv": str(output_tsv),
                    **stats,
                }
            )


if __name__ == "__main__":
    main()
