#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def normalize_chr(value: str) -> str:
    token = value.strip().upper().removeprefix("CHR")
    if token == "X":
        return "23"
    if token == "Y":
        return "24"
    return token


def main() -> None:
    parser = argparse.ArgumentParser(description="Build canonical ARIC SeqId resource table.")
    parser.add_argument("--seqid-table", required=True, help="ARIC seqid.txt metadata file")
    parser.add_argument("--sumstat-dir", required=True, help="Directory with SeqId_*.PHENO1.glm.linear files")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    seqid_table = Path(args.seqid_table)
    sumstat_dir = Path(args.sumstat_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with seqid_table.open() as fin, out_path.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = [
            "protein",
            "gene",
            "uniprot",
            "sumstat_file",
            "chr",
            "tss",
            "resource_source",
        ]
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            protein = row["seqid_in_sample"]
            writer.writerow(
                {
                    "protein": protein,
                    "gene": row["entrezgenesymbol"],
                    "uniprot": row["uniprot_id"],
                    "sumstat_file": str(sumstat_dir / f"{protein}.PHENO1.glm.linear"),
                    "chr": normalize_chr(row["chromosome_name"]),
                    "tss": row["transcription_start_site"],
                    "resource_source": "ARIC",
                }
            )


if __name__ == "__main__":
    main()
