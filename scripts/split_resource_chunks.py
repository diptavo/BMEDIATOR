#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Split a protein resource TSV into chunk TSVs.")
    parser.add_argument("--resource", required=True, help="Input resource TSV")
    parser.add_argument("--chunk-size", required=True, type=int, help="Proteins per chunk")
    parser.add_argument("--out-dir", required=True, help="Chunk output directory")
    parser.add_argument("--prefix", required=True, help="Chunk file prefix")
    args = parser.parse_args()

    resource_path = Path(args.resource)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with resource_path.open() as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    if not rows or not fieldnames:
        raise SystemExit("resource file has no rows")

    manifest_path = out_dir / f"{args.prefix}_chunk_manifest.tsv"
    with manifest_path.open("w", newline="") as mf:
        manifest_writer = csv.DictWriter(
            mf,
            fieldnames=["chunk_id", "chunk_label", "resource_subset", "n_rows"],
            delimiter="\t",
        )
        manifest_writer.writeheader()
        for idx in range(0, len(rows), args.chunk_size):
            chunk_rows = rows[idx : idx + args.chunk_size]
            chunk_id = idx // args.chunk_size + 1
            chunk_label = f"{args.prefix}_chunk_{chunk_id:03d}"
            chunk_path = out_dir / f"{chunk_label}.tsv"
            with chunk_path.open("w", newline="") as cf:
                writer = csv.DictWriter(cf, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerows(chunk_rows)
            manifest_writer.writerow(
                {
                    "chunk_id": chunk_id,
                    "chunk_label": chunk_label,
                    "resource_subset": str(chunk_path),
                    "n_rows": len(chunk_rows),
                }
            )


if __name__ == "__main__":
    main()
