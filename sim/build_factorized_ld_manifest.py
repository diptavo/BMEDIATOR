#!/usr/bin/env python3
"""Build a SLURM manifest for the factorized LD stress test."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from run_factorized_ld_stress import LD_CELLS, SCENARIOS


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("cell", "scenario", "replicate"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for cell in LD_CELLS:
            for scenario in SCENARIOS:
                for replicate in range(1, args.replicates + 1):
                    writer.writerow({"cell": cell, "scenario": scenario, "replicate": replicate})
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
