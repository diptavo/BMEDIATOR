#!/usr/bin/env python3
"""Aggregate array outputs from run_factorized_ld_task.py."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from run_factorized_ld_stress import LD_CELLS, SCENARIOS, write_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--replicates", type=int)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    if args.require_complete and args.replicates is None:
        raise SystemExit("--require-complete requires --replicates")
    paths = sorted(args.input.glob("*/*/rep_*/task_metrics.tsv"))
    if not paths:
        raise SystemExit(f"No task_metrics.tsv files under {args.input}")
    rows: list[dict[str, object]] = []
    for path in paths:
        with path.open(newline="", encoding="ascii") as handle:
            task_rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(task_rows) != 1:
            raise SystemExit(f"{path}: expected one row, found {len(task_rows)}")
        rows.extend(task_rows)
    keys = [(str(row["cell"]), str(row["scenario"]), int(row["replicate"])) for row in rows]
    if len(keys) != len(set(keys)):
        raise SystemExit("duplicate LD task identities found")
    if args.require_complete:
        expected = {
            (cell, scenario, replicate)
            for cell in LD_CELLS
            for scenario in SCENARIOS
            for replicate in range(1, args.replicates + 1)
        }
        missing = expected - set(keys)
        unexpected = set(keys) - expected
        if missing or unexpected:
            raise SystemExit(
                f"incomplete LD run: {len(missing)} missing and "
                f"{len(unexpected)} unexpected tasks"
            )
        binary_hashes = {str(row.get("binary_sha256", "")) for row in rows}
        if "" in binary_hashes or len(binary_hashes) != 1:
            raise SystemExit("LD tasks do not share one recorded binary hash")
    raw = args.input / "factorized_ld_stress.tsv"
    with raw.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    write_summary(args.input / "factorized_ld_stress_summary.tsv", rows)
    print(f"Summarized {len(rows)} LD stress tasks")


if __name__ == "__main__":
    main()
