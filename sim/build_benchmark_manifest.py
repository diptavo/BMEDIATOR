from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a SLURM manifest for BMEDIATOR simulation benchmarks.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--benchmark", choices=["classification", "calibration", "all"], default="all")
    parser.add_argument(
        "--replicates",
        type=int,
        help="Override the configured replicate count when building a validation manifest.",
    )
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.read_text())
    benchmarks = ["classification", "calibration"] if args.benchmark == "all" else [args.benchmark]
    rows = []
    for benchmark in benchmarks:
        bench_cfg = config[benchmark]
        for cell in bench_cfg["cells"]:
            replicates = args.replicates or int(
                cell.get("replicates", bench_cfg["replicates"])
            )
            for rep in range(1, replicates + 1):
                rows.append(
                    {
                        "benchmark": benchmark,
                        "cell": cell["name"],
                        "replicate": rep,
                    }
                )
    if args.shard_count < 1 or not 1 <= args.shard_index <= args.shard_count:
        raise ValueError("shard index must be between 1 and shard count")
    rows = rows[args.shard_index - 1 :: args.shard_count]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["benchmark", "cell", "replicate"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
