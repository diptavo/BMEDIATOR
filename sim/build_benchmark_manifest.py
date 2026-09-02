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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.read_text())
    benchmarks = ["classification", "calibration"] if args.benchmark == "all" else [args.benchmark]
    rows = []
    for benchmark in benchmarks:
        bench_cfg = config[benchmark]
        for cell in bench_cfg["cells"]:
            for rep in range(1, int(bench_cfg["replicates"]) + 1):
                rows.append(
                    {
                        "benchmark": benchmark,
                        "cell": cell["name"],
                        "replicate": rep,
                    }
                )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["benchmark", "cell", "replicate"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
