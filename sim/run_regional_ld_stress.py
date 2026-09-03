#!/usr/bin/env python3
"""Genotype-based stress test for BMEDIATOR regional LD resolution.

The shared and distinct scenarios test the implemented H4-versus-H3 gate.
The same_variant_pleiotropy scenario is intentionally observationally
equivalent to mediation and measures the documented exclusion-assumption limit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np


SCENARIOS = ("shared_mediation", "distinct_ld", "same_variant_pleiotropy")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path(__file__).resolve().parents[1] / "bmediator")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, default=None, help="Temporary input directory; defaults to OUTDIR/work.")
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--ld", type=float, nargs="+", default=(0.3, 0.6, 0.9))
    parser.add_argument("--variants", type=int, default=20)
    parser.add_argument("--reference-samples", type=int, default=512)
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--regional-prior-pp", type=float, default=1e-4)
    parser.add_argument("--regional-prior-outcome", type=float, default=1e-4)
    parser.add_argument("--regional-prior-shared", type=float, default=1e-8)
    parser.add_argument("--keep-inputs", action="store_true")
    return parser.parse_args()


def p_value(beta: float, se: float) -> float:
    return math.erfc(abs(beta / se) / math.sqrt(2.0))


def ar1_genotypes(rng: np.random.Generator, n: int, m: int, rho: float) -> np.ndarray:
    indices = np.arange(m)
    covariance = rho ** np.abs(indices[:, None] - indices[None, :])
    haplotypes = rng.multivariate_normal(np.zeros(m), covariance, size=2 * n) > 0.0
    return haplotypes.reshape(n, 2, m).sum(axis=1).astype(np.int8)


def independent_genotypes(rng: np.random.Generator, n: int, m: int) -> np.ndarray:
    return rng.binomial(2, 0.5, size=(n, m)).astype(np.int8)


def correlation_matrix(genotypes: np.ndarray) -> np.ndarray:
    return np.corrcoef(genotypes.astype(float), rowvar=False)


def write_plink(prefix: Path, genotypes: np.ndarray, variants: list[tuple[int, str, int]]) -> None:
    n_samples, n_variants = genotypes.shape
    with prefix.with_suffix(".fam").open("w", encoding="ascii") as handle:
        for idx in range(n_samples):
            handle.write(f"F{idx + 1} I{idx + 1} 0 0 0 -9\n")
    with prefix.with_suffix(".bim").open("w", encoding="ascii") as handle:
        for chrom, rsid, bp in variants:
            handle.write(f"{chrom} {rsid} 0 {bp} A G\n")

    code = {2: 0b00, 1: 0b10, 0: 0b11}
    with prefix.with_suffix(".bed").open("wb") as handle:
        handle.write(bytes((0x6C, 0x1B, 0x01)))
        for variant_idx in range(n_variants):
            dosage = genotypes[:, variant_idx]
            for start in range(0, n_samples, 4):
                byte = 0
                for offset, value in enumerate(dosage[start : start + 4]):
                    byte |= code[int(value)] << (2 * offset)
                handle.write(struct.pack("B", byte))


def write_sumstats(
    path: Path,
    variants: list[tuple[int, str, int]],
    frequencies: np.ndarray,
    effects: np.ndarray,
    standard_errors: np.ndarray,
) -> None:
    with path.open("w", encoding="ascii") as handle:
        handle.write("SNP A1 A2 FREQ BETA SE P CHR BP\n")
        for (chrom, rsid, bp), freq, beta, se in zip(variants, frequencies, effects, standard_errors):
            handle.write(
                f"{rsid} A G {freq:.8g} {beta:.8g} {se:.8g} "
                f"{p_value(float(beta), float(se)):.8g} {chrom} {bp}\n"
            )


def write_protein_gwas(
    path: Path,
    variants: list[tuple[int, str, int]],
    frequencies: np.ndarray,
    effects: np.ndarray,
    standard_errors: np.ndarray,
) -> None:
    header = "#CHROM\tPOS\tID\tREF\tALT\tA1\tA1_FREQ\tTEST\tBETA\tSE\tP\tERRCODE\n"
    with path.open("w", encoding="ascii") as handle:
        handle.write(header)
        for (chrom, rsid, bp), freq, beta, se in zip(variants, frequencies, effects, standard_errors):
            handle.write(
                f"{chrom}\t{bp}\t{rsid}\tG\tA\tA\t{freq:.8g}\tADD\t"
                f"{beta:.8g}\t{se:.8g}\t{p_value(float(beta), float(se)):.8g}\t.\n"
            )


def make_inputs(
    directory: Path,
    rng: np.random.Generator,
    scenario: str,
    rho: float,
    n_reference: int,
    n_cis: int,
) -> tuple[Path, float]:
    directory.mkdir(parents=True, exist_ok=True)
    cis_genotypes = ar1_genotypes(rng, n_reference, n_cis, rho)
    rf_genotypes = independent_genotypes(rng, n_reference, 2)
    genotypes = np.column_stack((cis_genotypes, rf_genotypes))
    cis_variants = [(1, f"cis{idx + 1}", 100_000 + idx * 1_000) for idx in range(n_cis)]
    variants = cis_variants + [(2, "rf1", 1_000_000), (2, "rf2", 3_000_000)]
    write_plink(directory / "ldref", genotypes, variants)

    regional_ld = correlation_matrix(cis_genotypes)
    protein_causal = n_cis // 3
    outcome_causal = protein_causal if scenario != "distinct_ld" else protein_causal + 1
    realized_ld = float(regional_ld[protein_causal, outcome_causal])
    protein_mean = 0.20 * regional_ld[:, protein_causal]
    outcome_mean = 0.14 * regional_ld[:, outcome_causal]

    protein_se = np.full(n_cis + 2, 0.02)
    outcome_se = np.full(n_cis + 2, 0.02)
    rf_se = np.full(n_cis + 2, 0.02)
    protein_effect = np.concatenate((protein_mean, np.array([0.12, 0.10])))
    outcome_effect = np.concatenate((outcome_mean, np.array([0.05, 0.04])))
    rf_effect = np.concatenate((np.zeros(n_cis), np.array([0.12, 0.10])))

    protein_effect += rng.normal(0.0, protein_se)
    outcome_effect += rng.normal(0.0, outcome_se)
    rf_effect += rng.normal(0.0, rf_se)
    # Keep the designed RF instruments strong in every replicate.
    rf_effect[-2:] = np.array([0.12, 0.10]) + rng.normal(0.0, 0.005, size=2)

    frequencies = genotypes.mean(axis=0) / 2.0
    write_sumstats(directory / "rf.txt", variants, frequencies, rf_effect, rf_se)
    write_sumstats(directory / "outcome.txt", variants, frequencies, outcome_effect, outcome_se)
    write_protein_gwas(directory / "protein.glm.linear", variants, frequencies, protein_effect, protein_se)
    (directory / "manifest.txt").write_text(f"P1 {directory / 'protein.glm.linear'}\n", encoding="ascii")
    (directory / "protein_info.txt").write_text(
        "PROTEIN GENE CHR START END\nP1 G1 1 106000 114000\n",
        encoding="ascii",
    )
    return directory / "ldref", realized_ld


def read_result(path: Path) -> dict[str, str]:
    with path.open(encoding="ascii") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"expected one result row in {path}, found {len(rows)}")
    return rows[0]


def run_one(
    binary: Path,
    directory: Path,
    ld_prefix: Path,
    regional_priors: tuple[float, float, float],
) -> dict[str, str]:
    output = directory / "bmediator"
    cmd = [
        str(binary),
        "--rf-sumstat", str(directory / "rf.txt"),
        "--protein-gwas-list", str(directory / "manifest.txt"),
        "--cancer-sumstat", str(directory / "outcome.txt"),
        "--protein-info", str(directory / "protein_info.txt"),
        "--bfile", str(ld_prefix),
        "--out", str(output),
        "--p-thresh-rf", "5e-6",
        "--p-thresh-cis", "5e-6",
        "--cis-window", "50",
        "--min-instruments", "1",
        "--heidi-off",
        "--no-steiger",
        "--max-cavi-iter", "100",
        "--regional-prior-pp", str(regional_priors[0]),
        "--regional-prior-outcome", str(regional_priors[1]),
        "--regional-prior-shared", str(regional_priors[2]),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return read_result(output.with_suffix(".mediation"))


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, ValueError):
        return math.nan


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    grouped: dict[tuple[str, float], list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault((str(row["scenario"]), float(row["target_ld"])), []).append(row)
    fields = [
        "scenario", "target_ld", "n", "mean_realized_ld", "mean_P_M1",
        "mean_P_mediator_ld_resolved", "ld_resolved_rate", "distinct_supported_rate",
    ]
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (scenario, target_ld), group in sorted(grouped.items()):
            n = len(group)
            writer.writerow(
                {
                    "scenario": scenario,
                    "target_ld": target_ld,
                    "n": n,
                    "mean_realized_ld": np.mean([float(row["realized_ld"]) for row in group]),
                    "mean_P_M1": np.mean([float(row["P_M1"]) for row in group]),
                    "mean_P_mediator_ld_resolved": np.mean(
                        [float(row["P_mediator_ld_resolved"]) for row in group]
                    ),
                    "ld_resolved_rate": np.mean(
                        [float(row["P_mediator_ld_resolved"]) > 0.5 for row in group]
                    ),
                    "distinct_supported_rate": np.mean(
                        [row["mediation_identifiability"] == "LD_DISTINCT_SUPPORTED" for row in group]
                    ),
                }
            )


def main() -> None:
    args = parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        raise SystemExit(f"binary not found: {binary}")
    if args.replicates < 1 or args.variants < 4 or args.reference_samples < 20:
        raise SystemExit("replicates >= 1, variants >= 4, and reference-samples >= 20 are required")
    if any(not (0.0 < value < 1.0) for value in args.ld):
        raise SystemExit("all --ld values must be between 0 and 1")

    args.outdir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "replicates": args.replicates,
        "target_ld": args.ld,
        "scenarios": SCENARIOS,
        "variants": args.variants,
        "reference_samples": args.reference_samples,
        "seed": args.seed,
        "workdir": str(args.workdir) if args.workdir is not None else str(args.outdir / "work"),
        "regional_prior_pp": args.regional_prior_pp,
        "regional_prior_outcome": args.regional_prior_outcome,
        "regional_prior_shared": args.regional_prior_shared,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
    }
    (args.outdir / "regional_ld_stress_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="ascii"
    )
    work = args.workdir if args.workdir is not None else args.outdir / "work"
    rows: list[dict[str, object]] = []
    total_runs = len(args.ld) * len(SCENARIOS) * args.replicates
    for ld_index, rho in enumerate(args.ld):
        for scenario_index, scenario in enumerate(SCENARIOS):
            for replicate in range(1, args.replicates + 1):
                seed = args.seed + ld_index * 1_000_000 + scenario_index * 100_000 + replicate
                rng = np.random.default_rng(seed)
                directory = work / f"ld_{rho:g}" / scenario / f"rep_{replicate:04d}"
                ld_prefix, realized_ld = make_inputs(
                    directory, rng, scenario, rho, args.reference_samples, args.variants
                )
                result = run_one(
                    binary,
                    directory,
                    ld_prefix,
                    (
                        args.regional_prior_pp,
                        args.regional_prior_outcome,
                        args.regional_prior_shared,
                    ),
                )
                rows.append(
                    {
                        "scenario": scenario,
                        "target_ld": rho,
                        "replicate": replicate,
                        "seed": seed,
                        "realized_ld": realized_ld,
                        "P_M1": as_float(result, "P_M1"),
                        "P_mediator_ld_resolved": as_float(result, "P_mediator_ld_resolved"),
                        "regional_PP_shared": as_float(result, "regional_PP_shared"),
                        "regional_PP_distinct": as_float(result, "regional_PP_distinct"),
                        "regional_shared_given_both": as_float(result, "regional_shared_given_both"),
                        "mediation_identifiability": result.get("mediation_identifiability", ""),
                    }
                )
                completed = len(rows)
                if completed % 100 == 0 or completed == total_runs:
                    print(f"Completed {completed} / {total_runs} runs", flush=True)

    raw_path = args.outdir / "regional_ld_stress.tsv"
    with raw_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary_path = args.outdir / "regional_ld_stress_summary.tsv"
    write_summary(summary_path, rows)
    if not args.keep_inputs:
        shutil.rmtree(work)
    print(f"Wrote {raw_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
