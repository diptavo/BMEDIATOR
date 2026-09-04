#!/usr/bin/env python3
"""Genotype-based validation of factorized mediation and the regional LD gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np


SCENARIOS = ("shared_mediation", "distinct_ld", "same_variant_pleiotropy", "no_second_stage")
LD_CELLS = {
    "matched_moderate": (0.60, 0.60),
    "matched_high": (0.90, 0.90),
    "mismatch_high_to_moderate": (0.90, 0.60),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path(__file__).resolve().parents[1] / "bmediator")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--cells", nargs="+", choices=tuple(LD_CELLS), default=tuple(LD_CELLS))
    parser.add_argument("--reference-samples", type=int, default=1000)
    parser.add_argument("--analysis-samples", type=int, default=4000)
    parser.add_argument("--blocks", type=int, default=3)
    parser.add_argument("--variants-per-block", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--keep-inputs", action="store_true")
    return parser.parse_args()


def p_value(beta: float, se: float) -> float:
    return math.erfc(abs(beta / se) / math.sqrt(2.0))


def ar1_genotypes(rng: np.random.Generator, n: int, m: int, rho: float) -> np.ndarray:
    index = np.arange(m)
    covariance = rho ** np.abs(index[:, None] - index[None, :])
    haplotypes = rng.multivariate_normal(np.zeros(m), covariance, size=2 * n) > 0.0
    return haplotypes.reshape(n, 2, m).sum(axis=1).astype(np.int8)


def block_genotypes(
    rng: np.random.Generator, n: int, blocks: int, variants_per_block: int, rho: float
) -> np.ndarray:
    return np.column_stack([
        ar1_genotypes(rng, n, variants_per_block, rho) for _ in range(blocks)
    ])


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
        for (chrom, rsid, bp), freq, beta, se in zip(
            variants, frequencies, effects, standard_errors
        ):
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
    with path.open("w", encoding="ascii") as handle:
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tA1\tA1_FREQ\tTEST\tBETA\tSE\tP\tERRCODE\n")
        for (chrom, rsid, bp), freq, beta, se in zip(
            variants, frequencies, effects, standard_errors
        ):
            handle.write(
                f"{chrom}\t{bp}\t{rsid}\tG\tA\tA\t{freq:.8g}\tADD\t"
                f"{beta:.8g}\t{se:.8g}\t{p_value(float(beta), float(se)):.8g}\t.\n"
            )


def sample_correlated_error(
    rng: np.random.Generator, correlation: np.ndarray, standard_error: float
) -> np.ndarray:
    covariance = standard_error * standard_error * correlation
    return rng.multivariate_normal(np.zeros(correlation.shape[0]), covariance)


def make_inputs(
    directory: Path,
    rng: np.random.Generator,
    scenario: str,
    target_rho: float,
    reference_rho: float,
    n_analysis: int,
    n_reference: int,
    blocks: int,
    variants_per_block: int,
) -> tuple[Path, float]:
    directory.mkdir(parents=True, exist_ok=True)
    target_cis = block_genotypes(rng, n_analysis, blocks, variants_per_block, target_rho)
    reference_cis = block_genotypes(rng, n_reference, blocks, variants_per_block, reference_rho)
    target_rf = rng.binomial(2, 0.5, size=(n_analysis, 6)).astype(np.int8)
    reference_rf = rng.binomial(2, 0.5, size=(n_reference, 6)).astype(np.int8)
    target_genotypes = np.column_stack((target_cis, target_rf))
    reference_genotypes = np.column_stack((reference_cis, reference_rf))

    cis_variants = []
    for block in range(blocks):
        for offset in range(variants_per_block):
            index = block * variants_per_block + offset
            cis_variants.append((1, f"cis{index + 1}", 100_000 + block * 100_000 + offset * 1_000))
    rf_variants = [(chrom, f"rf{chrom - 1}", 1_000_000) for chrom in range(2, 8)]
    variants = cis_variants + rf_variants
    write_plink(directory / "ldref", reference_genotypes, variants)

    target_ld = np.corrcoef(target_genotypes.astype(float), rowvar=False)
    cis_ld = target_ld[: len(cis_variants), : len(cis_variants)]
    protein_causal = [block * variants_per_block + variants_per_block // 2 for block in range(blocks)]
    outcome_causal = [index if scenario != "distinct_ld" else index + 1 for index in protein_causal]
    realized_ld = float(np.mean([cis_ld[left, right] for left, right in zip(protein_causal, outcome_causal)]))

    protein_joint = np.zeros(len(cis_variants))
    protein_sizes = np.array([0.20, -0.18, 0.16] + [0.14] * max(0, blocks - 3))[:blocks]
    protein_joint[protein_causal] = protein_sizes
    protein_cis_mean = cis_ld @ protein_joint

    beta2 = 0.55
    outcome_joint = np.zeros(len(cis_variants))
    if scenario in {"shared_mediation", "same_variant_pleiotropy"}:
        outcome_joint[protein_causal] = beta2 * protein_sizes
    elif scenario == "distinct_ld":
        outcome_joint[outcome_causal] = beta2 * protein_sizes
    outcome_cis_mean = cis_ld @ outcome_joint

    gamma_rf = np.array([0.20, -0.18, 0.17, -0.16, 0.15, -0.14])
    beta1 = 0.60
    beta3 = 0.10
    protein_rf_mean = beta1 * gamma_rf
    if scenario in {"shared_mediation", "same_variant_pleiotropy", "distinct_ld"}:
        outcome_rf_mean = beta2 * protein_rf_mean + beta3 * gamma_rf
    else:
        outcome_rf_mean = beta3 * gamma_rf

    rf_mean = np.concatenate((np.zeros(len(cis_variants)), gamma_rf))
    protein_mean = np.concatenate((protein_cis_mean, protein_rf_mean))
    outcome_mean = np.concatenate((outcome_cis_mean, outcome_rf_mean))
    rf_se_value, protein_se_value, outcome_se_value = 0.02, 0.02, 0.02
    rf_effect = rf_mean + sample_correlated_error(rng, target_ld, rf_se_value)
    protein_effect = protein_mean + sample_correlated_error(rng, target_ld, protein_se_value)
    outcome_effect = outcome_mean + sample_correlated_error(rng, target_ld, outcome_se_value)
    rf_effect[-6:] = gamma_rf + rng.normal(0.0, 0.004, size=6)

    standard_errors = np.full(len(variants), rf_se_value)
    frequencies = reference_genotypes.mean(axis=0) / 2.0
    write_sumstats(directory / "rf.txt", variants, frequencies, rf_effect, standard_errors)
    write_sumstats(directory / "outcome.txt", variants, frequencies, outcome_effect, np.full(len(variants), outcome_se_value))
    write_protein_gwas(
        directory / "protein.glm.linear",
        variants,
        frequencies,
        protein_effect,
        np.full(len(variants), protein_se_value),
    )
    (directory / "manifest.txt").write_text(f"P1 {directory / 'protein.glm.linear'}\n", encoding="ascii")
    (directory / "protein_info.txt").write_text(
        "PROTEIN GENE CHR START END\nP1 G1 1 90000 410000\n", encoding="ascii"
    )
    return directory / "ldref", realized_ld


def read_result(path: Path) -> dict[str, str]:
    with path.open(encoding="ascii") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"expected one result row in {path}, found {len(rows)}")
    return rows[0]


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


def run_one(binary: Path, directory: Path, ld_prefix: Path) -> dict[str, str]:
    output = directory / "bmediator"
    command = [
        str(binary),
        "--rf-sumstat", str(directory / "rf.txt"),
        "--protein-gwas-list", str(directory / "manifest.txt"),
        "--cancer-sumstat", str(directory / "outcome.txt"),
        "--protein-info", str(directory / "protein_info.txt"),
        "--bfile", str(ld_prefix),
        "--out", str(output),
        "--structural-method", "factorized",
        "--p-thresh-rf", "5e-6",
        "--p-thresh-cis", "5e-6",
        "--cis-window", "500",
        "--clump-r2", "0.05",
        "--min-instruments", "1",
        "--factor-min-set-a", "3",
        "--factor-min-set-b", "3",
        "--heidi-off",
        "--no-steiger",
        "--regional-signal-p", "5e-6",
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
    return read_result(output.with_suffix(".mediation"))


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault((str(row["cell"]), str(row["scenario"])), []).append(row)
    fields = [
        "cell", "scenario", "n", "mean_realized_causal_pair_ld", "mean_nA", "mean_nB",
        "mean_cross_set_max_r2",
        "finite_two_leg_rate", "stable_effect_rate", "numerical_failure_rate",
        "two_leg_p_rate", "two_leg_BF_rate", "two_leg_e_rate", "regional_shared_rate",
        "regional_single_shared_rate", "regional_distinct_rate",
        "mean_independent_shared_signals", "my_heterogeneity_BF_rate",
        "factor_supported_rate", "factor_unresolved_single_rate",
        "factor_unresolved_heterogeneity_rate", "factor_rejected_distinct_rate",
        "frequentist_supported_rate", "frequentist_rejected_distinct_rate",
        "ebh_supported_rate", "ebh_rejected_distinct_rate",
    ]
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (cell, scenario), group in sorted(grouped.items()):
            finite_two_leg = [
                math.isfinite(float(row["factor_conjunction_p"]))
                and math.isfinite(float(row["factor_min_log_BF"]))
                and math.isfinite(float(row["factor_log_e_mediation"]))
                for row in group
            ]
            writer.writerow({
                "cell": cell,
                "scenario": scenario,
                "n": len(group),
                "mean_realized_causal_pair_ld": np.mean([float(row["realized_causal_pair_ld"]) for row in group]),
                "mean_nA": np.mean([float(row["factor_nA"]) for row in group]),
                "mean_nB": np.mean([float(row["factor_nB"]) for row in group]),
                "mean_cross_set_max_r2": np.mean([
                    float(row["factor_cross_set_max_r2"]) for row in group
                ]),
                "finite_two_leg_rate": np.mean(finite_two_leg),
                "stable_effect_rate": np.mean([
                    math.isfinite(float(row["factor_beta1"]))
                    and math.isfinite(float(row["factor_beta2"]))
                    for row in group
                ]),
                "numerical_failure_rate": np.mean([
                    row["factor_two_stage_status"] == "NUMERICAL_FAILURE"
                    for row in group
                ]),
                "two_leg_p_rate": np.mean([float(row["factor_conjunction_p"]) <= 0.05 for row in group]),
                "two_leg_BF_rate": np.mean([float(row["factor_min_log_BF"]) >= math.log(10.0) for row in group]),
                "two_leg_e_rate": np.mean([float(row["factor_log_e_mediation"]) >= math.log(20.0) for row in group]),
                "regional_shared_rate": np.mean([row["mediation_identifiability"] == "OVERIDENTIFIED_SHARED_SIGNALS_ASSUMPTION_CONDITIONAL" for row in group]),
                "regional_single_shared_rate": np.mean([row["mediation_identifiability"] == "UNRESOLVED_SINGLE_SHARED_SIGNAL" for row in group]),
                "regional_distinct_rate": np.mean([row["mediation_identifiability"] == "LD_DISTINCT_SUPPORTED" for row in group]),
                "mean_independent_shared_signals": np.mean([float(row["regional_independent_shared_signals"]) for row in group]),
                "my_heterogeneity_BF_rate": np.mean([float(row["factor_log_BF_heterogeneity_MY"]) >= math.log(10.0) for row in group]),
                "factor_supported_rate": np.mean([row["factor_mediation_status"] == "SUPPORTED_EXCLUSION_RESTRICTION_CONDITIONAL" for row in group]),
                "factor_unresolved_single_rate": np.mean([row["factor_mediation_status"] == "UNRESOLVED_SINGLE_SHARED_SIGNAL" for row in group]),
                "factor_unresolved_heterogeneity_rate": np.mean([row["factor_mediation_status"] == "UNRESOLVED_MY_HETEROGENEITY" for row in group]),
                "factor_rejected_distinct_rate": np.mean([row["factor_mediation_status"] == "REJECTED_DISTINCT_REGIONAL_SIGNALS" for row in group]),
                "frequentist_supported_rate": np.mean([
                    row["factor_frequentist_status"] in {
                        "SUPPORTED_EXCLUSION_RESTRICTION_CONDITIONAL",
                        "SUPPORTED_SCALAR_DISPERSION_EXCLUSION_CONDITIONAL",
                    }
                    for row in group
                ]),
                "frequentist_rejected_distinct_rate": np.mean([row["factor_frequentist_status"] == "REJECTED_DISTINCT_REGIONAL_SIGNALS" for row in group]),
                "ebh_supported_rate": np.mean([row["factor_ebh_status"] == "SUPPORTED_EXCLUSION_RESTRICTION_CONDITIONAL" for row in group]),
                "ebh_rejected_distinct_rate": np.mean([row["factor_ebh_status"] == "REJECTED_DISTINCT_REGIONAL_SIGNALS" for row in group]),
            })


def main() -> None:
    args = parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        raise SystemExit(f"binary not found: {binary}")
    if args.replicates < 1 or args.blocks < 2 or args.variants_per_block < 4:
        raise SystemExit("replicates >= 1, blocks >= 2, and variants-per-block >= 4 are required")
    args.outdir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "replicates": args.replicates,
        "cells": {cell: LD_CELLS[cell] for cell in args.cells},
        "scenarios": SCENARIOS,
        "reference_samples": args.reference_samples,
        "analysis_samples": args.analysis_samples,
        "blocks": args.blocks,
        "variants_per_block": args.variants_per_block,
        "seed": args.seed,
    }
    (args.outdir / "factorized_ld_stress_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="ascii"
    )

    work = args.outdir / "work"
    rows: list[dict[str, object]] = []
    for cell_index, cell in enumerate(args.cells):
        target_rho, reference_rho = LD_CELLS[cell]
        for scenario_index, scenario in enumerate(SCENARIOS):
            for replicate in range(1, args.replicates + 1):
                seed = args.seed + cell_index * 1_000_000 + scenario_index * 100_000 + replicate
                rng = np.random.default_rng(seed)
                directory = work / cell / scenario / f"rep_{replicate:04d}"
                ld_prefix, realized_ld = make_inputs(
                    directory,
                    rng,
                    scenario,
                    target_rho,
                    reference_rho,
                    args.analysis_samples,
                    args.reference_samples,
                    args.blocks,
                    args.variants_per_block,
                )
                result = run_one(binary, directory, ld_prefix)
                rows.append({
                    "cell": cell,
                    "scenario": scenario,
                    "replicate": replicate,
                    "seed": seed,
                    "target_rho": target_rho,
                    "reference_rho": reference_rho,
                    "realized_causal_pair_ld": realized_ld,
                    "factor_nA": as_float(result, "factor_nA"),
                    "factor_nB": as_float(result, "factor_nB"),
                    "factor_cross_set_max_r2": as_float(
                        result, "factor_cross_set_max_r2"
                    ),
                    "factor_beta1": as_float(result, "factor_beta1"),
                    "factor_beta2": as_float(result, "factor_beta2"),
                    "factor_conjunction_p": as_float(result, "factor_conjunction_p"),
                    "factor_min_log_BF": as_float(result, "factor_min_log_BF"),
                    "factor_log_e_mediation": as_float(result, "factor_log_e_mediation"),
                    "factor_log_BF_heterogeneity_MY": as_float(
                        result, "factor_log_BF_heterogeneity_MY"
                    ),
                    "regional_PP_shared": as_float(result, "regional_PP_shared"),
                    "regional_PP_distinct": as_float(result, "regional_PP_distinct"),
                    "regional_independent_shared_signals": as_float(
                        result, "regional_independent_shared_signals"
                    ),
                    "mediation_identifiability": result.get("mediation_identifiability", ""),
                    "factor_two_stage_status": result.get("factor_two_stage_status", ""),
                    "factor_mediation_status": result.get("factor_mediation_status", ""),
                    "factor_frequentist_status": result.get("factor_frequentist_status", ""),
                    "factor_ebh_status": result.get("factor_ebh_status", ""),
                })

    raw_path = args.outdir / "factorized_ld_stress.tsv"
    with raw_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    write_summary(args.outdir / "factorized_ld_stress_summary.tsv", rows)
    if not args.keep_inputs:
        shutil.rmtree(work)
    print(f"Wrote {raw_path} ({len(rows)} runs)")


if __name__ == "__main__":
    main()
