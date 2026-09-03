#!/usr/bin/env python3
"""Generate a deterministic, tiny full-mode fixture including a PLINK panel."""

from __future__ import annotations

import math
import random
import struct
import sys
from pathlib import Path


VARIANTS = [
    (1, "s1", 90_000),
    (1, "s2", 100_000),
    (1, "s3", 110_000),
    (1, "s4", 120_000),
    (2, "d1", 190_000),
    (2, "d2", 200_000),
    (2, "d3", 210_000),
    (2, "d4", 220_000),
    (3, "r1", 300_000),
    (3, "r2", 400_000),
]


def two_sided_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def correlated_copy(rng: random.Random, source: list[int], flip_rate: float) -> list[int]:
    result = []
    for dosage in source:
        alleles = [1] * dosage + [0] * (2 - dosage)
        result.append(sum(1 - allele if rng.random() < flip_rate else allele for allele in alleles))
    return result


def correlation(x: list[int], y: list[int]) -> float:
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    centered_x = [value - mean_x for value in x]
    centered_y = [value - mean_y for value in y]
    covariance = sum(a * b for a, b in zip(centered_x, centered_y))
    scale = math.sqrt(sum(a * a for a in centered_x) * sum(b * b for b in centered_y))
    return covariance / scale


def write_plink(prefix: Path) -> dict[tuple[str, str], float]:
    n_samples = 512
    rng = random.Random(20260817)
    independent = lambda: [
            int(rng.random() < 0.5) + int(rng.random() < 0.5)
            for _ in range(n_samples)
        ]
    s1 = independent()
    d1 = independent()
    genotypes = [
        s1,
        correlated_copy(rng, s1, 0.18),
        independent(),
        independent(),
        d1,
        correlated_copy(rng, d1, 0.18),
        correlated_copy(rng, d1, 0.25),
        None,
        independent(),
        independent(),
    ]
    genotypes[7] = correlated_copy(rng, genotypes[6], 0.18)

    with prefix.with_suffix(".fam").open("w", encoding="ascii") as handle:
        for idx in range(n_samples):
            handle.write(f"F{idx + 1} I{idx + 1} 0 0 0 -9\n")

    with prefix.with_suffix(".bim").open("w", encoding="ascii") as handle:
        for chrom, rsid, bp in VARIANTS:
            handle.write(f"{chrom} {rsid} 0 {bp} A G\n")

    code_for_dosage = {2: 0b00, 1: 0b10, 0: 0b11}
    with prefix.with_suffix(".bed").open("wb") as handle:
        handle.write(bytes((0x6C, 0x1B, 0x01)))
        for dosage in genotypes:
            for start in range(0, n_samples, 4):
                byte = 0
                for offset, value in enumerate(dosage[start : start + 4]):
                    byte |= code_for_dosage[value] << (2 * offset)
                handle.write(struct.pack("B", byte))

    rsids = [rsid for _, rsid, _ in VARIANTS]
    return {
        (left, right): correlation(genotypes[i], genotypes[j])
        for i, left in enumerate(rsids)
        for j, right in enumerate(rsids)
    }


def write_standard_sumstats(
    path: Path,
    effects: dict[str, tuple[float, float]],
    excluded: set[str] | None = None,
) -> None:
    lookup = {rsid: (chrom, bp) for chrom, rsid, bp in VARIANTS}
    excluded = excluded or set()
    with path.open("w", encoding="ascii") as handle:
        handle.write("SNP A1 A2 FREQ BETA SE P CHR BP\n")
        for _, rsid, _ in VARIANTS:
            if rsid in excluded:
                continue
            beta, se = effects.get(rsid, (0.0, 0.02))
            chrom, bp = lookup[rsid]
            handle.write(
                f"{rsid} A G 0.5 {beta:.6g} {se:.6g} "
                f"{two_sided_p(beta / se):.8g} {chrom} {bp}\n"
            )


def write_protein_gwas(
    path: Path,
    rsids: list[str],
    effects: dict[str, tuple[float, float]],
) -> None:
    lookup = {rsid: (chrom, bp) for chrom, rsid, bp in VARIANTS}
    columns = [
        "#CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ",
        "TEST", "BETA", "SE", "P", "ERRCODE",
    ]
    with path.open("w", encoding="ascii") as handle:
        handle.write("\t".join(columns) + "\n")
        for rsid in rsids:
            beta, se = effects[rsid]
            chrom, bp = lookup[rsid]
            row = [
                str(chrom), str(bp), rsid, "G", "A", "A", "0.5", "ADD",
                f"{beta:.6g}", f"{se:.6g}",
                f"{two_sided_p(beta / se):.8g}", ".",
            ]
            handle.write("\t".join(row) + "\n")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: generate_full_mode_fixture.py OUTPUT_DIR")
    output = Path(sys.argv[1])
    output.mkdir(parents=True, exist_ok=True)
    ld = write_plink(output / "ldref")

    rf_effects = {"r1": (0.12, 0.015), "r2": (0.10, 0.015)}
    outcome_effects = {
        snp: (0.18 * ld[(snp, "s1")] + 0.16 * ld[(snp, "s3")], 0.02)
        for snp in ("s1", "s2", "s3", "s4")
    }
    outcome_effects.update({
        "d1": (0.18 * ld[("d1", "d3")], 0.02),
        "d2": (0.18 * ld[("d2", "d3")], 0.02),
        "d3": (0.18 * ld[("d3", "d3")], 0.02),
        "d4": (0.18 * ld[("d4", "d3")], 0.02),
        "r1": (0.05, 0.015), "r2": (0.04, 0.015),
    })
    write_standard_sumstats(output / "rf.txt", rf_effects)
    write_standard_sumstats(output / "rf_missing_s3.txt", rf_effects, {"s3"})
    write_standard_sumstats(
        output / "rf_no_cis.txt",
        rf_effects,
        {"s1", "s2", "s3", "s4", "d1", "d2", "d3", "d4"},
    )
    write_standard_sumstats(output / "outcome.txt", outcome_effects)

    shared_effects = {
        snp: (0.20 * ld[(snp, "s1")] + 0.17 * ld[(snp, "s3")], 0.02)
        for snp in ("s1", "s2", "s3", "s4")
    }
    shared_effects.update({
        "r1": (0.12, 0.02), "r2": (0.10, 0.02),
    })
    distinct_effects = {
        "d1": (0.20 * ld[("d1", "d1")], 0.02),
        "d2": (0.20 * ld[("d2", "d1")], 0.02),
        "d3": (0.20 * ld[("d3", "d1")], 0.02),
        "d4": (0.20 * ld[("d4", "d1")], 0.02),
        "r1": (0.12, 0.02), "r2": (0.10, 0.02),
    }
    write_protein_gwas(
        output / "shared.glm.linear", ["s1", "s2", "s3", "s4", "r1", "r2"],
        shared_effects,
    )
    write_protein_gwas(
        output / "distinct.glm.linear", ["d1", "d2", "d3", "d4", "r1", "r2"],
        distinct_effects,
    )

    with (output / "manifest.txt").open("w", encoding="ascii") as handle:
        handle.write(f"P_SHARED {output / 'shared.glm.linear'}\n")
        handle.write(f"P_DISTINCT {output / 'distinct.glm.linear'}\n")
    with (output / "manifest_single.txt").open("w", encoding="ascii") as handle:
        handle.write(f"P_SHARED {output / 'shared.glm.linear'}\n")
    with (output / "protein_info.txt").open("w", encoding="ascii") as handle:
        handle.write("PROTEIN GENE CHR START END\n")
        handle.write("P_SHARED G_SHARED 1 100000 110000\n")
        handle.write("P_DISTINCT G_DISTINCT 2 200000 210000\n")


if __name__ == "__main__":
    main()
