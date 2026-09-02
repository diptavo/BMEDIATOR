#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import os


MISSING = {"", ".", "NA", "NaN", "nan", None}


def open_text(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def norm_chrom(value):
    return str(value).replace("chr", "").replace("CHR", "")


def is_snp(a1, a2):
    return len(a1) == 1 and len(a2) == 1 and a1 in "ACGT" and a2 in "ACGT"


def complement(a):
    return {"A": "T", "T": "A", "C": "G", "G": "C"}[a]


def allele_key(a1, a2):
    return tuple(sorted((a1.upper(), a2.upper())))


def parse_float(value):
    if value in MISSING:
        raise ValueError("missing numeric value")
    return float(value)


def finite_number(value):
    return not math.isnan(value) and not math.isinf(value)


def parse_p_from_lp(value):
    lp = parse_float(value)
    if lp <= 0:
        return 1.0
    if lp > 323:
        return 1e-323
    return 10.0 ** (-lp)


def load_bim(path):
    by_rsid = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split()
            if len(fields) < 6:
                continue
            chrom, rsid, _, bp, a1, a2 = fields[:6]
            if rsid == "." or rsid in by_rsid:
                continue
            by_rsid[rsid] = (norm_chrom(chrom), int(bp), allele_key(a1, a2))
    return by_rsid


def alleles_match(a1, a2, bim_key):
    key = allele_key(a1, a2)
    comp_key = allele_key(complement(a1), complement(a2))
    return bim_key in (key, comp_key)


def bladder_rows(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                yield {
                    "snp": row["rsid"],
                    "a1": row["allele1"].upper(),
                    "a2": row["allele2"].upper(),
                    "freq": parse_float(row["freq1"]),
                    "beta": parse_float(row["effect"]),
                    "se": parse_float(row["std_err"]),
                    "p": parse_float(row["p_value"]),
                    "chr": norm_chrom(row["chromosome"]),
                    "bp": int(row["position"]),
                }
            except (KeyError, ValueError):
                continue


def parse_vcf_sample(format_keys, sample):
    values = sample.split(":")
    parsed = dict(zip(format_keys, values))
    return {key: value.split(",")[0] for key, value in parsed.items()}


def ieu_vcf_rows(path):
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, rsid, ref, alt = fields[:5]
            if "," in alt:
                continue
            fmt = fields[8].split(":")
            sample = parse_vcf_sample(fmt, fields[9])
            try:
                yield {
                    "snp": rsid,
                    "a1": alt.upper(),
                    "a2": ref.upper(),
                    "freq": parse_float(sample["AF"]),
                    "beta": parse_float(sample["ES"]),
                    "se": parse_float(sample["SE"]),
                    "p": parse_p_from_lp(sample["LP"]),
                    "chr": norm_chrom(chrom),
                    "bp": int(pos),
                }
            except (KeyError, ValueError):
                continue


def reader_for(source):
    if source == "bladder_csv":
        return bladder_rows
    if source == "ieu_vcf":
        return ieu_vcf_rows
    raise ValueError(f"unsupported source: {source}")


def write_bmediator(source, input_path, output_path, bim_path, min_af):
    bim = load_bim(bim_path) if bim_path else None
    stats = {
        "total_rows": 0,
        "kept_rows": 0,
        "removed_non_snp": 0,
        "removed_af": 0,
        "removed_nonfinite": 0,
        "removed_no_rsid_match": 0,
    }
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w") as out:
        out.write("SNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tCHR\tBP\n")
        for row in reader_for(source)(input_path):
            stats["total_rows"] += 1
            snp = row["snp"]
            a1 = row["a1"]
            a2 = row["a2"]
            if not snp.startswith("rs") or not is_snp(a1, a2):
                stats["removed_non_snp"] += 1
                continue
            freq = float(row["freq"])
            beta = float(row["beta"])
            se = float(row["se"])
            p = float(row["p"])
            if not all(finite_number(x) for x in [freq, beta, se, p]) or se <= 0 or p <= 0 or p > 1:
                stats["removed_nonfinite"] += 1
                continue
            if not (min_af < freq < 1.0 - min_af):
                stats["removed_af"] += 1
                continue
            chrom = row["chr"]
            bp = int(row["bp"])
            if bim is not None:
                bim_row = bim.get(snp)
                if not bim_row or not alleles_match(a1, a2, bim_row[2]):
                    stats["removed_no_rsid_match"] += 1
                    continue
                chrom = bim_row[0]
                bp = bim_row[1]
            out.write(f"{snp}\t{a1}\t{a2}\t{freq:.6g}\t{beta:.6g}\t{se:.6g}\t{p:.6g}\t{chrom}\t{bp}\n")
            stats["kept_rows"] += 1
    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, choices=["bladder_csv", "ieu_vcf"])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bim", default="")
    parser.add_argument("--min-af", type=float, default=0.01)
    args = parser.parse_args()
    stats = write_bmediator(args.source, args.input, args.output, args.bim, args.min_af)
    print(f"source\t{args.source}")
    print(f"input\t{args.input}")
    print(f"output\t{args.output}")
    for key, value in stats.items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
