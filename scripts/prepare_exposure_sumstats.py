#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import os
import sys
import tarfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, os.path.join(REPO_DIR, ".deps"))
from pyliftover import LiftOver


def open_text(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def open_tar_member_text(path):
    tf = tarfile.open(path, "r:gz")
    member = next((m for m in tf.getmembers() if m.isfile()), None)
    if member is None:
        tf.close()
        raise ValueError(f"no regular file found in tar archive: {path}")
    fh = tf.extractfile(member)
    if fh is None:
        tf.close()
        raise ValueError(f"failed to extract {member.name} from {path}")
    return tf, fh


def is_snp(a1, a2):
    return len(a1) == 1 and len(a2) == 1 and a1 in "ACGT" and a2 in "ACGT"


def complement(a):
    return {"A": "T", "T": "A", "C": "G", "G": "C"}[a]


def allele_key(a1, a2):
    return tuple(sorted((a1.upper(), a2.upper())))


def parse_p_from_neglog10(value):
    x = float(value)
    if x <= 0:
        return 1.0
    return 10.0 ** (-x)


def finite_number(x):
    return x is not None and not math.isnan(x) and not math.isinf(x)


def parse_float(value):
    if value in ("", "NA", None):
        raise ValueError("missing numeric value")
    return float(value)


def panukb_rows(path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                if row.get("low_confidence_EUR", "").lower() == "true":
                    continue
                if row.get("beta_EUR") in ("", "NA") or row.get("se_EUR") in ("", "NA"):
                    continue
                if row.get("neglog10_pval_EUR") in ("", "NA"):
                    continue
                a1 = row["alt"].upper()
                a2 = row["ref"].upper()
                freq = parse_float(row["af_EUR"])
                beta = parse_float(row["beta_EUR"])
                se = parse_float(row["se_EUR"])
                p = parse_p_from_neglog10(row["neglog10_pval_EUR"])
                yield {
                    "snp": f'{row["chr"]}:{row["pos"]}:{a2}:{a1}',
                    "a1": a1,
                    "a2": a2,
                    "freq": freq,
                    "beta": beta,
                    "se": se,
                    "p": p,
                    "chr": row["chr"],
                    "bp": int(row["pos"]),
                }
            except (KeyError, ValueError):
                continue


def giant_rows(path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                yield {
                    "snp": row["MarkerName"],
                    "a1": row["Allele1"].upper(),
                    "a2": row["Allele2"].upper(),
                    "freq": parse_float(row["FreqAllele1HapMapCEU"]),
                    "beta": parse_float(row["b"]),
                    "se": parse_float(row["se"]),
                    "p": parse_float(row["p"]),
                    "chr": row["Chr"],
                    "bp": int(row["Pos"]),
                }
            except (KeyError, ValueError):
                continue


def giant_legacy_rows(path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                yield {
                    "snp": row["SNP"],
                    "a1": row["A1"].upper(),
                    "a2": row["A2"].upper(),
                    "freq": parse_float(row["Freq1.Hapmap"]),
                    "beta": parse_float(row["b"]),
                    "se": parse_float(row["se"]),
                    "p": parse_float(row["p"]),
                    "chr": "",
                    "bp": 0,
                }
            except (KeyError, ValueError):
                continue


def giant_tar_rows(path):
    tf, raw_fh = open_tar_member_text(path)
    try:
        f = (line.decode("utf-8") for line in raw_fh)
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                chrom = row.get("Chr", "")
                bp = row.get("Pos", "")
                yield {
                    "snp": row["MarkerName"],
                    "a1": row["Allele1"].upper(),
                    "a2": row["Allele2"].upper(),
                    "freq": parse_float(row["FreqAllele1HapMapCEU"]),
                    "beta": parse_float(row["b"]),
                    "se": parse_float(row["se"]),
                    "p": parse_float(row["p"]),
                    "chr": chrom,
                    "bp": int(bp) if bp not in ("", "NA") else 0,
                }
            except (KeyError, ValueError):
                continue
    finally:
        raw_fh.close()
        tf.close()


def glgc_rows(path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                yield {
                    "snp": row["rsID"],
                    "a1": row["ALT"].upper(),
                    "a2": row["REF"].upper(),
                    "freq": parse_float(row["POOLED_ALT_AF"]),
                    "beta": parse_float(row["EFFECT_SIZE"]),
                    "se": parse_float(row["SE"]),
                    "p": parse_float(row["pvalue"]),
                    "chr": row["CHROM"],
                    "bp": int(row["POS_b37"]),
                }
            except (KeyError, ValueError):
                continue


def gwas_catalog_harmonized_rows(path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                snp = row.get("rsid") or row.get("rs_id") or ""
                if not snp:
                    continue
                yield {
                    "snp": snp,
                    "a1": row["effect_allele"].upper(),
                    "a2": row["other_allele"].upper(),
                    "freq": parse_float(row["effect_allele_frequency"]),
                    "beta": parse_float(row["beta"]),
                    "se": parse_float(row["standard_error"]),
                    "p": parse_float(row["p_value"]),
                    "chr": row["chromosome"],
                    "bp": int(row["base_pair_location"]),
                }
            except (KeyError, ValueError):
                continue


def finngen_rows(path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                rsids = row["rsids"].split(",")
                snp = next((x for x in rsids if x.startswith("rs")), "")
                if not snp:
                    continue
                yield {
                    "snp": snp,
                    "a1": row["alt"].upper(),
                    "a2": row["ref"].upper(),
                    "freq": parse_float(row["af_alt"]),
                    "beta": parse_float(row["beta"]),
                    "se": parse_float(row["sebeta"]),
                    "p": parse_float(row["pval"]),
                    "chr": row["#chrom"],
                    "bp": int(row["pos"]),
                }
            except (KeyError, ValueError):
                continue


def ckdgen_space_rows(path):
    with open_text(path) as f:
        header = f.readline().strip().split()
        for line in f:
            fields = line.strip().split()
            if len(fields) != len(header):
                continue
            row = dict(zip(header, fields))
            try:
                yield {
                    "snp": row["RSID"],
                    "a1": row["Allele1"].upper(),
                    "a2": row["Allele2"].upper(),
                    "freq": parse_float(row["Freq1"]),
                    "beta": parse_float(row["Effect"]),
                    "se": parse_float(row["StdErr"]),
                    "p": parse_float(row["P-value"]),
                    "chr": row["Chr"],
                    "bp": int(row["Pos_b37"]),
                }
            except (KeyError, ValueError):
                continue


def ckdgen_csv_rows(path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            try:
                snp = row["rsID"]
                chrom = ""
                bp = 0
                if not snp.startswith("rs"):
                    parts = snp.split(",")[0].split(":")
                    if len(parts) >= 2:
                        chrom = parts[0]
                        bp = int(parts[1])
                yield {
                    "snp": snp,
                    "a1": row["allele1"].upper(),
                    "a2": row["allele2"].upper(),
                    "freq": parse_float(row["freqA1"]),
                    "beta": parse_float(row["beta"]),
                    "se": parse_float(row["se"]),
                    "p": parse_float(row["pval"]),
                    "chr": chrom,
                    "bp": bp,
                }
            except (KeyError, ValueError):
                continue


def bmediator_rows(path):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                yield {
                    "snp": row["SNP"],
                    "a1": row["A1"].upper(),
                    "a2": row["A2"].upper(),
                    "freq": parse_float(row["FREQ"]),
                    "beta": parse_float(row["BETA"]),
                    "se": parse_float(row["SE"]),
                    "p": parse_float(row["P"]),
                    "chr": row["CHR"],
                    "bp": int(row["BP"]),
                }
            except (KeyError, ValueError):
                continue


def choose_reader(source):
    if source == "panukb":
        return panukb_rows
    if source == "giant":
        return giant_rows
    if source == "giant_legacy":
        return giant_legacy_rows
    if source == "giant_tar":
        return giant_tar_rows
    if source == "glgc":
        return glgc_rows
    if source == "gwas_catalog_harmonized":
        return gwas_catalog_harmonized_rows
    if source == "finngen":
        return finngen_rows
    if source == "ckdgen_space":
        return ckdgen_space_rows
    if source == "ckdgen_csv":
        return ckdgen_csv_rows
    if source == "bmediator":
        return bmediator_rows
    raise ValueError(f"unknown source: {source}")


def load_bim_map(path, want_coord=False, want_rsid=False):
    coord_map = {} if want_coord else None
    rsid_map = {} if want_rsid else None
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split()
            if len(fields) < 6:
                continue
            chrom, rsid, _, bp, a1, a2 = fields[:6]
            if want_coord:
                key = (chrom.replace("chr", ""), int(bp), allele_key(a1, a2))
                if key not in coord_map:
                    coord_map[key] = rsid
            if want_rsid and rsid not in rsid_map:
                rsid_map[rsid] = (chrom.replace("chr", ""), int(bp), allele_key(a1, a2))
    return coord_map, rsid_map


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--source",
        required=True,
        choices=[
            "panukb",
            "giant",
            "giant_legacy",
            "giant_tar",
            "glgc",
            "gwas_catalog_harmonized",
            "finngen",
            "ckdgen_space",
            "ckdgen_csv",
            "bmediator",
        ],
    )
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--chain", default=os.path.join(REPO_DIR, "generated", "hg19ToHg38.over.chain.gz"))
    ap.add_argument("--bim", default="")
    ap.add_argument("--liftover", action="store_true")
    ap.add_argument("--min-af", type=float, default=0.01)
    args = ap.parse_args()

    reader = choose_reader(args.source)
    lo = LiftOver(args.chain) if args.liftover else None
    if args.bim:
        coord_needed = args.source in ("panukb", "ckdgen_csv")
        rsid_needed = args.source in ("giant", "glgc", "bmediator")
        rsid_needed = rsid_needed or args.source in (
            "giant_legacy",
            "giant_tar",
            "gwas_catalog_harmonized",
            "finngen",
            "ckdgen_space",
            "ckdgen_csv",
        )
        coord_bim_map, rsid_bim_map = load_bim_map(args.bim, want_coord=coord_needed, want_rsid=rsid_needed)
    else:
        coord_bim_map, rsid_bim_map = (None, None)

    total = 0
    kept = 0
    nonsnp = 0
    af_fail = 0
    nonfinite = 0
    unmapped = 0
    multi = 0
    rsid_miss = 0

    with open(args.output, "w") as out:
        out.write("SNP\tA1\tA2\tFREQ\tBETA\tSE\tP\tCHR\tBP\n")
        for row in reader(args.input):
            total += 1
            a1 = row["a1"]
            a2 = row["a2"]
            if not is_snp(a1, a2):
                nonsnp += 1
                continue
            freq = float(row["freq"])
            beta = float(row["beta"])
            se = float(row["se"])
            p = float(row["p"])
            if not all(finite_number(x) for x in [freq, beta, se, p]) or se <= 0 or p <= 0 or p > 1:
                nonfinite += 1
                continue
            if not (args.min_af < freq < 1.0 - args.min_af):
                af_fail += 1
                continue

            chrom = str(row["chr"]).replace("chr", "")
            bp = int(row["bp"])
            if lo is not None:
                lifted = lo.convert_coordinate(f"chr{chrom}", bp - 1)
                if not lifted:
                    unmapped += 1
                    continue
                if len(lifted) != 1:
                    multi += 1
                    continue
                chrom = lifted[0][0].replace("chr", "")
                bp = int(lifted[0][1]) + 1

            snp = row["snp"]
            key = allele_key(a1, a2)
            comp_key = allele_key(complement(a1), complement(a2))

            if rsid_bim_map is not None and snp.startswith("rs"):
                bim_row = rsid_bim_map.get(snp)
                if not bim_row or bim_row[2] not in (key, comp_key):
                    rsid_miss += 1
                    continue
                chrom = bim_row[0]
                bp = bim_row[1]
            elif coord_bim_map is not None:
                coord_key = (chrom, bp, key)
                comp_coord_key = (chrom, bp, comp_key)
                snp = coord_bim_map.get(coord_key, "") or coord_bim_map.get(comp_coord_key, "")
                if not snp or snp == ".":
                    rsid_miss += 1
                    continue

            out.write(
                f'{snp}\t{a1}\t{a2}\t{freq:.6g}\t{beta:.6g}\t{se:.6g}\t{p:.6g}\t{chrom}\t{bp}\n'
            )
            kept += 1

    print(f"source\t{args.source}")
    print(f"input\t{args.input}")
    print(f"output\t{args.output}")
    print(f"total_rows\t{total}")
    print(f"kept_rows\t{kept}")
    print(f"removed_non_snp\t{nonsnp}")
    print(f"removed_af\t{af_fail}")
    print(f"removed_nonfinite\t{nonfinite}")
    if lo is not None:
        print(f"removed_unmapped\t{unmapped}")
        print(f"removed_multimapped\t{multi}")
    if coord_bim_map is not None or rsid_bim_map is not None:
        print(f"removed_no_rsid_match\t{rsid_miss}")


if __name__ == "__main__":
    main()
