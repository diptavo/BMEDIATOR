#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(EnsDb.Hsapiens.v86)
  library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: build_ukbb_oid_resource.R <ukbb_dir> <out.tsv>", call. = FALSE)
}

input_source <- normalizePath(args[[1]], mustWork = TRUE)
out_tsv <- args[[2]]
cohort_label <- if (length(args) >= 3) args[[3]] else basename(input_source)

if (dir.exists(input_source)) {
  dirs <- list.dirs(input_source, full.names = TRUE, recursive = FALSE)
} else {
  dirs <- readLines(input_source, warn = FALSE)
  dirs <- dirs[nzchar(dirs)]
}
stopifnot(length(dirs) > 0)

dir_dt <- data.table(sumstat_dir = dirs)
dir_dt[, dir_name := basename(sumstat_dir)]
dir_dt[, c("gene", "uniprot", "protein", "version", "panel") := tstrsplit(dir_name, "_", fixed = TRUE, keep = 1:5)]
dir_dt[, panel := sub("^[^_]+_[^_]+_[^_]+_[^_]+_", "", dir_name)]
dir_dt <- dir_dt[!is.na(protein) & protein != ""]

edb <- EnsDb.Hsapiens.v86
gene_keys <- sort(unique(dir_dt$gene))
anno_raw <- as.data.table(AnnotationDbi::select(
  edb,
  keys = gene_keys,
  keytype = "GENENAME",
  columns = c("GENENAME", "GENEID", "SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQSTRAND")
))

std_chr <- c(as.character(1:22), "X", "Y")
anno_raw <- anno_raw[!is.na(SEQNAME) & SEQNAME %in% std_chr]
anno_raw[, strand_num := as.integer(SEQSTRAND)]
anno_raw[, tss := fifelse(strand_num == -1L, as.integer(GENESEQEND), as.integer(GENESEQSTART))]
anno_raw[, annotation_matches := .N, by = GENENAME]
anno_raw[, chr_rank := match(SEQNAME, std_chr)]
setorder(anno_raw, GENENAME, chr_rank, GENESEQSTART, GENEID)
anno_one <- anno_raw[, .SD[1], by = GENENAME]
setnames(anno_one, "GENENAME", "gene")

resource <- merge(dir_dt, anno_one, by = "gene", all.x = TRUE)
resource[, annotation_status := fifelse(is.na(SEQNAME), "missing", fifelse(annotation_matches > 1L, "multi_match_first", "unique_match"))]

out_dt <- resource[, .(
  protein,
  gene,
  uniprot,
  version,
  panel,
  sumstat_dir,
  chr = fifelse(SEQNAME == "X", "23", fifelse(SEQNAME == "Y", "24", SEQNAME)),
  tss,
  gene_start = as.integer(GENESEQSTART),
  gene_end = as.integer(GENESEQEND),
  strand = strand_num,
  gene_id = GENEID,
  annotation_matches = fifelse(is.na(annotation_matches), 0L, as.integer(annotation_matches)),
  annotation_status,
  resource_source = cohort_label
)]

setorder(out_dt, protein)
fwrite(out_dt, out_tsv, sep = "\t", quote = FALSE, na = "")
