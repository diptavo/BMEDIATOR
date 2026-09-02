#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(SomaScan.db)
  library(EnsDb.Hsapiens.v86)
})

sumstat_dir <- "/data/BB_Bioinformatics/ProjectData/DeCode_protein_sumstat/sumstats"
out_dir <- "/data/Dutta_lab/BMEDIATOR/generated/decode_resource"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(sumstat_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
if (length(files) == 0L) stop("No deCODE sumstat files found in: ", sumstat_dir)

collapse_unique <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  paste(x, collapse = ";")
}

first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  x[[1]]
}

base <- basename(files)
parts <- tstrsplit(base, "_", fixed = TRUE)
seqid_raw <- paste(parts[[1]], parts[[2]], sep = "_")
parsed_symbol <- parts[[3]]
prefix_nchar <- nchar(seqid_raw) + nchar(parsed_symbol) + 2L
alias_raw <- substring(base, prefix_nchar + 1L)
alias_clean <- sub("\\.txt\\.gz$", "", alias_raw)

meta <- data.table(
  sumstat_file = files,
  file_basename = base,
  seqid_raw = seqid_raw,
  probeid = gsub("_", "-", seqid_raw, fixed = TRUE),
  seqid = paste0("SeqId_", seqid_raw),
  parsed_symbol = parsed_symbol,
  parsed_alias = alias_clean
)

probeids <- unique(meta$probeid)

ann <- as.data.table(select(
  SomaScan.db,
  keys = probeids,
  keytype = "PROBEID",
  columns = c("SYMBOL", "GENENAME", "UNIPROT", "MAP", "ENSEMBL", "ENTREZID")
))

full_name_map <- as.list(SomaScanTARGETFULLNAME[probeids])
full_name_dt <- data.table(
  probeid = names(full_name_map),
  target_full_name = vapply(full_name_map, first_nonempty, character(1))
)

ann_agg <- ann[, .(
  package_symbol = first_nonempty(SYMBOL),
  gene_name = first_nonempty(GENENAME),
  uniprot_ids = collapse_unique(UNIPROT),
  cytoband = first_nonempty(MAP),
  ensembl_gene_ids = collapse_unique(ENSEMBL),
  entrez_ids = collapse_unique(ENTREZID),
  n_somascan_rows = .N
), by = PROBEID]
setnames(ann_agg, "PROBEID", "probeid")

res <- merge(meta, ann_agg, by = "probeid", all.x = TRUE)
res <- merge(res, full_name_dt, by = "probeid", all.x = TRUE)
res[, symbol := fifelse(!is.na(package_symbol) & nzchar(package_symbol), package_symbol, parsed_symbol)]

coord_keys <- unique(na.omit(res$symbol))
coord_dt <- data.table()
if (length(coord_keys) > 0L) {
  coord_raw <- as.data.table(select(
    EnsDb.Hsapiens.v86,
    keys = coord_keys,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "GENEID", "GENENAME", "SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQSTRAND")
  ))
  std_chr <- c(as.character(1:22), "X", "Y", "MT")
  coord_raw <- coord_raw[SEQNAME %in% std_chr]
  setorder(coord_raw, SYMBOL, SEQNAME, GENESEQSTART, GENEID)
  coord_dt <- coord_raw[, .(
    gene_id_hg38 = first_nonempty(GENEID),
    gene_name_hg38 = first_nonempty(GENENAME),
    chromosome_hg38 = first_nonempty(SEQNAME),
    gene_start_hg38 = suppressWarnings(as.integer(first_nonempty(GENESEQSTART))),
    gene_end_hg38 = suppressWarnings(as.integer(first_nonempty(GENESEQEND))),
    strand_hg38 = suppressWarnings(as.integer(first_nonempty(SEQSTRAND))),
    n_coord_rows = .N
  ), by = SYMBOL]
  coord_dt[, tss_hg38 := fifelse(strand_hg38 == -1L, gene_end_hg38, gene_start_hg38)]
  setnames(coord_dt, "SYMBOL", "symbol")
}

res <- merge(res, coord_dt, by = "symbol", all.x = TRUE)
res[, mapping_source := "SomaScan.db + EnsDb.Hsapiens.v86"]
setcolorder(res, c(
  "seqid", "seqid_raw", "probeid", "sumstat_file", "file_basename",
  "parsed_symbol", "parsed_alias", "package_symbol", "symbol",
  "gene_name", "target_full_name", "gene_name_hg38", "uniprot_ids", "entrez_ids",
  "ensembl_gene_ids", "gene_id_hg38", "cytoband",
  "chromosome_hg38", "gene_start_hg38", "gene_end_hg38", "strand_hg38", "tss_hg38",
  "n_somascan_rows", "n_coord_rows", "mapping_source"
))

out_tsv <- file.path(out_dir, "decode_seqid_to_gene_uniprot_tss_hg38.tsv")
out_rds <- file.path(out_dir, "decode_seqid_to_gene_uniprot_tss_hg38.rds")
out_summary <- file.path(out_dir, "decode_seqid_mapping_summary.txt")
out_session <- file.path(out_dir, "decode_seqid_mapping_sessionInfo.txt")
out_missing_tss <- file.path(out_dir, "decode_seqid_missing_hg38_tss.tsv")
out_missing_uniprot <- file.path(out_dir, "decode_seqid_missing_uniprot.tsv")

fwrite(res, out_tsv, sep = "\t", na = "NA")
saveRDS(res, out_rds)
fwrite(res[is.na(tss_hg38)], out_missing_tss, sep = "\t", na = "NA")
fwrite(res[is.na(uniprot_ids) | !nzchar(uniprot_ids)], out_missing_uniprot, sep = "\t", na = "NA")

summary_lines <- c(
  sprintf("n_sumstat_files\t%d", nrow(meta)),
  sprintf("n_unique_seqids\t%d", uniqueN(meta$seqid)),
  sprintf("n_mapped_symbol\t%d", sum(!is.na(res$symbol) & nzchar(res$symbol))),
  sprintf("n_mapped_target_full_name\t%d", sum(!is.na(res$target_full_name) & nzchar(res$target_full_name))),
  sprintf("n_mapped_uniprot\t%d", sum(!is.na(res$uniprot_ids) & nzchar(res$uniprot_ids))),
  sprintf("n_mapped_hg38_tss\t%d", sum(!is.na(res$tss_hg38))),
  sprintf("n_missing_symbol\t%d", sum(is.na(res$symbol) | !nzchar(res$symbol))),
  sprintf("n_missing_hg38_tss\t%d", sum(is.na(res$tss_hg38))),
  sprintf("n_missing_uniprot\t%d", sum(is.na(res$uniprot_ids) | !nzchar(res$uniprot_ids)))
)
writeLines(summary_lines, out_summary)
writeLines(capture.output(sessionInfo()), out_session)

cat("Wrote:\n")
cat(out_tsv, "\n", sep = "")
cat(out_rds, "\n", sep = "")
cat(out_summary, "\n", sep = "")
cat(out_session, "\n", sep = "")
cat(out_missing_tss, "\n", sep = "")
cat(out_missing_uniprot, "\n", sep = "")
