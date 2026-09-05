#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else
    "research/joint_graph_v0_2_1_family_manifest.tsv"
scenarios <- c(
    "baseline", "rare", "composite_null", "mixed", "strong_ld",
    "ld_mismatch", "scale_uncertainty", "undeclared_overlap",
    "weak_paths", "aligned_sensitivity"
)
manifest <- expand.grid(
    replicate = seq_len(50L),
    scenario = scenarios,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
manifest <- manifest[order(match(manifest$scenario, scenarios),
                           manifest$replicate), c("scenario", "replicate")]
manifest$proteins <- 100L
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(manifest, output, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", nrow(manifest), "held-out family tasks to", output, "\n")
