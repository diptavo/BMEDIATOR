#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) {
    stop("usage: build_joint_graph_v021_family_manifest.R [OUTPUT [REPLICATES [PROTEINS]]]")
}
output <- if (length(args) >= 1L) args[[1]] else
    "research/joint_graph_v0_2_1_family_manifest.tsv"
replicates <- if (length(args) >= 2L) as.integer(args[[2]]) else 50L
proteins <- if (length(args) >= 3L) as.integer(args[[3]]) else 100L
if (!is.finite(replicates) || replicates < 1L ||
    !is.finite(proteins) || proteins < 1L || proteins > 100L) {
    stop("REPLICATES and PROTEINS are invalid")
}
scenarios <- c(
    "baseline", "rare", "composite_null", "mixed", "strong_ld",
    "ld_mismatch", "scale_uncertainty", "undeclared_overlap",
    "weak_paths", "aligned_sensitivity"
)
manifest <- expand.grid(
    replicate = seq_len(replicates),
    scenario = scenarios,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
manifest <- manifest[order(match(manifest$scenario, scenarios),
                           manifest$replicate), c("scenario", "replicate")]
manifest$proteins <- proteins
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(manifest, output, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", nrow(manifest), "held-out family tasks to", output, "\n")
