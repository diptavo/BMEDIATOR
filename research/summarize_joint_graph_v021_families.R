#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("usage: summarize_joint_graph_v021_families.R INPUT_ROOT OUTPUT.tsv")
}
paths <- list.files(args[[1]], pattern = "summary\\.tsv$", recursive = TRUE,
                    full.names = TRUE)
if (!length(paths)) stop("no family summary files found")
families <- do.call(rbind, lapply(paths, read.delim, stringsAsFactors = FALSE,
                                 check.names = FALSE))
expected <- 10L * 50L
if (nrow(families) != expected) {
    stop("expected ", expected, " family summaries, found ", nrow(families))
}

scenario_order <- c(
    "baseline", "rare", "composite_null", "mixed", "strong_ld",
    "ld_mismatch", "scale_uncertainty", "undeclared_overlap",
    "weak_paths", "aligned_sensitivity"
)
summaries <- do.call(rbind, lapply(scenario_order, function(name) {
    x <- families[families$scenario == name, ]
    data.frame(
        scenario = name,
        families = nrow(x),
        complete_families = sum(x$successful == x$proteins),
        failed_proteins = sum(x$proteins - x$successful),
        max_failed_proteins_per_family = max(x$proteins - x$successful),
        mean_bfdr05_discoveries = mean(x$bfdr05_discoveries),
        mean_bfdr05_false_discoveries = mean(x$bfdr05_false_discoveries),
        families_with_bfdr05_discovery = mean(x$bfdr05_discoveries > 0),
        mean_bfdr05_fdr = mean(x$bfdr05_fdp),
        mean_bfdr05_power = mean(x$bfdr05_power, na.rm = TRUE),
        mean_pp80_fdr = mean(x$pp80_fdp),
        mean_pp80_power = mean(x$pp80_power, na.rm = TRUE),
        median_protein_seconds = median(x$median_seconds),
        stringsAsFactors = FALSE
    )
}))

criterion <- function(name, expression) {
    row <- summaries[summaries$scenario == name, , drop = FALSE]
    isTRUE(with(row, eval(parse(text = expression))))
}
decisions <- data.frame(
    criterion = c(
        "baseline", "rare", "composite_null", "mixed", "strong_ld",
        "confirmatory_completion"
    ),
    passed = c(
        criterion("baseline", "mean_bfdr05_fdr <= 0.05 & mean_bfdr05_power >= 0.70"),
        criterion("rare", "mean_bfdr05_fdr <= 0.05 & mean_bfdr05_power >= 0.50"),
        criterion("composite_null", paste0(
            "mean_bfdr05_false_discoveries <= 0.05 & ",
            "families_with_bfdr05_discovery <= 0.05"
        )),
        criterion("mixed", "mean_bfdr05_fdr <= 0.05 & mean_bfdr05_power >= 0.60"),
        criterion("strong_ld", "mean_bfdr05_fdr <= 0.05 & mean_bfdr05_power >= 0.65"),
        all(summaries$complete_families[
            summaries$scenario %in% c("baseline", "rare", "composite_null",
                                      "mixed", "strong_ld")
        ] >= 49L)
    ),
    stringsAsFactors = FALSE
)

dir.create(dirname(args[[2]]), recursive = TRUE, showWarnings = FALSE)
write.table(summaries, args[[2]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(decisions, sub("\\.tsv$", "_decision.tsv", args[[2]]), sep = "\t",
            quote = FALSE, row.names = FALSE)
incomplete <- families[families$successful != families$proteins, ]
write.table(incomplete, sub("\\.tsv$", "_incomplete.tsv", args[[2]]), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(summaries, row.names = FALSE)
print(decisions, row.names = FALSE)
if (!all(decisions$passed)) {
    message("One or more frozen held-out criteria failed.")
} else {
    message("All frozen held-out family criteria passed.")
}
