#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("usage: summarize_joint_graph_numerics.R INPUT_ROOT OUTPUT.tsv")
}
paths <- list.files(args[[1]], pattern = "^proteins\\.tsv$", recursive = TRUE,
                    full.names = TRUE)
if (!length(paths)) stop("no protein result files found")
results <- do.call(rbind, lapply(
    paths, read.delim, stringsAsFactors = FALSE, check.names = FALSE
))
required <- c(
    "scenario", "success", "diagnostic", "elapsed_seconds",
    "estimated_quadrature_posterior_error", "max_quadrature_order",
    "posterior_aware_refinements", "sparse_grid_states",
    "max_sparse_grid_level", "max_sparse_grid_cancellation",
    "max_tensor_sparse_difference"
)
missing <- setdiff(required, names(results))
if (length(missing)) stop("missing numerical columns: ", paste(missing, collapse = ", "))

finite_quantile <- function(x, probability) {
    x <- x[is.finite(x)]
    if (length(x)) unname(quantile(x, probability, names = FALSE)) else NA_real_
}
summarize_one <- function(x, scenario) {
    successful <- x[x$success, , drop = FALSE]
    data.frame(
        scenario = scenario,
        proteins = nrow(x),
        successful = nrow(successful),
        failures = sum(!x$success),
        failure_rate = mean(!x$success),
        posterior_refined = sum(successful$posterior_aware_refinements > 0),
        sparse_grid_fits = sum(successful$sparse_grid_states > 0),
        sparse_level_10_or_higher = sum(successful$max_sparse_grid_level >= 10),
        sparse_level_12 = sum(successful$max_sparse_grid_level >= 12),
        max_sparse_cancellation = finite_quantile(
            successful$max_sparse_grid_cancellation, 1.00
        ),
        max_tensor_sparse_difference = finite_quantile(
            successful$max_tensor_sparse_difference, 1.00
        ),
        median_posterior_error = finite_quantile(
            successful$estimated_quadrature_posterior_error, 0.50
        ),
        p95_posterior_error = finite_quantile(
            successful$estimated_quadrature_posterior_error, 0.95
        ),
        max_posterior_error = finite_quantile(
            successful$estimated_quadrature_posterior_error, 1.00
        ),
        median_seconds = finite_quantile(successful$elapsed_seconds, 0.50),
        p95_seconds = finite_quantile(successful$elapsed_seconds, 0.95),
        max_seconds = finite_quantile(successful$elapsed_seconds, 1.00),
        stringsAsFactors = FALSE
    )
}
scenario_order <- unique(results$scenario)
summary <- do.call(rbind, c(
    lapply(scenario_order, function(name) {
        summarize_one(results[results$scenario == name, , drop = FALSE], name)
    }),
    list(summarize_one(results, "ALL"))
))

dir.create(dirname(args[[2]]), recursive = TRUE, showWarnings = FALSE)
write.table(summary, args[[2]], sep = "\t", quote = FALSE, row.names = FALSE,
            na = "NA")
write.table(results[!results$success, ], sub("\\.tsv$", "_failures.tsv", args[[2]]),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
print(summary, row.names = FALSE)
