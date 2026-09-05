#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
    stop("usage: run_joint_graph_manifest.R BINARY MANIFEST.tsv OUTPUT_PREFIX [CORES]")
}
binary <- normalizePath(args[[1]])
manifest_path <- normalizePath(args[[2]])
output_prefix <- args[[3]]
cores <- if (length(args) == 4L) as.integer(args[[4]]) else 1L
if (!is.finite(cores) || cores < 1L) stop("CORES must be a positive integer")

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
required <- c("protein", "input", "ld")
missing <- setdiff(required, names(manifest))
if (length(missing)) stop("manifest is missing columns: ", paste(missing, collapse = ", "))
if (!"options" %in% names(manifest)) manifest$options <- ""
if (!nrow(manifest)) stop("manifest contains no proteins")
if (any(!nzchar(manifest$protein)) || anyDuplicated(manifest$protein)) {
    stop("protein identifiers must be nonempty and unique")
}

manifest_directory <- dirname(manifest_path)
resolve_path <- function(path) {
    if (!nzchar(path) || is.na(path)) return("")
    if (!grepl("^/", path)) path <- file.path(manifest_directory, path)
    normalizePath(path, mustWork = TRUE)
}
manifest$input <- vapply(manifest$input, resolve_path, character(1))
manifest$ld <- vapply(manifest$ld, resolve_path, character(1))
manifest$options <- vapply(manifest$options, resolve_path, character(1))

run_one <- function(index) {
    temporary <- tempfile(paste0("bmediator_joint_", index, "_"))
    dir.create(temporary)
    on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
    result_path <- file.path(temporary, "result.tsv")
    stderr_path <- file.path(temporary, "stderr.txt")
    command <- c("--input", manifest$input[[index]], "--ld", manifest$ld[[index]],
                 "--out", result_path)
    if (nzchar(manifest$options[[index]])) {
        command <- c(command, "--options", manifest$options[[index]])
    }
    command <- shQuote(command)
    elapsed <- system.time(status <- system2(
        binary, command, stdout = FALSE, stderr = stderr_path
    ))[["elapsed"]]
    diagnostic <- if (file.exists(stderr_path)) {
        paste(readLines(stderr_path, warn = FALSE), collapse = " ")
    } else {
        ""
    }
    success <- identical(status, 0L) && file.exists(result_path)
    if (!success) {
        return(list(success = FALSE, protein = manifest$protein[[index]],
                    exit_status = as.integer(status), diagnostic = diagnostic,
                    elapsed_seconds = elapsed))
    }
    result <- read.delim(result_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
    if (nrow(result) != 1L || !identical(result$model_version, "JG-0.2.3")) {
        return(list(success = FALSE, protein = manifest$protein[[index]],
                    exit_status = 1L, diagnostic = "invalid result schema or model version",
                    elapsed_seconds = elapsed))
    }
    result <- cbind(data.frame(protein = manifest$protein[[index]],
                               elapsed_seconds = elapsed,
                               stringsAsFactors = FALSE), result)
    list(success = TRUE, result = result)
}

indices <- seq_len(nrow(manifest))
if (.Platform$OS.type == "unix" && cores > 1L) {
    fitted <- parallel::mclapply(indices, run_one, mc.cores = cores,
                                mc.preschedule = FALSE)
} else {
    fitted <- lapply(indices, run_one)
}

successful <- vapply(fitted, function(x) isTRUE(x$success), logical(1))
failures <- lapply(fitted[!successful], function(x) {
    data.frame(protein = x$protein, exit_status = x$exit_status,
               diagnostic = x$diagnostic, elapsed_seconds = x$elapsed_seconds,
               stringsAsFactors = FALSE)
})
failure_table <- if (length(failures)) do.call(rbind, failures) else data.frame(
    protein = character(), exit_status = integer(), diagnostic = character(),
    elapsed_seconds = numeric(), stringsAsFactors = FALSE
)

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
write.table(failure_table, paste0(output_prefix, ".failures.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
if (!any(successful)) stop("all joint-model analyses failed; see failure table")

results <- do.call(rbind, lapply(fitted[successful], `[[`, "result"))
results$posterior_lfdr <- 1 - results$PP_two_path
order_index <- order(results$posterior_lfdr, results$protein)
cumulative_fdr <- cumsum(results$posterior_lfdr[order_index]) / seq_along(order_index)
acceptable <- which(cumulative_fdr <= 0.05)
results$selected_bfdr05 <- FALSE
results$posterior_rank <- NA_integer_
results$posterior_rank[order_index] <- seq_along(order_index)
if (length(acceptable)) {
    results$selected_bfdr05[order_index[seq_len(max(acceptable))]] <- TRUE
}
results <- results[match(manifest$protein[manifest$protein %in% results$protein],
                         results$protein), ]
write.table(results, paste0(output_prefix, ".joint.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("Completed", nrow(results), "proteins;", nrow(failure_table),
    "failed;", sum(results$selected_bfdr05), "selected at posterior FDR 5%\n")
