#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop(paste(
        "usage: replay_joint_graph_failures.R",
        "BINARY OLD_RESULT_ROOT OUTPUT.tsv SEED_BASE"
    ))
}
binary <- normalizePath(args[[1]])
result_root <- normalizePath(args[[2]])
output <- args[[3]]
seed_base <- as.integer(args[[4]])
if (!is.finite(seed_base) || seed_base < 1L) stop("SEED_BASE is invalid")

root <- normalizePath(file.path(dirname(binary), ".."))
source(file.path(root, "research", "joint_graph_v02_simulation.R"))

paths <- list.files(result_root, pattern = "proteins\\.tsv$", recursive = TRUE,
                    full.names = TRUE)
if (!length(paths)) stop("no historical protein result files found")
historical <- do.call(rbind, lapply(
    paths, read.delim, stringsAsFactors = FALSE, check.names = FALSE
))
failed <- historical[!historical$success,
                     c("scenario", "replicate", "protein", "truth")]
if (!nrow(failed)) stop("historical result set contains no failed proteins")

scenario_names <- c(
    "baseline", "rare", "composite_null", "mixed", "strong_ld",
    "ld_mismatch", "scale_uncertainty", "undeclared_overlap",
    "weak_paths", "aligned_sensitivity"
)
type_arguments <- function(type) {
    switch(
        type,
        null = list(),
        mediation = list(a = 0.40, b = 0.40),
        weak_mediation = list(a = 0.20, b = 0.20, se = 0.05),
        xm_only = list(a = 0.60),
        my_only = list(b = 0.40),
        sparse = list(a = 0.60, lambda = 0.70, q = 0.35),
        directional = list(a = 0.60, eta = 0.40),
        mediation_sparse = list(blocks = 30L, a = 0.40, b = 0.40,
                                lambda = 0.70, q = 0.35),
        aligned = list(a = 0.60, lambda = 0.70, q = 1),
        stop("unknown protein type: ", type)
    )
}

temporary <- tempfile("jg025_failure_replay_")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
rows <- vector("list", nrow(failed))
for (i in seq_len(nrow(failed))) {
    x <- failed[i, ]
    scenario_index <- match(x$scenario, scenario_names)
    if (is.na(scenario_index)) stop("unknown scenario: ", x$scenario)
    protein_index <- as.integer(sub("^P", "", x$protein))
    family_seed <- seed_base + scenario_index * 100000L +
        as.integer(x$replicate) * 1000L
    simulation_arguments <- modifyList(
        list(seed = family_seed + protein_index, blocks = 20L,
             variants_per_block = 2L, ld_rho = 0.40),
        type_arguments(x$truth)
    )
    if (x$scenario %in% c("strong_ld", "ld_mismatch")) {
        simulation_arguments$ld_rho <- 0.70
    }
    if (x$scenario == "scale_uncertainty") {
        set.seed(family_seed + 500000L + protein_index)
        simulation_arguments$reported_variance_scale <- exp(rnorm(1, sd = 0.25))
    }
    if (x$scenario == "undeclared_overlap") {
        simulation_arguments$sampling_rho <- c(xm = 0.40, xy = 0.30, my = 0.20)
        simulation_arguments$reported_sampling_rho <- c(xm = 0, xy = 0, my = 0)
    }
    fixture <- do.call(jg02_simulate, simulation_arguments)
    if (x$scenario == "ld_mismatch") {
        fixture$ld <- jg02_block_ld(fixture$data$ld_block, 0.50)
    }
    input <- file.path(temporary, "input.tsv")
    ld_path <- file.path(temporary, "ld.tsv")
    result_path <- file.path(temporary, "result.tsv")
    stderr_path <- file.path(temporary, "stderr.txt")
    unlink(c(result_path, stderr_path))
    jg02_write_fixture(fixture, input, ld_path)
    timing <- system.time(status <- system2(
        binary, c(input, ld_path, result_path), stdout = FALSE,
        stderr = stderr_path
    ))
    success <- identical(status, 0L) && file.exists(result_path)
    diagnostic <- if (file.exists(stderr_path)) {
        paste(readLines(stderr_path, warn = FALSE), collapse = " ")
    } else ""
    values <- setNames(rep(NA_real_, 10L), c(
        "PP_two_path", "estimated_quadrature_posterior_error",
        "max_relevant_quadrature_difference", "max_quadrature_order",
        "posterior_aware_refinements", "sparse_grid_states",
        "max_sparse_grid_level", "max_sparse_grid_cancellation",
        "max_tensor_sparse_difference", "states_regularized"
    ))
    model_version <- NA_character_
    if (success) {
        result <- read.delim(result_path, check.names = FALSE)
        model_version <- result$model_version
        values[] <- unlist(result[names(values)], use.names = FALSE)
    }
    rows[[i]] <- data.frame(
        scenario = x$scenario, replicate = x$replicate, protein = x$protein,
        truth = x$truth, success = success, exit_status = as.integer(status),
        model_version = model_version, as.list(values),
        elapsed_seconds = timing[["elapsed"]], diagnostic = diagnostic,
        check.names = FALSE, stringsAsFactors = FALSE
    )
}

result <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(result, output, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "NA")
cat("Replayed", nrow(result), "historical failures;",
    sum(result$success), "now reportable\n")
if (!all(result$success)) quit(status = 2L)
if (any(result$model_version != "JG-0.2.7") ||
    any(result$estimated_quadrature_posterior_error > 0.01)) {
    stop("one or more replayed results violated JG-0.2.7 safeguards")
}
