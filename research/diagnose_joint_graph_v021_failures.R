#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (!length(script_argument)) stop("cannot locate diagnostic script")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
    stop(paste(
        "usage: diagnose_joint_graph_v021_failures.R",
        "BINARY RESULT_ROOT OUTPUT.tsv [SEED_BASE]"
    ))
}
binary <- normalizePath(args[[1]])
result_root <- normalizePath(args[[2]])
output <- args[[3]]
seed_base <- if (length(args) == 4L) as.integer(args[[4]]) else 20400000L
if (!is.finite(seed_base) || seed_base < 1L) stop("SEED_BASE is invalid")
root_candidates <- c(dirname(binary), file.path(dirname(binary), ".."))
support_paths <- file.path(
    root_candidates, "research", "joint_graph_v02_simulation.R"
)
support_path <- support_paths[file.exists(support_paths)][1]
if (is.na(support_path)) {
    stop("cannot locate research/joint_graph_v02_simulation.R from BINARY")
}
source(normalizePath(support_path))

paths <- list.files(result_root, pattern = "proteins\\.tsv$", recursive = TRUE,
                    full.names = TRUE)
if (length(paths) != 500L) stop("expected 500 protein result files, found ", length(paths))
results <- do.call(rbind, lapply(paths, read.delim, stringsAsFactors = FALSE,
                                check.names = FALSE))
failed <- results[!results$success, c("scenario", "replicate", "protein", "truth")]

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

temporary <- tempfile("jg021_failure_replay_")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
diagnostics <- vector("list", nrow(failed))
for (i in seq_len(nrow(failed))) {
    x <- failed[i, ]
    scenario_index <- match(x$scenario, scenario_names)
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
    status <- system2(binary, c(input, ld_path, result_path), stdout = FALSE,
                      stderr = stderr_path)
    replay_succeeded <- identical(status, 0L) && file.exists(result_path)
    message <- if (file.exists(stderr_path)) {
        paste(readLines(stderr_path, warn = FALSE), collapse = " ")
    } else {
        ""
    }
    relaxed_status <- NA_integer_
    max_relevant_quadrature_difference <- NA_real_
    posterior_error <- NA_real_
    pp_two_path <- NA_real_
    if (replay_succeeded) {
        replay_result <- read.delim(result_path, check.names = FALSE)
        max_relevant_quadrature_difference <- if (
            "max_relevant_quadrature_difference" %in% names(replay_result)
        ) replay_result$max_relevant_quadrature_difference else
            replay_result$max_relevant_evidence_difference
        if ("estimated_quadrature_posterior_error" %in% names(replay_result)) {
            posterior_error <- replay_result$estimated_quadrature_posterior_error
        }
        pp_two_path <- replay_result$PP_two_path
    }
    numerical_match <- regexec(
        paste0("posterior_error=([^,]+), ",
               "max_relevant_quadrature_difference=([^,]+),"),
        message
    )
    numerical_values <- regmatches(message, numerical_match)[[1]]
    if (length(numerical_values) == 3L) {
        posterior_error <- as.numeric(numerical_values[[2]])
        max_relevant_quadrature_difference <-
            as.numeric(numerical_values[[3]])
    }
    if (!identical(status, 0L) &&
        !is.finite(posterior_error) &&
        grepl("adaptive (evidence diagnostic|quadrature did not converge)", message)) {
        options_path <- file.path(temporary, "diagnostic_options.tsv")
        writeLines(c("max_evidence_discrepancy\t100",
                     "max_quadrature_posterior_error\t100"), options_path)
        unlink(result_path)
        relaxed_status <- system2(
            binary, c(input, ld_path, result_path, options_path),
            stdout = FALSE, stderr = FALSE
        )
        if (identical(relaxed_status, 0L) && file.exists(result_path)) {
            relaxed <- read.delim(result_path, check.names = FALSE)
            max_relevant_quadrature_difference <-
                relaxed$max_relevant_quadrature_difference
            posterior_error <- relaxed$estimated_quadrature_posterior_error
            pp_two_path <- relaxed$PP_two_path
        }
    }
    diagnostics[[i]] <- data.frame(
        scenario = x$scenario,
        replicate = x$replicate,
        protein = x$protein,
        truth = x$truth,
        replay_status = as.integer(status),
        replay_succeeded = replay_succeeded,
        relaxed_status = relaxed_status,
        max_relevant_quadrature_difference =
            max_relevant_quadrature_difference,
        estimated_quadrature_posterior_error = posterior_error,
        relaxed_PP_two_path = pp_two_path,
        diagnostic = message,
        stringsAsFactors = FALSE
    )
}

diagnostics <- if (length(diagnostics)) do.call(rbind, diagnostics) else failed
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(diagnostics, output, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Replayed", nrow(failed), "failed protein analyses\n")
