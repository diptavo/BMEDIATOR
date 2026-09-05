#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (!length(script_argument)) stop("cannot locate task script")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(root, "research", "joint_graph_v02_simulation.R"))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(5L, 6L)) {
    stop(paste(
        "usage: run_joint_graph_v021_family_task.R",
        "BINARY SCENARIO REPLICATE OUTDIR N_PROTEINS [SEED_BASE]"
    ))
}
binary <- normalizePath(args[[1]])
scenario_name <- args[[2]]
replicate <- as.integer(args[[3]])
outdir <- args[[4]]
n_proteins <- as.integer(args[[5]])
seed_base <- if (length(args) == 6L) as.integer(args[[6]]) else 20400000L
if (!is.finite(replicate) || replicate < 1L ||
    !is.finite(n_proteins) || n_proteins < 1L || n_proteins > 100L ||
    !is.finite(seed_base) || seed_base < 1L) {
    stop("replicate, N_PROTEINS, and SEED_BASE are invalid")
}

scenario_names <- c(
    "baseline", "rare", "composite_null", "mixed", "strong_ld",
    "ld_mismatch", "scale_uncertainty", "undeclared_overlap",
    "weak_paths", "aligned_sensitivity"
)
scenario_index <- match(scenario_name, scenario_names)
if (is.na(scenario_index)) stop("unknown family scenario: ", scenario_name)

composition <- switch(
    scenario_name,
    baseline = c(mediation = 10, null = 90),
    rare = c(mediation = 2, null = 98),
    composite_null = c(null = 20, xm_only = 20, my_only = 20,
                       sparse = 20, directional = 20),
    mixed = c(mediation = 10, mediation_sparse = 10, sparse = 20,
              directional = 20, null = 40),
    strong_ld = c(mediation = 10, null = 90),
    ld_mismatch = c(mediation = 10, null = 90),
    scale_uncertainty = c(mediation = 10, null = 90),
    undeclared_overlap = c(mediation = 10, null = 90),
    weak_paths = c(weak_mediation = 10, null = 90),
    aligned_sensitivity = c(aligned = 20, null = 80)
)
types <- rep(names(composition), composition)
family_seed <- seed_base + scenario_index * 100000L + replicate * 1000L
set.seed(family_seed)
types <- sample(types, length(types), replace = FALSE)
types <- types[seq_len(n_proteins)]

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
        stop("unknown protein type")
    )
}

replicate_directory <- file.path(
    outdir, scenario_name, sprintf("rep_%04d", replicate)
)
dir.create(replicate_directory, recursive = TRUE, showWarnings = FALSE)
temporary <- tempfile("jg021_family_")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)

rows <- vector("list", n_proteins)
for (protein in seq_len(n_proteins)) {
    type <- types[[protein]]
    simulation_arguments <- modifyList(
        list(seed = family_seed + protein, blocks = 20L,
             variants_per_block = 2L, ld_rho = 0.40),
        type_arguments(type)
    )
    if (scenario_name %in% c("strong_ld", "ld_mismatch")) {
        simulation_arguments$ld_rho <- 0.70
    }
    if (scenario_name == "scale_uncertainty") {
        set.seed(family_seed + 500000L + protein)
        simulation_arguments$reported_variance_scale <- exp(rnorm(1, sd = 0.25))
    }
    if (scenario_name == "undeclared_overlap") {
        simulation_arguments$sampling_rho <- c(xm = 0.40, xy = 0.30, my = 0.20)
        simulation_arguments$reported_sampling_rho <- c(xm = 0, xy = 0, my = 0)
    }
    fixture <- do.call(jg02_simulate, simulation_arguments)
    if (scenario_name == "ld_mismatch") {
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
    } else {
        ""
    }
    result_columns <- c(
        "PP_XM", "PP_global_MY", "PP_sparse_P", "PP_directional_P",
        "PP_any_P", "PP_two_path", "max_relevant_evidence_difference",
        "max_relevant_quadrature_difference",
        "estimated_quadrature_posterior_error", "max_quadrature_order",
        "posterior_aware_refinements", "states_regularized"
    )
    if (success) {
        result <- read.delim(result_path, check.names = FALSE,
                             stringsAsFactors = FALSE)
        if (!identical(result$model_version, "JG-0.2.6")) {
            stop("unexpected joint model version: ", result$model_version)
        }
        values <- result[result_columns]
    } else {
        values <- as.data.frame(as.list(setNames(
            rep(NA_real_, length(result_columns)), result_columns
        )))
    }
    rows[[protein]] <- data.frame(
        scenario = scenario_name,
        replicate = replicate,
        protein = sprintf("P%03d", protein),
        truth = type,
        true_mediator = type %in% c("mediation", "mediation_sparse", "weak_mediation"),
        success = success,
        exit_status = as.integer(status),
        diagnostic = diagnostic,
        values,
        elapsed_seconds = timing[["elapsed"]],
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}
protein_results <- do.call(rbind, rows)

bayesian_fdr_selection <- function(pp, alpha = 0.05) {
    if (any(!is.finite(pp))) return(rep(NA, length(pp)))
    selected <- rep(FALSE, length(pp))
    ordered <- order(1 - pp)
    acceptable <- which(cumsum(1 - pp[ordered]) / seq_along(ordered) <= alpha)
    if (length(acceptable)) selected[ordered[seq_len(max(acceptable))]] <- TRUE
    selected
}

protein_results$selected_bfdr05 <- bayesian_fdr_selection(protein_results$PP_two_path)
protein_results$selected_pp80 <- if (all(protein_results$success)) {
    protein_results$PP_two_path >= 0.80
} else {
    rep(NA, nrow(protein_results))
}

metric <- function(selected) {
    if (anyNA(selected)) {
        return(c(discoveries = NA, true_discoveries = NA,
                 false_discoveries = NA, fdp = NA, power = NA))
    }
    discoveries <- sum(selected)
    true_discoveries <- sum(selected & protein_results$true_mediator)
    false_discoveries <- discoveries - true_discoveries
    total_true <- sum(protein_results$true_mediator)
    c(
        discoveries = discoveries,
        true_discoveries = true_discoveries,
        false_discoveries = false_discoveries,
        fdp = if (discoveries > 0) false_discoveries / discoveries else 0,
        power = if (total_true > 0) true_discoveries / total_true else NA
    )
}
bfdr <- metric(protein_results$selected_bfdr05)
pp80 <- metric(protein_results$selected_pp80)
summary <- data.frame(
    model_version = "JG-0.2.6",
    scenario = scenario_name,
    replicate = replicate,
    family_seed = family_seed,
    proteins = n_proteins,
    true_mediators = sum(protein_results$true_mediator),
    successful = sum(protein_results$success),
    family_complete = all(protein_results$success),
    bfdr05_discoveries = bfdr[["discoveries"]],
    bfdr05_true_discoveries = bfdr[["true_discoveries"]],
    bfdr05_false_discoveries = bfdr[["false_discoveries"]],
    bfdr05_fdp = bfdr[["fdp"]],
    bfdr05_power = bfdr[["power"]],
    pp80_discoveries = pp80[["discoveries"]],
    pp80_true_discoveries = pp80[["true_discoveries"]],
    pp80_false_discoveries = pp80[["false_discoveries"]],
    pp80_fdp = pp80[["fdp"]],
    pp80_power = pp80[["power"]],
    median_seconds = median(protein_results$elapsed_seconds),
    stringsAsFactors = FALSE
)

write.table(protein_results, file.path(replicate_directory, "proteins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(summary, file.path(replicate_directory, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Completed", scenario_name, "replicate", replicate, "with",
    summary$successful, "successful proteins\n")
