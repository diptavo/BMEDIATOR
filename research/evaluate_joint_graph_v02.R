#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (!length(script_argument)) stop("cannot locate evaluation script")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(root, "research", "joint_graph_v02_simulation.R"))

trailing <- commandArgs(trailingOnly = TRUE)
output <- if (length(trailing) >= 1L) trailing[[1]] else
    file.path(root, "research", "joint_graph_v0_2_development_summary.tsv")
replicates_per_cell <- if (length(trailing) >= 2L) as.integer(trailing[[2]]) else 50L
requested_cores <- if (length(trailing) >= 3L) as.integer(trailing[[3]]) else
    min(8L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
if (!is.finite(replicates_per_cell) || replicates_per_cell < 1L) {
    stop("replicates must be positive")
}

binary <- file.path(root, "build", "joint_graph_v02_cli")
if (!file.exists(binary)) stop("build/joint_graph_v02_cli is missing")

scenario <- function(name, pattern = "diagnostic", profile = "default",
                     seed_group = NA_integer_, ...) {
    list(name = name, pattern = pattern, profile = profile,
         seed_group = seed_group, arguments = list(...))
}
scenarios <- list(
    scenario("matched_null", "null"),
    scenario("moderate_mediation", "mediation", a = 0.40, b = 0.40),
    scenario("weak_mediation", a = 0.20, b = 0.20, se = 0.05),
    scenario("sparse_pleiotropy", "sparse", a = 0.60,
             lambda = 0.70, q = 0.35),
    scenario("directional_pleiotropy", "directional", a = 0.60, eta = 0.40),
    scenario("two_pleiotropy", "two_pleiotropy", a = 0.60,
             lambda = 0.70, q = 0.35, eta = 0.40),
    scenario("mediation_sparse", "mediation_sparse", a = 0.40, b = 0.40,
             lambda = 0.70, q = 0.35),
    scenario("declared_overlap_null", "null",
             sampling_rho = c(xm = 0.40, xy = 0.30, my = 0.20)),
    scenario("undeclared_overlap_null",
             sampling_rho = c(xm = 0.40, xy = 0.30, my = 0.20),
             reported_sampling_rho = c(xm = 0, xy = 0, my = 0)),
    scenario("ld_null", "null", ld_rho = 0.70),
    scenario("ld_mediation", "mediation", ld_rho = 0.70, a = 0.40, b = 0.40),
    scenario("high_scale_declared_null", "null", true_variance_scale = 4),
    scenario("high_scale_misspecified_null", true_variance_scale = 4,
             reported_variance_scale = 1),
    scenario("low_scale_mediation", true_variance_scale = 0.25,
             a = 0.40, b = 0.40),
    scenario("few_blocks_sparse", blocks = 6L, variants_per_block = 3L,
             a = 0.60, lambda = 0.70, q = 0.35),
    scenario("near_aligned_pleiotropy", a = 0.60, lambda = 0.70, q = 0.90),
    scenario("exact_aligned_pleiotropy", a = 0.60, lambda = 0.70, q = 1),
    scenario("orientation_70pct", a = 0.60, eta = 0.40,
             orientation_accuracy = 0.70),
    scenario("skeptical_prior_mediation", "mediation", "skeptical", 2L,
             a = 0.40, b = 0.40),
    scenario("diffuse_prior_null", "null", "diffuse", 1L)
)

profiles <- list(
    skeptical = c(
        pi_xm = 0.10, pi_my = 0.05, pi_sparse = 0.05,
        pi_directional = 0.05, prior_sd_a = 0.50, prior_sd_b = 0.50,
        prior_sd_lambda = 0.50, prior_sd_eta = 0.50
    ),
    diffuse = c(
        pi_xm = 0.50, pi_my = 0.30, pi_sparse = 0.20,
        pi_directional = 0.20, prior_sd_a = 1.0, prior_sd_b = 1.0,
        prior_sd_lambda = 1.0, prior_sd_eta = 1.0
    )
)

temporary <- tempfile("jg02_validation_")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
profile_paths <- list()
for (name in names(profiles)) {
    path <- file.path(temporary, paste0(name, ".options.tsv"))
    write.table(data.frame(key = names(profiles[[name]]), value = profiles[[name]]),
                path, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = FALSE)
    profile_paths[[name]] <- path
}

matches_pattern <- function(pattern, result) {
    if (pattern == "diagnostic") return(NA)
    if (pattern == "null") {
        return(result$PP_two_path <= 0.20 && result$PP_global_MY <= 0.20)
    }
    if (pattern == "mediation") {
        return(result$PP_two_path >= 0.80 && result$PP_any_P <= 0.20)
    }
    if (pattern == "sparse") {
        return(result$PP_sparse_P >= 0.80 && result$PP_global_MY <= 0.20)
    }
    if (pattern == "directional") {
        return(result$PP_directional_P >= 0.80 && result$PP_global_MY <= 0.20)
    }
    if (pattern == "two_pleiotropy") {
        return(result$PP_sparse_P >= 0.70 && result$PP_directional_P >= 0.70 &&
               result$PP_global_MY <= 0.30)
    }
    result$PP_two_path >= 0.70 && result$PP_sparse_P >= 0.70
}

jobs <- do.call(rbind, lapply(seq_along(scenarios), function(index) {
    seed_group <- if (is.na(scenarios[[index]]$seed_group)) index else
        scenarios[[index]]$seed_group
    data.frame(scenario_index = index, replicate = seq_len(replicates_per_cell),
               seed = 20264000L + 1000L * seed_group + seq_len(replicates_per_cell))
}))

run_job <- function(job_index) {
    job <- jobs[job_index, ]
    specification <- scenarios[[job$scenario_index]]
    defaults <- list(blocks = 20L, variants_per_block = 2L, ld_rho = 0.40)
    simulation_arguments <- modifyList(defaults, specification$arguments)
    simulation_arguments$seed <- job$seed
    fixture <- do.call(jg02_simulate, simulation_arguments)
    input <- tempfile("jg02_", tmpdir = temporary, fileext = ".tsv")
    ld_path <- tempfile("jg02_ld_", tmpdir = temporary, fileext = ".tsv")
    result_path <- tempfile("jg02_result_", tmpdir = temporary, fileext = ".tsv")
    on.exit(unlink(c(input, ld_path, result_path)), add = TRUE)
    jg02_write_fixture(fixture, input, ld_path)
    command_arguments <- c(input, ld_path, result_path)
    if (specification$profile != "default") {
        command_arguments <- c(command_arguments,
                               profile_paths[[specification$profile]])
    }
    timing <- system.time(status <- system2(binary, command_arguments,
                                            stdout = FALSE, stderr = FALSE))
    if (!identical(status, 0L) || !file.exists(result_path)) {
        return(data.frame(
            scenario = specification$name, pattern = specification$pattern,
            profile = specification$profile, replicate = job$replicate,
            seed = job$seed, success = FALSE, pattern_correct = FALSE,
            PP_XM = NA, PP_global_MY = NA, PP_sparse_P = NA,
            PP_directional_P = NA, PP_any_P = NA, PP_two_path = NA,
            max_adaptive_laplace_difference = NA,
            elapsed_seconds = timing[["elapsed"]], stringsAsFactors = FALSE
        ))
    }
    result <- read.delim(result_path, check.names = FALSE,
                         stringsAsFactors = FALSE)
    data.frame(
        scenario = specification$name, pattern = specification$pattern,
        profile = specification$profile, replicate = job$replicate,
        seed = job$seed, success = TRUE,
        pattern_correct = matches_pattern(specification$pattern, result),
        PP_XM = result$PP_XM, PP_global_MY = result$PP_global_MY,
        PP_sparse_P = result$PP_sparse_P,
        PP_directional_P = result$PP_directional_P,
        PP_any_P = result$PP_any_P, PP_two_path = result$PP_two_path,
        max_adaptive_laplace_difference = result$max_adaptive_laplace_difference,
        elapsed_seconds = timing[["elapsed"]], stringsAsFactors = FALSE
    )
}

cores <- min(requested_cores, nrow(jobs))
if (.Platform$OS.type == "unix" && cores > 1L) {
    rows <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = cores,
                              mc.preschedule = FALSE)
} else {
    rows <- lapply(seq_len(nrow(jobs)), run_job)
}
replicate_results <- do.call(rbind, rows)

summaries <- do.call(rbind, lapply(scenarios, function(specification) {
    x <- replicate_results[replicate_results$scenario == specification$name, ]
    evaluated <- x$success & !is.na(x$pattern_correct)
    valid <- x$success
    data.frame(
        scenario = specification$name,
        pattern = specification$pattern,
        profile = specification$profile,
        replicates = nrow(x),
        successful = sum(valid),
        correct = if (any(evaluated)) sum(x$pattern_correct[evaluated]) else NA,
        correct_rate = if (any(evaluated)) mean(x$pattern_correct[evaluated]) else NA,
        high_two_path_rate = mean(x$PP_two_path[valid] >= 0.80),
        mean_PP_XM = mean(x$PP_XM[valid]),
        mean_PP_global_MY = mean(x$PP_global_MY[valid]),
        mean_PP_sparse_P = mean(x$PP_sparse_P[valid]),
        mean_PP_directional_P = mean(x$PP_directional_P[valid]),
        mean_PP_two_path = mean(x$PP_two_path[valid]),
        max_evidence_difference = max(x$max_adaptive_laplace_difference[valid]),
        median_elapsed_seconds = median(x$elapsed_seconds),
        stringsAsFactors = FALSE
    )
}))
row.names(summaries) <- NULL
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(summaries, output, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "NA")
write.table(replicate_results, sub("\\.tsv$", "_replicates.tsv", output),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
print(summaries, row.names = FALSE)

if (replicates_per_cell >= 50L) {
    required <- !is.na(summaries$correct_rate)
    failed <- summaries$scenario[
        summaries$successful < replicates_per_cell |
        (required & summaries$correct_rate < 0.80)
    ]
    if (length(failed)) {
        message("JG-0.2 prespecified failures: ", paste(failed, collapse = ", "))
    } else {
        message("All JG-0.2 prespecified development cells passed.")
    }
} else {
    message("Smoke run only; full acceptance requires 50 replicates per cell.")
}
