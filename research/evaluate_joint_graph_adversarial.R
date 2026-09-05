#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (!length(script_argument)) stop("cannot locate evaluation script")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(root, "research", "joint_graph_reference.R"))

trailing <- commandArgs(trailingOnly = TRUE)
output <- if (length(trailing) >= 1L) trailing[[1]] else
    file.path(root, "research", "joint_graph_v0_1_adversarial_summary.tsv")
replicates_per_cell <- if (length(trailing) >= 2L) as.integer(trailing[[2]]) else 50L
requested_cores <- if (length(trailing) >= 3L) as.integer(trailing[[3]]) else
    min(8L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
if (!is.finite(replicates_per_cell) || replicates_per_cell < 1L) {
    stop("replicates must be a positive integer")
}
if (!is.finite(requested_cores) || requested_cores < 1L) {
    stop("cores must be a positive integer")
}

cpp_binary <- file.path(root, "build", "joint_graph_cli")
if (!file.exists(cpp_binary)) stop("build/joint_graph_cli is missing")

path_value <- sqrt(2) * 0.70 * 0.8162878828589647
residual_value <- sqrt(2) * 0.175 * 0.8162878828589647

scenario <- function(name, status, pattern = "diagnostic", a = 0, b = 0,
                     lambda = 0, q = 0, se = 0.02,
                     counts = c(A = 60L, B = 60L, C = 60L),
                     variance_scale = 1, ld_rho = 0, overlap_rho = 0,
                     directional_d = 0, outlier_rate = 0,
                     outlier_shift = 0, swap_role_fraction = 0) {
    list(
        name = name, status = status, pattern = pattern, a = a, b = b,
        c_path = residual_value, lambda = lambda, q = q, se = se,
        counts = counts, variance_scale = variance_scale, ld_rho = ld_rho,
        overlap_rho = overlap_rho, directional_d = directional_d,
        outlier_rate = outlier_rate, outlier_shift = outlier_shift,
        swap_role_fraction = swap_role_fraction
    )
}

scenarios <- list(
    scenario("matched_null", "matched", "null"),
    scenario("matched_mediation", "matched", "mediation",
             a = path_value, b = path_value),
    scenario("matched_pleiotropy", "matched", "pleiotropy",
             a = path_value, lambda = path_value, q = 0.35),
    scenario("matched_coexistence", "matched", "coexistence",
             a = path_value, b = path_value, lambda = path_value, q = 0.35),
    scenario("moderate_mediation", "numerical_stress", "mediation",
             a = 0.40, b = 0.40),
    scenario("weak_mediation", "weak_information",
             a = 0.20, b = 0.20, se = 0.05),
    scenario("few_B_mediation", "weak_information", "mediation",
             a = path_value, b = path_value,
             counts = c(A = 60L, B = 8L, C = 0L)),
    scenario("no_M_anchor_mediation", "identification_boundary",
             a = path_value, b = path_value,
             counts = c(A = 60L, B = 0L, C = 0L)),
    scenario("sparse_pleiotropy", "misspecified_mixture",
             a = path_value, lambda = path_value, q = 0.05),
    scenario("near_aligned_pleiotropy", "near_boundary",
             a = path_value, lambda = path_value, q = 0.95),
    scenario("exact_aligned_pleiotropy", "nonidentified_boundary",
             a = path_value, lambda = path_value, q = 1),
    scenario("directional_d_pleiotropy", "misspecified_mean",
             a = path_value, directional_d = path_value),
    scenario("sparse_outcome_outliers", "misspecified_tails",
             a = path_value, outlier_rate = 0.05, outlier_shift = 0.80),
    scenario("overlap_null", "undeclared_sample_overlap",
             overlap_rho = 0.50),
    scenario("ld_null", "ignored_ld", ld_rho = 0.80),
    scenario("ld_mediation", "ignored_ld",
             a = path_value, b = path_value, ld_rho = 0.80),
    scenario("low_variance_null", "scale_misspecification",
             variance_scale = 0.25),
    scenario("high_variance_null", "scale_misspecification",
             variance_scale = 4),
    scenario("low_variance_mediation", "scale_misspecification",
             a = path_value, b = path_value, variance_scale = 0.25),
    scenario("role_contamination_pleiotropy", "role_misspecification",
             a = path_value, lambda = path_value, q = 0.35,
             swap_role_fraction = 0.20)
)

ar1_factor <- function(n, rho) {
    if (rho == 0) return(diag(n))
    indices <- seq_len(n)
    covariance <- rho^abs(outer(indices, indices, "-"))
    t(chol(covariance))
}

simulate_adversarial <- function(specification, seed) {
    set.seed(seed)
    roles <- rep(names(specification$counts), specification$counts)
    n <- length(roles)
    if (!n) stop("scenario contains no variants")
    config <- jg_default_config()
    vx <- unname(config$vx[roles]) * specification$variance_scale
    vm <- unname(config$vm[roles]) * specification$variance_scale
    vy <- config$vy * specification$variance_scale
    ld_factor <- ar1_factor(n, specification$ld_rho)

    correlated_normal <- function(sd) {
        as.vector(ld_factor %*% rnorm(n)) * sd
    }
    g <- correlated_normal(1) * sqrt(vx)
    d <- correlated_normal(1) * sqrt(vm)
    direct_y <- correlated_normal(sqrt(vy))
    h <- if (specification$q > 0) rbinom(n, 1, specification$q) else rep(0, n)
    true_m <- specification$a * g + d
    true_y <- specification$c_path * g + specification$b * true_m +
        h * specification$lambda * d + direct_y

    if (specification$directional_d != 0) {
        true_y <- true_y + specification$directional_d * abs(d)
    }
    if (specification$outlier_rate > 0) {
        outlier <- rbinom(n, 1, specification$outlier_rate)
        true_y <- true_y + outlier * specification$outlier_shift
    }

    trait_correlation <- matrix(specification$overlap_rho, 3, 3)
    diag(trait_correlation) <- 1
    trait_factor <- t(chol(trait_correlation))
    independent_error <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    error <- ld_factor %*% independent_error %*% t(trait_factor)
    error <- error * specification$se

    reported_roles <- roles
    if (specification$swap_role_fraction > 0) {
        a_indices <- which(roles == "A")
        b_indices <- which(roles == "B")
        count <- floor(min(length(a_indices), length(b_indices)) *
                       specification$swap_role_fraction)
        if (count > 0) {
            swap_a <- sample(a_indices, count)
            swap_b <- sample(b_indices, count)
            reported_roles[swap_a] <- "B"
            reported_roles[swap_b] <- "A"
        }
    }

    data.frame(
        variant = sprintf("v%04d", seq_len(n)), role = reported_roles,
        beta_x = g + error[, 1], se_x = specification$se,
        beta_m = true_m + error[, 2], se_m = specification$se,
        beta_y = true_y + error[, 3], se_y = specification$se,
        stringsAsFactors = FALSE
    )
}

matches_pattern <- function(pattern, result) {
    if (pattern == "diagnostic") return(NA)
    if (pattern == "null") {
        return(result$PP_XM <= 0.20 && result$PP_global_MY <= 0.20)
    }
    if (pattern == "mediation") {
        return(result$PP_XM >= 0.80 && result$PP_global_MY >= 0.80 &&
               result$PP_nonaligned_P <= 0.20)
    }
    if (pattern == "pleiotropy") {
        return(result$PP_XM >= 0.80 && result$PP_global_MY <= 0.20 &&
               result$PP_nonaligned_P >= 0.80)
    }
    result$PP_XM >= 0.70 && result$PP_global_MY >= 0.70 &&
        result$PP_nonaligned_P >= 0.70
}

jobs <- do.call(rbind, lapply(seq_along(scenarios), function(index) {
    data.frame(
        scenario_index = index,
        replicate = seq_len(replicates_per_cell),
        seed = 20262000L + 1000L * index + seq_len(replicates_per_cell)
    )
}))

run_job <- function(job_index) {
    job <- jobs[job_index, ]
    specification <- scenarios[[job$scenario_index]]
    data <- simulate_adversarial(specification, job$seed)
    input <- tempfile("jg_adversarial_", fileext = ".tsv")
    result_path <- tempfile("jg_adversarial_result_", fileext = ".tsv")
    on.exit(unlink(c(input, result_path)), add = TRUE)
    write.table(data, input, sep = "\t", quote = FALSE, row.names = FALSE)
    timing <- system.time(status <- system2(cpp_binary, c(input, result_path)))
    if (!identical(status, 0L)) stop("C++ evaluator failed for ", specification$name)
    result <- read.delim(result_path, check.names = FALSE,
                         stringsAsFactors = FALSE)
    if (!identical(
        result$identification_scope,
        "CONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY"
    )) stop("identification scope was not preserved")
    data.frame(
        scenario = specification$name,
        status = specification$status,
        pattern = specification$pattern,
        replicate = job$replicate,
        seed = job$seed,
        n_variants = nrow(data),
        true_a = specification$a,
        true_b = specification$b,
        true_lambda = specification$lambda,
        true_q = specification$q,
        PP_XM = result$PP_XM,
        PP_global_MY = result$PP_global_MY,
        PP_nonaligned_P = result$PP_nonaligned_P,
        PP_two_path = result$PP_two_path,
        PP_two_path_plus_P = result$PP_two_path_plus_P,
        pattern_correct = matches_pattern(specification$pattern, result),
        high_two_path = result$PP_two_path >= 0.80,
        elapsed_seconds = timing[["elapsed"]],
        stringsAsFactors = FALSE
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
    evaluated <- !is.na(x$pattern_correct)
    data.frame(
        scenario = specification$name,
        status = specification$status,
        pattern = specification$pattern,
        replicates = nrow(x),
        correct = if (any(evaluated)) sum(x$pattern_correct[evaluated]) else NA,
        correct_rate = if (any(evaluated)) mean(x$pattern_correct[evaluated]) else NA,
        high_two_path_rate = mean(x$high_two_path),
        mean_PP_XM = mean(x$PP_XM),
        mean_PP_global_MY = mean(x$PP_global_MY),
        mean_PP_nonaligned_P = mean(x$PP_nonaligned_P),
        mean_PP_two_path = mean(x$PP_two_path),
        median_elapsed_seconds = median(x$elapsed_seconds),
        stringsAsFactors = FALSE
    )
}))
row.names(summaries) <- NULL

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
replicate_output <- sub("\\.tsv$", "_replicates.tsv", output)
write.table(summaries, output, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "NA")
write.table(replicate_results, replicate_output, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
print(summaries, row.names = FALSE)

if (replicates_per_cell >= 50L) {
    tested <- !is.na(summaries$correct_rate)
    failed <- summaries$scenario[tested & summaries$correct_rate < 0.80]
    if (length(failed)) {
        message("Prespecified pattern failures: ", paste(failed, collapse = ", "))
    } else {
        message("All prespecified pattern cells passed.")
    }
} else {
    message("Smoke run only: acceptance requires at least 50 replicates per cell.")
}
