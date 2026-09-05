#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (!length(script_argument)) stop("cannot locate test script")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(root, "research", "joint_graph_v02_simulation.R"))
scale_environment <- new.env(parent = globalenv())
sys.source(file.path(root, "research", "estimate_joint_graph_v02_scales.R"),
           envir = scale_environment)

binary <- file.path(root, "build", "joint_graph_v02_cli")
if (!file.exists(binary)) stop("build/joint_graph_v02_cli is missing")
output_directory <- file.path(root, "build", "test", "joint_graph_v02")
unlink(output_directory, recursive = TRUE)
dir.create(output_directory, recursive = TRUE)

scenarios <- list(
    null = list(seed = 20263001),
    moderate_mediation = list(seed = 20263002, a = 0.40, b = 0.40),
    sparse_pleiotropy = list(seed = 20263003, a = 0.60,
                             lambda = 0.70, q = 0.35),
    directional_pleiotropy = list(seed = 20263004, a = 0.60, eta = 0.40),
    uncertain_directional = list(seed = 20263009, a = 0.60, eta = 0.40,
                                 orientation_accuracy = 0.70),
    overlap_null = list(seed = 20263005,
                        sampling_rho = c(xm = 0.40, xy = 0.30, my = 0.20)),
    ld_mediation = list(seed = 20263006, ld_rho = 0.70,
                        a = 0.40, b = 0.40),
    curved_strong_ld = list(seed = 20915038, ld_rho = 0.70,
                            a = 0.40, b = 0.40),
    small_weight_curvature_null = list(seed = 41401056)
)

rows <- vector("list", length(scenarios))
for (index in seq_along(scenarios)) {
    name <- names(scenarios)[index]
    arguments <- c(list(blocks = 20L, variants_per_block = 2L),
                   scenarios[[index]])
    fixture <- do.call(jg02_simulate, arguments)
    input <- file.path(output_directory, paste0(name, ".tsv"))
    ld_path <- file.path(output_directory, paste0(name, ".ld.tsv"))
    output <- file.path(output_directory, paste0(name, ".result.tsv"))
    jg02_write_fixture(fixture, input, ld_path)
    fixed <- c(a = 0.31, b = -0.22, c_path = 0.08,
               lambda = 0.27, q = 0.35, eta = -0.11)
    r_loglik <- do.call(jg02_loglik, c(list(data = fixture$data, ld = fixture$ld),
                                       as.list(fixed)))
    cpp_loglik <- as.numeric(system2(
        binary,
        c("--loglik", input, ld_path, format(fixed, digits = 17)),
        stdout = TRUE
    ))
    if (abs(r_loglik - cpp_loglik) > 1e-8) {
        stop(name, " fixed likelihood R/C++ difference is ",
             abs(r_loglik - cpp_loglik))
    }
    if (name == "sparse_pleiotropy") {
        fixed_without_q <- fixed[names(fixed) != "q"]
        q_grid <- seq(0.01, 0.99, length.out = 99)
        grid_loglik <- vapply(q_grid, function(q) {
            do.call(jg02_loglik, c(
                list(data = fixture$data, ld = fixture$ld),
                as.list(fixed_without_q), list(q = q)
            ))
        }, numeric(1))
        center <- max(grid_loglik)
        integrand <- Vectorize(function(q) {
            exp(do.call(jg02_loglik, c(
                list(data = fixture$data, ld = fixture$ld),
                as.list(fixed_without_q), list(q = q)
            )) - center) * dbeta(q, 2, 2)
        })
        r_integrated <- log(integrate(
            integrand, 0, 1, rel.tol = 1e-9, subdivisions = 500L
        )$value) + center
        cpp_integrated <- as.numeric(system2(
            binary,
            c("--loglik-integrated-q", input, ld_path,
              format(fixed_without_q, digits = 17)),
            stdout = TRUE
        ))
        if (abs(r_integrated - cpp_integrated) > 1e-7) {
            stop("exact q marginalization R/C++ difference is ",
                 abs(r_integrated - cpp_integrated))
        }
    }
    timing <- system.time(status <- system2(binary, c(input, ld_path, output)))
    if (!identical(status, 0L)) stop("JG-0.2 failed for ", name)
    result <- read.delim(output, check.names = FALSE,
                         stringsAsFactors = FALSE)
    if (!identical(result$model_version, "JG-0.2.4")) stop("wrong model version")
    if (result$estimated_quadrature_posterior_error > 0.01 + 1e-10 ||
        result$max_relevant_quadrature_difference > 1 + 1e-10) {
        stop("successfully reported posterior failed quadrature stability")
    }
    if (!identical(
        result$identification_scope,
        "CONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY"
    )) stop("identification scope was lost")
    if (result$states_converged < 16L) stop(name, " optimizer did not converge")
    rows[[index]] <- data.frame(
        scenario = name,
        PP_XM = result$PP_XM,
        PP_global_MY = result$PP_global_MY,
        PP_sparse_P = result$PP_sparse_P,
        PP_directional_P = result$PP_directional_P,
        PP_two_path = result$PP_two_path,
        states_regularized = result$states_regularized,
        estimated_quadrature_posterior_error =
            result$estimated_quadrature_posterior_error,
        elapsed_seconds = timing[["elapsed"]],
        stringsAsFactors = FALSE
    )
}

results <- do.call(rbind, rows)
print(results, row.names = FALSE)
write.table(results, file.path(output_directory, "results.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

row_for <- function(name) results[results$scenario == name, , drop = FALSE]
failures <- character()
require_result <- function(condition, message) {
    if (!isTRUE(condition)) failures <<- c(failures, message)
}
null <- row_for("null")
mediation <- row_for("moderate_mediation")
sparse <- row_for("sparse_pleiotropy")
directional <- row_for("directional_pleiotropy")
uncertain_directional <- row_for("uncertain_directional")
overlap <- row_for("overlap_null")
ld_mediation <- row_for("ld_mediation")
curved_strong_ld <- row_for("curved_strong_ld")
small_weight_curvature_null <- row_for("small_weight_curvature_null")
require_result(null$PP_two_path < 0.20, "null produced two-path support")
require_result(mediation$PP_two_path > 0.80,
               "off-grid moderate mediation was not recovered")
require_result(mediation$PP_sparse_P < 0.20,
               "off-grid mediation was aliased with sparse pleiotropy")
require_result(sparse$PP_sparse_P > 0.70 && sparse$PP_global_MY < 0.30,
               "sparse pleiotropy was not separated from the global path")
require_result(directional$PP_directional_P > 0.70 &&
               directional$PP_global_MY < 0.30,
               "directional pleiotropy was not separated from the global path")
require_result(uncertain_directional$PP_directional_P > 0.70 &&
               uncertain_directional$PP_global_MY < 0.30,
               "uncertain orientation was not integrated correctly")
require_result(overlap$PP_two_path < 0.20,
               "declared overlap null produced two-path support")
require_result(ld_mediation$PP_two_path > 0.80,
               "signed-LD mediation was not recovered")
require_result(curved_strong_ld$PP_two_path > 0.80,
               "curved strong-LD mediation was not reportable")
require_result(small_weight_curvature_null$PP_two_path < 0.20,
               "small-weight curved null produced two-path support")

scale_fixture <- jg02_simulate(
    seed = 20263007,
    blocks = 300L,
    variants_per_block = 1L,
    ld_rho = 0,
    a = 0.40,
    b = 0.30,
    c_path = 0.10,
    sampling_rho = c(xm = 0.20, xy = 0.10, my = 0.15)
)
scales <- scale_environment$jg02_estimate_scales(scale_fixture$data)
truth <- jg02_default_variances(scales$role)
require_result(all(abs(log(scales$v_x / truth$v_x)) < log(3)),
               "external v_x estimates exceeded the factor-three tolerance")
require_result(all(abs(log(scales$v_m / truth$v_m)) < log(3)),
               "external v_m estimates exceeded the factor-three tolerance")
require_result(all(abs(log(scales$v_y / truth$v_y)) < log(3)),
               "external v_y estimate exceeded the factor-three tolerance")

bad_blocks <- jg02_simulate(seed = 20263008, blocks = 1L,
                            variants_per_block = 4L, ld_rho = 0.70)
bad_blocks$data$ld_block <- paste0("wrong", seq_len(nrow(bad_blocks$data)))
bad_input <- file.path(output_directory, "bad_blocks.tsv")
bad_ld <- file.path(output_directory, "bad_blocks.ld.tsv")
bad_output <- file.path(output_directory, "bad_blocks.result.tsv")
jg02_write_fixture(bad_blocks, bad_input, bad_ld)
bad_status <- suppressWarnings(system2(binary, c(bad_input, bad_ld, bad_output),
                                       stdout = FALSE, stderr = FALSE))
require_result(!identical(bad_status, 0L),
               "cross-block LD violation did not fail closed")

underidentified <- jg02_simulate(seed = 20263009, blocks = 4L,
                                 variants_per_block = 3L)
underidentified$data$role <- "C"
underidentified$data$role[underidentified$data$ld_block == "block001"] <- "A"
underidentified$data$role[underidentified$data$ld_block == "block002"] <- "B"
under_input <- file.path(output_directory, "underidentified.tsv")
under_ld <- file.path(output_directory, "underidentified.ld.tsv")
under_output <- file.path(output_directory, "underidentified.result.tsv")
jg02_write_fixture(underidentified, under_input, under_ld)
under_status <- suppressWarnings(system2(
    binary, c(under_input, under_ld, under_output), stdout = FALSE, stderr = FALSE
))
require_result(!identical(under_status, 0L),
               "single-block A/B evidence did not fail closed")

singular <- jg02_simulate(seed = 20263010, blocks = 6L,
                          variants_per_block = 2L)
singular$ld[1, 2] <- singular$ld[2, 1] <- 1
singular_input <- file.path(output_directory, "singular_ld.tsv")
singular_ld <- file.path(output_directory, "singular_ld.ld.tsv")
singular_output <- file.path(output_directory, "singular_ld.result.tsv")
jg02_write_fixture(singular, singular_input, singular_ld)
singular_status <- suppressWarnings(system2(
    binary, c(singular_input, singular_ld, singular_output),
    stdout = FALSE, stderr = FALSE
))
require_result(!identical(singular_status, 0L),
               "singular within-block LD did not fail closed")

duplicate_options <- file.path(output_directory, "duplicate_options.tsv")
writeLines(c("pi_xm\t0.20", "pi_xm\t0.30"), duplicate_options)
duplicate_option_status <- suppressWarnings(system2(
    binary, c(file.path(output_directory, "null.tsv"),
              file.path(output_directory, "null.ld.tsv"),
              file.path(output_directory, "duplicate_options.result.tsv"),
              duplicate_options), stdout = FALSE, stderr = FALSE
))
require_result(!identical(duplicate_option_status, 0L),
               "duplicate option was silently accepted")

unsafe_options <- file.path(output_directory, "unsafe_options.tsv")
writeLines(c("max_quadrature_posterior_error\t100", "min_role_blocks\t2"),
           unsafe_options)
unsafe_option_status <- suppressWarnings(system2(
    binary,
    c(file.path(output_directory, "null.tsv"),
      file.path(output_directory, "null.ld.tsv"),
      file.path(output_directory, "unsafe_options.result.tsv"), unsafe_options),
    stdout = FALSE, stderr = FALSE
))
require_result(!identical(unsafe_option_status, 0L),
               "release safeguards were relaxed through the options file")

strict_options <- file.path(output_directory, "strict_options.tsv")
strict_stderr <- file.path(output_directory, "strict_options.stderr.txt")
writeLines("max_quadrature_posterior_error\t1e-12", strict_options)
strict_option_status <- suppressWarnings(system2(
    binary,
    c(file.path(output_directory, "null.tsv"),
      file.path(output_directory, "null.ld.tsv"),
      file.path(output_directory, "strict_options.result.tsv"), strict_options),
    stdout = FALSE, stderr = strict_stderr
))
strict_message <- paste(readLines(strict_stderr, warn = FALSE), collapse = " ")
require_result(!identical(strict_option_status, 0L) &&
               grepl("posterior_error=", strict_message, fixed = TRUE) &&
               grepl("max_relevant_quadrature_difference=", strict_message,
                     fixed = TRUE),
               "suppressed fit did not report actionable numerical diagnostics")

input_lines <- readLines(file.path(output_directory, "null.tsv"))
input_fields <- strsplit(input_lines, "\t", fixed = TRUE)
duplicate_input <- vapply(input_fields, function(fields) {
    paste(c(fields[[1]], fields), collapse = "\t")
}, character(1))
duplicate_input_path <- file.path(output_directory, "duplicate_input_column.tsv")
writeLines(duplicate_input, duplicate_input_path)
duplicate_input_status <- suppressWarnings(system2(
    binary,
    c(duplicate_input_path, file.path(output_directory, "null.ld.tsv"),
      file.path(output_directory, "duplicate_input_column.result.tsv")),
    stdout = FALSE, stderr = FALSE
))
require_result(!identical(duplicate_input_status, 0L),
               "duplicate input column was silently accepted")
if (length(failures)) {
    stop("JG-0.2 acceptance failed:\n- ", paste(failures, collapse = "\n- "))
}
cat("JG-0.2.4 adaptive, LD, overlap, scale, and directional tests passed.\n")
