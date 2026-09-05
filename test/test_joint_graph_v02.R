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
                        a = 0.40, b = 0.40)
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
    timing <- system.time(status <- system2(binary, c(input, ld_path, output)))
    if (!identical(status, 0L)) stop("JG-0.2 failed for ", name)
    result <- read.delim(output, check.names = FALSE,
                         stringsAsFactors = FALSE)
    if (!identical(result$model_version, "JG-0.2.1")) stop("wrong model version")
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
if (length(failures)) {
    stop("JG-0.2 acceptance failed:\n- ", paste(failures, collapse = "\n- "))
}
cat("JG-0.2.1 adaptive, LD, overlap, scale, and directional tests passed.\n")
