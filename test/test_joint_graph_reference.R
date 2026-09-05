#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (!length(script_argument)) stop("cannot locate test script")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(root, "research", "joint_graph_reference.R"))

cpp_binary <- file.path(root, "build", "joint_graph_cli")
if (!file.exists(cpp_binary)) stop("build/joint_graph_cli is missing")
output_directory <- file.path(root, "build", "test", "joint_graph")
unlink(output_directory, recursive = TRUE)
dir.create(output_directory, recursive = TRUE)

# The nonzero values are interior seven-node Gauss-Hermite points under the
# frozen priors, preventing numerical grid error from being mistaken for a
# graph-identification failure in these deterministic implementation fixtures.
path_value <- sqrt(2) * 0.70 * 0.8162878828589647
residual_value <- sqrt(2) * 0.175 * 0.8162878828589647
fixtures <- list(
    null = list(seed = 20260910, a = 0, b = 0, c_path = residual_value,
                lambda = 0, q = 0),
    mediation = list(seed = 20260911, a = path_value, b = path_value,
                     c_path = residual_value, lambda = 0, q = 0),
    pleiotropy = list(seed = 20260912, a = path_value, b = 0,
                     c_path = residual_value, lambda = path_value, q = 0.35),
    coexistence = list(seed = 20260913, a = path_value, b = path_value,
                      c_path = residual_value, lambda = path_value, q = 0.35)
)

rows <- vector("list", length(fixtures))
for (index in seq_along(fixtures)) {
    name <- names(fixtures)[index]
    fixture <- fixtures[[index]]
    data <- jg_simulate(
        seed = fixture$seed,
        n_per_role = 60,
        a = fixture$a,
        b = fixture$b,
        c_path = fixture$c_path,
        lambda = fixture$lambda,
        q = fixture$q
    )
    input <- file.path(output_directory, paste0(name, ".tsv"))
    r_output <- file.path(output_directory, paste0(name, ".r.tsv"))
    cpp_output <- file.path(output_directory, paste0(name, ".cpp.tsv"))
    write.table(data, input, sep = "\t", quote = FALSE, row.names = FALSE)

    r_fit <- jg_fit(data)
    jg_write_result(r_fit, r_output)
    elapsed <- system.time({
        status <- system2(cpp_binary, c(input, cpp_output))
    })[["elapsed"]]
    if (!identical(status, 0L)) stop("C++ evaluator failed for ", name)

    r_row <- read.delim(r_output, check.names = FALSE)
    cpp_row <- read.delim(cpp_output, check.names = FALSE)
    if (!identical(
        r_row$identification_scope,
        "CONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY"
    ) || !identical(r_row$identification_scope, cpp_row$identification_scope)) {
        stop(name, " did not preserve the aligned-pleiotropy scope label")
    }
    numeric_columns <- setdiff(
        names(r_row), c("model_version", "identification_scope")
    )
    maximum_difference <- max(abs(
        unlist(r_row[numeric_columns], use.names = FALSE) -
        unlist(cpp_row[numeric_columns], use.names = FALSE)
    ))
    if (maximum_difference > 1e-8) {
        stop(name, " R/C++ maximum absolute difference is ", maximum_difference)
    }

    rows[[index]] <- data.frame(
        fixture = name,
        seed = fixture$seed,
        n_variants = nrow(data),
        true_a = fixture$a,
        true_b = fixture$b,
        true_c = fixture$c_path,
        true_lambda = fixture$lambda,
        true_q = fixture$q,
        PP_XM = r_fit[["PP_XM"]],
        PP_global_MY = r_fit[["PP_global_MY"]],
        PP_nonaligned_P = r_fit[["PP_nonaligned_P"]],
        PP_two_path = r_fit[["PP_two_path"]],
        PP_two_path_plus_P = r_fit[["PP_two_path_plus_P"]],
        r_cpp_max_abs_diff = maximum_difference,
        cpp_elapsed_seconds = elapsed,
        stringsAsFactors = FALSE
    )
}

results <- do.call(rbind, rows)
write.table(results, file.path(output_directory, "development_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

row_for <- function(name) results[results$fixture == name, , drop = FALSE]
null <- row_for("null")
mediation <- row_for("mediation")
pleiotropy <- row_for("pleiotropy")
coexistence <- row_for("coexistence")

# The two parameterizations in the identification proof must have exactly the
# same observation likelihood. This is a required nonidentification result,
# not a failure to optimize the model.
equivalence_data <- jg_simulate(
    seed = 20260914,
    n_per_role = 20,
    a = path_value,
    b = path_value,
    c_path = residual_value
)
mediation_loglik <- jg_loglik(
    equivalence_data, path_value, path_value, residual_value, 0, 0,
    jg_default_config()
)
aligned_loglik <- jg_loglik(
    equivalence_data, path_value, 0,
    residual_value + path_value^2, path_value, 1,
    jg_default_config()
)
if (abs(mediation_loglik - aligned_loglik) > 1e-10) {
    stop("aligned-pleiotropy likelihood equivalence was not preserved")
}

failures <- character()
require_condition <- function(condition, message) {
    if (!isTRUE(condition)) failures <<- c(failures, message)
}
require_condition(null$PP_XM <= 0.20, "null PP_XM exceeded 0.20")
require_condition(null$PP_global_MY <= 0.20,
                  "null PP_global_MY exceeded 0.20")
require_condition(mediation$PP_XM >= 0.80,
                  "mediation PP_XM was below 0.80")
require_condition(mediation$PP_global_MY >= 0.80,
                  "mediation PP_global_MY was below 0.80")
require_condition(mediation$PP_nonaligned_P <= 0.20,
                  "mediation PP_nonaligned_P exceeded 0.20")
require_condition(pleiotropy$PP_XM >= 0.80,
                  "pleiotropy PP_XM was below 0.80")
require_condition(pleiotropy$PP_global_MY <= 0.20,
                  "pleiotropy PP_global_MY exceeded 0.20")
require_condition(pleiotropy$PP_nonaligned_P >= 0.80,
                  "pleiotropy PP_nonaligned_P was below 0.80")
require_condition(coexistence$PP_XM >= 0.70,
                  "coexistence PP_XM was below 0.70")
require_condition(coexistence$PP_global_MY >= 0.70,
                  "coexistence PP_global_MY was below 0.70")
require_condition(coexistence$PP_nonaligned_P >= 0.70,
                  "coexistence PP_nonaligned_P was below 0.70")

print(results[, c("fixture", "PP_XM", "PP_global_MY", "PP_nonaligned_P",
                  "PP_two_path", "r_cpp_max_abs_diff",
                  "cpp_elapsed_seconds")], row.names = FALSE)
if (length(failures)) {
    stop("JG-0.1 development acceptance failed:\n- ",
         paste(failures, collapse = "\n- "))
}
cat("JG-0.1 R/C++ reference and mechanism tests passed.\n")
