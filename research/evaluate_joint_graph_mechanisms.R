#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (!length(script_argument)) stop("cannot locate evaluation script")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(root, "research", "joint_graph_reference.R"))

trailing <- commandArgs(trailingOnly = TRUE)
output <- if (length(trailing)) trailing[[1]] else
    file.path(root, "research", "joint_graph_v0_1_development_sweep.tsv")
cpp_binary <- file.path(root, "build", "joint_graph_cli")
if (!file.exists(cpp_binary)) stop("build/joint_graph_cli is missing")

path_value <- sqrt(2) * 0.70 * 0.8162878828589647
residual_value <- sqrt(2) * 0.175 * 0.8162878828589647
scenarios <- list(
    null = list(offset = 0, a = 0, b = 0, c_path = residual_value,
                lambda = 0, q = 0),
    mediation = list(offset = 100, a = path_value, b = path_value,
                     c_path = residual_value, lambda = 0, q = 0),
    pleiotropy = list(offset = 200, a = path_value, b = 0,
                     c_path = residual_value, lambda = path_value, q = 0.35),
    coexistence = list(offset = 300, a = path_value, b = path_value,
                      c_path = residual_value, lambda = path_value, q = 0.35)
)

is_correct <- function(name, result) {
    if (name == "null") {
        return(result$PP_XM <= 0.20 && result$PP_global_MY <= 0.20)
    }
    if (name == "mediation") {
        return(result$PP_XM >= 0.80 && result$PP_global_MY >= 0.80 &&
               result$PP_nonaligned_P <= 0.20)
    }
    if (name == "pleiotropy") {
        return(result$PP_XM >= 0.80 && result$PP_global_MY <= 0.20 &&
               result$PP_nonaligned_P >= 0.80)
    }
    result$PP_XM >= 0.70 && result$PP_global_MY >= 0.70 &&
        result$PP_nonaligned_P >= 0.70
}

temporary <- tempfile("joint_graph_sweep_")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
replicate_rows <- list()
position <- 0L
for (name in names(scenarios)) {
    scenario <- scenarios[[name]]
    for (replicate in seq_len(50)) {
        seed <- 20261000 + replicate + scenario$offset
        data <- jg_simulate(
            seed = seed,
            n_per_role = 60,
            a = scenario$a,
            b = scenario$b,
            c_path = scenario$c_path,
            lambda = scenario$lambda,
            q = scenario$q
        )
        input <- file.path(temporary, "input.tsv")
        result_path <- file.path(temporary, "result.tsv")
        write.table(data, input, sep = "\t", quote = FALSE, row.names = FALSE)
        timing <- system.time(status <- system2(cpp_binary, c(input, result_path)))
        if (!identical(status, 0L)) stop("C++ evaluator failed")
        result <- read.delim(result_path, check.names = FALSE)
        position <- position + 1L
        replicate_rows[[position]] <- data.frame(
            scenario = name,
            replicate = replicate,
            seed = seed,
            PP_XM = result$PP_XM,
            PP_global_MY = result$PP_global_MY,
            PP_nonaligned_P = result$PP_nonaligned_P,
            PP_two_path = result$PP_two_path,
            correct = is_correct(name, result),
            elapsed_seconds = timing[["elapsed"]],
            stringsAsFactors = FALSE
        )
    }
}

replicates <- do.call(rbind, replicate_rows)
summaries <- do.call(rbind, lapply(split(replicates, replicates$scenario), function(x) {
    data.frame(
        scenario = x$scenario[[1]],
        replicates = nrow(x),
        correct = sum(x$correct),
        correct_rate = mean(x$correct),
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
write.table(summaries, output, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(replicates, sub("\\.tsv$", "_replicates.tsv", output),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summaries, row.names = FALSE)
if (any(summaries$correct_rate < 0.80)) {
    failed <- summaries$scenario[summaries$correct_rate < 0.80]
    stop("JG-0.1 development sweep failed for: ", paste(failed, collapse = ", "))
}
cat("JG-0.1 development replication sweep passed.\n")
