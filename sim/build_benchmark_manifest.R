#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- normalizePath(sub(file_arg, "", script_args[grep(file_arg, script_args)][1]))
source(file.path(dirname(script_path), "benchmark_lib.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) default else args[[idx + 1L]]
}

config_path <- get_arg("--config")
out_path <- get_arg("--out")
benchmark <- get_arg("--benchmark", "all")
if (is.null(config_path) || is.null(out_path)) {
  stop("Usage: Rscript sim/build_benchmark_manifest.R --config <json> --out <tsv> [--benchmark classification|calibration|all]", call. = FALSE)
}

config <- load_config(config_path)
benchmarks <- if (identical(benchmark, "all")) c("classification", "calibration") else benchmark
rows <- list()
ri <- 1L
for (bench in benchmarks) {
  bench_cfg <- config[[bench]]
  for (cell in bench_cfg[["cells"]]) {
    for (rep in seq_len(as.integer(bench_cfg[["replicates"]]))) {
      rows[[ri]] <- data.frame(benchmark = bench, cell = cell[["name"]], replicate = rep, stringsAsFactors = FALSE)
      ri <- ri + 1L
    }
  }
}

ensure_dir(dirname(out_path))
manifest <- do.call(rbind, rows)
write_table(out_path, manifest, c("benchmark", "cell", "replicate"), "\t")
cat(sprintf("Wrote %s\n", out_path))
