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
benchmark <- get_arg("--benchmark")
cell_name <- get_arg("--cell")
replicate_id <- as.integer(get_arg("--replicate"))
outdir <- get_arg("--outdir")
binary <- get_arg("--binary", file.path(dirname(dirname(script_path)), "bmediator"))
seed_arg <- get_arg("--seed")

if (any(vapply(list(config_path, benchmark, cell_name, replicate_id, outdir), is.null, logical(1)))) {
  stop("Usage: Rscript sim/run_single_benchmark_task.R --config <json> --benchmark <classification|calibration> --cell <name> --replicate <int> --outdir <dir> [--binary <path>] [--seed <int>]", call. = FALSE)
}

config <- load_config(config_path)
ensure_dir(outdir)
writeLines(readLines(config_path), con = file.path(outdir, "config.json"))

global_cfg <- config[["global"]]
bench_cfg <- config[[benchmark]]
cells <- bench_cfg[["cells"]]
cell_index <- which(vapply(cells, function(x) identical(x[["name"]], cell_name), logical(1)))
if (length(cell_index) != 1L) {
  stop(sprintf("Could not find cell '%s'", cell_name), call. = FALSE)
}
cell_cfg <- cells[[cell_index]]
scenarios <- benchmark_scenarios(config, benchmark)
alias_map <- scenario_alias_map(config, benchmark)
cell_scenario_set <- cell_scenarios(config, benchmark, cell_cfg)

base_seed <- if (is.null(seed_arg)) as.integer(global_cfg[["seed"]]) else as.integer(seed_arg)
set.seed(base_seed + cell_index * 10000L + replicate_id)

if (identical(benchmark, "classification")) {
  scenario_sequence <- classification_sequence(as.integer(bench_cfg[["proteins_per_scenario"]]), cell_scenario_set)
} else {
  scenario_sequence <- calibration_sequence(as.integer(bench_cfg[["proteins_per_replicate"]]), bench_cfg[["scenario_mix"]], cell_scenario_set)
}

proteins <- generate_proteins(benchmark, cell_cfg, scenario_sequence)
files <- write_dataset(outdir, benchmark, cell_name, replicate_id, proteins, global_cfg = global_cfg)
out_prefix <- file.path(dirname(files[["truth"]]), "bmediator")
result_paths <- run_bmediator(binary, files, out_prefix, analysis_config(global_cfg, cell_cfg))

truth <- read_truth(files[["truth"]])
results <- read_results(result_paths[["mediation"]])
merged <- compute_benchmark_metrics(attach_truth(results, truth), alias_map = alias_map)
merged$benchmark <- benchmark
merged$cell <- cell_name
merged$replicate <- replicate_id

task_metrics <- file.path(dirname(files[["truth"]]), "task_metrics.tsv")
write_table(task_metrics, merged, names(merged), "\t")
cat(sprintf("Wrote %s\n", task_metrics))
