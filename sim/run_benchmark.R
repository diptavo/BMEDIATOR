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
has_flag <- function(flag) flag %in% args

config_path <- get_arg("--config")
benchmark <- get_arg("--benchmark", "all")
outdir <- get_arg("--outdir")
binary <- get_arg("--binary", file.path(dirname(dirname(script_path)), "bmediator"))
seed_arg <- get_arg("--seed")
dry_run <- has_flag("--dry-run")

if (is.null(config_path) || is.null(outdir)) {
  stop("Usage: Rscript sim/run_benchmark.R --config <json> --outdir <dir> [--benchmark classification|calibration|all] [--binary <path>] [--seed <int>] [--dry-run]", call. = FALSE)
}

config <- load_config(config_path)
ensure_dir(outdir)
writeLines(readLines(config_path), con = file.path(outdir, "config.json"))
global_cfg <- config[["global"]]
requested <- if (identical(benchmark, "all")) c("classification", "calibration") else benchmark
base_seed <- if (is.null(seed_arg)) as.integer(global_cfg[["seed"]]) else as.integer(seed_arg)

for (bench in requested) {
  cfg <- config[[bench]]
  ensure_dir(file.path(outdir, "summary", bench))
  scenarios <- benchmark_scenarios(config, bench)
  alias_map <- scenario_alias_map(config, bench)
  summary_scenarios <- unique(vapply(scenarios, effective_scenario, character(1), alias_map = alias_map))

  all_metrics <- list()
  all_class <- list()
  all_conf <- list()
  all_cal <- list()
  all_bfdr <- list()
  all_fdr_power <- list()
  mi <- cli <- cfi <- cai <- bfi <- fpi <- 1L

  for (cell_idx in seq_along(cfg[["cells"]])) {
    cell <- cfg[["cells"]][[cell_idx]]
    cell_name <- cell[["name"]]
    cell_scenario_set <- cell_scenarios(config, bench, cell)
    for (replicate_id in seq_len(as.integer(cfg[["replicates"]]))) {
      set.seed(base_seed + cell_idx * 10000L + replicate_id)
      if (identical(bench, "classification")) {
        scenario_sequence <- classification_sequence(as.integer(cfg[["proteins_per_scenario"]]), cell_scenario_set)
      } else {
        scenario_sequence <- calibration_sequence(as.integer(cfg[["proteins_per_replicate"]]), cfg[["scenario_mix"]], cell_scenario_set)
      }

      proteins <- generate_proteins(bench, cell, scenario_sequence)
      files <- write_dataset(outdir, bench, cell_name, replicate_id, proteins, global_cfg = global_cfg)
      if (dry_run) {
        next
      }

      out_prefix <- file.path(dirname(files[["truth"]]), "bmediator")
      result_paths <- run_bmediator(binary, files, out_prefix, analysis_config(global_cfg, cell))
      truth <- read_truth(files[["truth"]])
      results <- read_results(result_paths[["mediation"]])
      merged <- compute_benchmark_metrics(attach_truth(results, truth), alias_map = alias_map)
      merged$benchmark <- bench
      merged$cell <- cell_name
      merged$replicate <- replicate_id
      all_metrics[[mi]] <- merged
      mi <- mi + 1L
      fdr_power <- summarize_fdr_power(merged)
      fdr_power$benchmark <- bench
      fdr_power$cell <- cell_name
      fdr_power$replicate <- replicate_id
      all_fdr_power[[fpi]] <- fdr_power
      fpi <- fpi + 1L

      if (identical(bench, "classification")) {
        class_summary <- summarize_classification(merged, scenarios = summary_scenarios)
        class_summary$by_scenario$cell <- cell_name
        class_summary$by_scenario$replicate <- replicate_id
        class_summary$confusion$cell <- cell_name
        class_summary$confusion$replicate <- replicate_id
        all_class[[cli]] <- class_summary$by_scenario
        all_conf[[cfi]] <- class_summary$confusion
        cli <- cli + 1L
        cfi <- cfi + 1L
      } else {
        cal_summary <- summarize_calibration(merged)
        cal_summary$calibration$cell <- cell_name
        cal_summary$calibration$replicate <- replicate_id
        cal_summary$bfdr$cell <- cell_name
        cal_summary$bfdr$replicate <- replicate_id
        all_cal[[cai]] <- cal_summary$calibration
        all_bfdr[[bfi]] <- cal_summary$bfdr
        cai <- cai + 1L
        bfi <- bfi + 1L
      }
    }
  }

  if (dry_run) {
    next
  }

  summary_dir <- file.path(outdir, "summary", bench)
  if (length(all_metrics)) {
    metrics <- do.call(rbind, all_metrics)
    write_table(file.path(summary_dir, "protein_level_metrics.tsv"), metrics, names(metrics), "\t")
  }
  if (length(all_class)) {
    class_rows <- do.call(rbind, all_class)
    write_table(file.path(summary_dir, "classification_by_scenario.tsv"), class_rows, names(class_rows), "\t")
  }
  if (length(all_conf)) {
    conf_rows <- do.call(rbind, all_conf)
    write_table(file.path(summary_dir, "classification_confusion.tsv"), conf_rows, names(conf_rows), "\t")
  }
  if (length(all_cal)) {
    cal_rows <- do.call(rbind, all_cal)
    write_table(file.path(summary_dir, "calibration_bins.tsv"), cal_rows, names(cal_rows), "\t")
  }
  if (length(all_bfdr)) {
    bfdr_rows <- do.call(rbind, all_bfdr)
    write_table(file.path(summary_dir, "calibration_bfdr.tsv"), bfdr_rows, names(bfdr_rows), "\t")
  }
  if (length(all_fdr_power)) {
    fdr_power_rows <- do.call(rbind, all_fdr_power)
    write_table(file.path(summary_dir, "fdr_power_by_threshold.tsv"), fdr_power_rows, names(fdr_power_rows), "\t")
  }
}
