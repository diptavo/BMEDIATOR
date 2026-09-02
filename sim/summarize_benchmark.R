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

outdir <- get_arg("--outdir")
rebuild <- has_flag("--rebuild")
if (is.null(outdir)) {
  stop("Usage: Rscript sim/summarize_benchmark.R --outdir <dir> [--rebuild]", call. = FALSE)
}

rebuild_summary_tables <- function(outdir) {
  ensure_dir(file.path(outdir, "summary"))
  config_path <- file.path(outdir, "config.json")
  config <- if (file.exists(config_path)) load_config(config_path) else NULL

  for (benchmark in c("classification", "calibration")) {
    benchmark_root <- file.path(outdir, benchmark)
    if (!dir.exists(benchmark_root)) {
      next
    }
    alias_map <- if (!is.null(config)) scenario_alias_map(config, benchmark) else character()
    summary_scenarios <- NULL
    if (!is.null(config)) {
      scenarios <- benchmark_scenarios(config, benchmark)
      summary_scenarios <- unique(vapply(scenarios, effective_scenario, character(1), alias_map = alias_map))
    }

    rep_dirs <- list.dirs(benchmark_root, recursive = TRUE, full.names = TRUE)
    rep_dirs <- rep_dirs[grepl("/rep_[0-9]+$", rep_dirs)]
    task_frames <- list()
    class_frames <- list()
    conf_frames <- list()
    cal_frames <- list()
    bfdr_frames <- list()
    fdr_power_frames <- list()
    ti <- cli <- cfi <- cai <- bfi <- fpi <- 1L

    for (rep_dir in rep_dirs) {
      truth_path <- file.path(rep_dir, "truth.tsv")
      mediation_path <- file.path(rep_dir, "bmediator.mediation")
      if (file.exists(truth_path) && file.exists(mediation_path)) {
        truth <- read_truth(truth_path)
        results <- read_results(mediation_path)
        merged <- compute_benchmark_metrics(attach_truth(results, truth), alias_map = alias_map)
      } else {
        next
      }
      merged$benchmark <- benchmark
      merged$cell <- basename(dirname(rep_dir))
      merged$replicate <- as.integer(sub("rep_", "", basename(rep_dir)))
      task_frames[[ti]] <- merged
      ti <- ti + 1L
      fdr_power <- summarize_fdr_power(merged)
      fdr_power$benchmark <- benchmark
      fdr_power$cell <- merged$cell[[1]]
      fdr_power$replicate <- merged$replicate[[1]]
      fdr_power_frames[[fpi]] <- fdr_power
      fpi <- fpi + 1L

      if (identical(benchmark, "classification")) {
        class_summary <- summarize_classification(merged, scenarios = summary_scenarios)
        class_summary$by_scenario$cell <- merged$cell[[1]]
        class_summary$by_scenario$replicate <- merged$replicate[[1]]
        class_summary$confusion$cell <- merged$cell[[1]]
        class_summary$confusion$replicate <- merged$replicate[[1]]
        class_frames[[cli]] <- class_summary$by_scenario
        conf_frames[[cfi]] <- class_summary$confusion
        cli <- cli + 1L
        cfi <- cfi + 1L
      } else {
        cal_summary <- summarize_calibration(merged)
        cal_summary$calibration$cell <- merged$cell[[1]]
        cal_summary$calibration$replicate <- merged$replicate[[1]]
        cal_summary$bfdr$cell <- merged$cell[[1]]
        cal_summary$bfdr$replicate <- merged$replicate[[1]]
        cal_frames[[cai]] <- cal_summary$calibration
        bfdr_frames[[bfi]] <- cal_summary$bfdr
        cai <- cai + 1L
        bfi <- bfi + 1L
      }
    }

    summary_dir <- ensure_dir(file.path(outdir, "summary", benchmark))
    if (length(task_frames)) {
      rows <- do.call(rbind, task_frames)
      write_table(file.path(summary_dir, "protein_level_metrics.tsv"), rows, names(rows), "\t")
    }
    if (length(class_frames)) {
      rows <- do.call(rbind, class_frames)
      write_table(file.path(summary_dir, "classification_by_scenario.tsv"), rows, names(rows), "\t")
    }
    if (length(conf_frames)) {
      rows <- do.call(rbind, conf_frames)
      write_table(file.path(summary_dir, "classification_confusion.tsv"), rows, names(rows), "\t")
    }
    if (length(cal_frames)) {
      rows <- do.call(rbind, cal_frames)
      write_table(file.path(summary_dir, "calibration_bins.tsv"), rows, names(rows), "\t")
    }
    if (length(bfdr_frames)) {
      rows <- do.call(rbind, bfdr_frames)
      write_table(file.path(summary_dir, "calibration_bfdr.tsv"), rows, names(rows), "\t")
    }
    if (length(fdr_power_frames)) {
      rows <- do.call(rbind, fdr_power_frames)
      write_table(file.path(summary_dir, "fdr_power_by_threshold.tsv"), rows, names(rows), "\t")
    }
  }
}

plot_classification <- function(summary_dir) {
  by_path <- file.path(summary_dir, "classification", "classification_by_scenario.tsv")
  conf_path <- file.path(summary_dir, "classification", "classification_confusion.tsv")
  if (!file.exists(by_path) || !file.exists(conf_path)) {
    return(invisible(NULL))
  }

  by_rows <- read_table(by_path, "\t")
  conf_rows <- read_table(conf_path, "\t")
  if (nrow(by_rows) == 0 || nrow(conf_rows) == 0) {
    return(invisible(NULL))
  }
  scenarios <- unique(c(conf_rows$true_scenario, conf_rows$pred_scenario, by_rows$true_scenario))
  scenarios <- scenarios[!is.na(scenarios) & nzchar(scenarios)]
  cells <- unique(by_rows$cell)

  agg <- stats::aggregate(accuracy ~ cell + true_scenario, data = by_rows, FUN = mean)
  mat <- matrix(0, nrow = length(scenarios), ncol = length(cells), dimnames = list(scenarios, cells))
  for (i in seq_len(nrow(agg))) {
    mat[agg$true_scenario[[i]], agg$cell[[i]]] <- agg$accuracy[[i]]
  }

  grDevices::png(file.path(summary_dir, "classification", "classification_accuracy.png"), width = 1200, height = 700, res = 160)
  op <- graphics::par(mar = c(5, 5, 3, 2))
  graphics::barplot(
    t(mat),
    beside = TRUE,
    ylim = c(0, 1),
    col = grDevices::hcl.colors(length(cells), "Set 2"),
    ylab = "Classification accuracy",
    xlab = "True scenario",
    names.arg = scenarios,
    main = "BMEDIATOR classification benchmark"
  )
  graphics::legend("topright", legend = cells, fill = grDevices::hcl.colors(length(cells), "Set 2"), bty = "n", cex = 0.8)
  graphics::par(op)
  grDevices::dev.off()

  agg_conf <- stats::aggregate(count ~ true_scenario + pred_scenario, data = conf_rows, FUN = sum)
  pivot <- matrix(0, nrow = length(scenarios), ncol = length(scenarios), dimnames = list(scenarios, scenarios))
  for (i in seq_len(nrow(agg_conf))) {
    pivot[agg_conf$true_scenario[[i]], agg_conf$pred_scenario[[i]]] <- agg_conf$count[[i]]
  }

  grDevices::png(file.path(summary_dir, "classification", "classification_confusion.png"), width = 900, height = 800, res = 160)
  op <- graphics::par(mar = c(5, 5, 3, 4))
  graphics::image(x = seq_along(scenarios), y = seq_along(scenarios), z = t(pivot[nrow(pivot):1, , drop = FALSE]), col = grDevices::hcl.colors(12, "Blues 3"), axes = FALSE, xlab = "Predicted", ylab = "True", main = "Classification confusion")
  graphics::axis(1, at = seq_along(scenarios), labels = scenarios)
  graphics::axis(2, at = seq_along(scenarios), labels = rev(scenarios))
  for (i in seq_along(scenarios)) {
    for (j in seq_along(scenarios)) {
      graphics::text(j, length(scenarios) - i + 1, labels = as.integer(pivot[i, j]), cex = 0.9)
    }
  }
  graphics::par(op)
  grDevices::dev.off()
}

plot_calibration <- function(summary_dir) {
  bins_path <- file.path(summary_dir, "calibration", "calibration_bins.tsv")
  bfdr_path <- file.path(summary_dir, "calibration", "calibration_bfdr.tsv")
  if (!file.exists(bins_path) || !file.exists(bfdr_path)) {
    return(invisible(NULL))
  }

  bin_rows <- read_table(bins_path, "\t")
  bfdr_rows <- read_table(bfdr_path, "\t")
  if (nrow(bin_rows) == 0 || nrow(bfdr_rows) == 0) {
    return(invisible(NULL))
  }
  cells <- unique(bin_rows$cell)

  grDevices::png(file.path(summary_dir, "calibration", "calibration_curve.png"), width = 900, height = 800, res = 160)
  op <- graphics::par(mar = c(5, 5, 3, 2))
  graphics::plot(c(0, 1), c(0, 1), type = "n", xlab = "Mean predicted P(M1)", ylab = "Observed fraction of true M1", main = "Calibration of posterior P(M1)")
  graphics::abline(0, 1, lty = 2, col = "#666666")
  cols <- grDevices::hcl.colors(length(cells), "Dark 3")
  for (i in seq_along(cells)) {
    cell <- cells[[i]]
    cur <- bin_rows[bin_rows$cell == cell, , drop = FALSE]
    cur <- cur[order(cur$bin_mid), , drop = FALSE]
    graphics::lines(cur$mean_pred, cur$observed_m1, type = "b", pch = 16, col = cols[[i]], lwd = 2)
  }
  graphics::legend("topleft", legend = cells, col = cols, lty = 1, pch = 16, bty = "n", cex = 0.8)
  graphics::par(op)
  grDevices::dev.off()

  rank_bins <- unique(bfdr_rows$rank_bin)
  bfdr_plot <- stats::aggregate(
    empirical_fdr ~ cell + rank_bin,
    data = bfdr_rows,
    FUN = mean
  )
  cells <- unique(bfdr_plot$cell)
  grDevices::png(file.path(summary_dir, "calibration", "bfdr_comparison.png"), width = 1200, height = 700, res = 160)
  op <- graphics::par(mar = c(8, 5, 3, 2))
  plot_df <- reshape(
    bfdr_plot[, c("cell", "rank_bin", "empirical_fdr")],
    idvar = "rank_bin",
    timevar = "cell",
    direction = "wide"
  )
  plot_mat <- as.matrix(plot_df[, -1, drop = FALSE])
  rownames(plot_mat) <- plot_df$rank_bin
  bp <- graphics::barplot(t(plot_mat), beside = TRUE, col = grDevices::hcl.colors(ncol(plot_mat), "Set 2"), ylim = c(0, 1), ylab = "FDR", main = "Estimated vs empirical Bayesian FDR by rank bin")
  graphics::axis(1, at = colMeans(bp), labels = rank_bins, las = 2)
  graphics::legend("topright", legend = sub("^empirical_fdr\\.", "", colnames(plot_mat)), fill = grDevices::hcl.colors(ncol(plot_mat), "Set 2"), bty = "n", cex = 0.8)
  graphics::par(op)
  grDevices::dev.off()
}

if (rebuild) {
  rebuild_summary_tables(outdir)
}

summary_dir <- file.path(outdir, "summary")
plot_classification(summary_dir)
plot_calibration(summary_dir)
