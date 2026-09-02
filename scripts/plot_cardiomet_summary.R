args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript plot_cardiomet_summary.R <summary_dir> <plot_dir>")
}

summary_dir <- normalizePath(args[1], mustWork = TRUE)
plot_dir <- args[2]
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(ggplot2))

exposure_order <- c(
  "BMI_GIANT", "BMI_MVP", "WHR_GIANT", "WHRadjBMI_GIANT",
  "SBP_KEATON", "SBP_MVP", "DBP_KEATON", "DBP_MVP", "PP_KEATON",
  "HDL_GLGC", "LDL_GLGC", "TG_GLGC", "TC_GLGC"
)

vascular_outcomes <- c(
  "CHD_FinnGen", "MI_FinnGen", "Stroke_FinnGen",
  "HeartFail_FinnGen", "PADProxy_FinnGen"
)

kidney_outcomes <- c(
  "CKD_CKDGen", "Microalbuminuria_CKDGen", "UACR_CKDGen",
  "RenFail_FinnGen", "KidneyStones_FinnGen", "eGFRCrea_CKDGen", "eGFRCys_CKDGen"
)

run_summary <- read.delim(file.path(summary_dir, "run_summary.tsv"), sep = "\t", check.names = FALSE)
fdr5_hits <- read.delim(file.path(summary_dir, "fdr5_hits.tsv"), sep = "\t", check.names = FALSE)
gene_recurrence <- read.delim(file.path(summary_dir, "gene_recurrence.tsv"), sep = "\t", check.names = FALSE)

heatmap_df <- function(domain_name, outcomes) {
  df <- subset(run_summary, domain == domain_name)
  full <- expand.grid(
    exposure = exposure_order,
    outcome = outcomes,
    stringsAsFactors = FALSE
  )
  merged <- merge(full, df[, c("exposure", "outcome", "fdr5_hits")], by = c("exposure", "outcome"), all.x = TRUE)
  merged$fdr5_hits[is.na(merged$fdr5_hits)] <- 0
  merged$exposure <- factor(merged$exposure, levels = rev(exposure_order))
  merged$outcome <- factor(merged$outcome, levels = outcomes)
  merged
}

plot_heatmap <- function(df, title, out_path) {
  p <- ggplot(df, aes(x = outcome, y = exposure, fill = fdr5_hits)) +
    geom_tile(color = "#f2f2f2", linewidth = 0.4) +
    geom_text(aes(label = fdr5_hits), size = 3.0, colour = "#222222") +
    scale_fill_gradient(low = "#fff7ec", high = "#d7301f") +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = "FDR<5%\nmediators"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  ggsave(out_path, plot = p, width = 9, height = max(4.5, 0.38 * nrow(unique(df["exposure"])) + 1.8), dpi = 220)
}

vascular_df <- heatmap_df("vascular", vascular_outcomes)
kidney_df <- heatmap_df("kidney", kidney_outcomes)

plot_heatmap(
  vascular_df,
  "Vascular outcomes: FDR<5% mediator counts",
  file.path(plot_dir, "vascular_fdr5_heatmap_r.png")
)

plot_heatmap(
  kidney_df,
  "Kidney outcomes: FDR<5% mediator counts",
  file.path(plot_dir, "kidney_fdr5_heatmap_r.png")
)

recurrence_path <- file.path(plot_dir, "top_recurrent_genes_fdr5_r.png")
if (nrow(gene_recurrence) == 0) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No FDR<5% mediators detected", size = 6) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = "Most recurrent FDR<5% mediator genes") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(recurrence_path, plot = p, width = 8, height = 4.5, dpi = 220)
} else {
  top <- head(gene_recurrence, 20)
  top$gene <- factor(top$gene, levels = rev(top$gene))
  p <- ggplot(top, aes(x = fdr5_run_count, y = gene)) +
    geom_col(fill = "#2c7fb8") +
    geom_text(aes(label = fdr5_run_count), hjust = -0.15, size = 3) +
    labs(
      title = "Most recurrent FDR<5% mediator genes",
      x = "Number of exposure-outcome runs",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(recurrence_path, plot = p, width = 8.5, height = 6, dpi = 220)
}
