#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("usage: prepare_joint_graph_input.R HARMONIZED.tsv SCALES.tsv OUTPUT.tsv")
}
data <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
scales <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)

data_columns <- c(
    "variant", "ld_block", "role", "beta_x", "se_x", "beta_m", "se_m",
    "beta_y", "se_y", "orientation", "orientation_probability",
    "rho_xm", "rho_xy", "rho_my"
)
scale_columns <- c("role", "v_x", "v_m", "v_y")
missing_data <- setdiff(data_columns, names(data))
missing_scales <- setdiff(scale_columns, names(scales))
if (length(missing_data)) {
    stop("harmonized table is missing columns: ", paste(missing_data, collapse = ", "))
}
if (length(missing_scales)) {
    stop("scale table is missing columns: ", paste(missing_scales, collapse = ", "))
}
if (!nrow(data) || any(!nzchar(data$variant)) || anyDuplicated(data$variant)) {
    stop("variant identifiers must be nonempty and unique")
}
if (any(!data$role %in% c("A", "B", "C")) ||
    any(!scales$role %in% c("A", "B", "C")) || anyDuplicated(scales$role)) {
    stop("role values must be unique A/B/C entries in the scale table")
}
if (!all(c("A", "B", "C") %in% scales$role)) {
    stop("scale table must contain A, B, and C")
}
numeric_data <- setdiff(data_columns, c("variant", "ld_block", "role"))
if (any(!vapply(data[numeric_data], function(x) all(is.finite(x)), logical(1)))) {
    stop("harmonized table contains nonfinite numeric values")
}
if (any(!vapply(scales[c("v_x", "v_m", "v_y")],
                function(x) all(is.finite(x)), logical(1))) ||
    any(scales$v_x < 0 | scales$v_m < 0 | scales$v_y <= 0)) {
    stop("scale variances are invalid")
}
if (any(data$se_x <= 0 | data$se_m <= 0 | data$se_y <= 0)) {
    stop("standard errors must be positive")
}
if (any(!data$orientation %in% c(-1, 1)) ||
    any(data$orientation_probability < 0.5 |
        data$orientation_probability > 1)) {
    stop("orientation inputs are invalid")
}
role_blocks <- tapply(data$ld_block, data$role, function(x) length(unique(x)))
a_blocks <- if ("A" %in% names(role_blocks)) role_blocks[["A"]] else 0L
b_blocks <- if ("B" %in% names(role_blocks)) role_blocks[["B"]] else 0L
if (a_blocks < 3L || b_blocks < 3L) {
    stop("at least three independent A-role and B-role blocks are required")
}

scale_index <- match(data$role, scales$role)
data$v_x <- scales$v_x[scale_index]
data$v_m <- scales$v_m[scale_index]
data$v_y <- scales$v_y[scale_index]
output_columns <- c(
    "variant", "ld_block", "role", "beta_x", "se_x", "beta_m", "se_m",
    "beta_y", "se_y", "v_x", "v_m", "v_y", "orientation",
    "orientation_probability", "rho_xm", "rho_xy", "rho_my"
)
dir.create(dirname(args[[3]]), recursive = TRUE, showWarnings = FALSE)
write.table(data[output_columns], args[[3]], sep = "\t", quote = FALSE,
            row.names = FALSE)
cat("Prepared", nrow(data), "joint-model variants\n")
