#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (!length(script_argument)) stop("cannot locate test script")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."))
binary <- file.path(root, "bmediator-joint")
fixture <- file.path(root, "build", "test", "joint_graph_v02")
output <- file.path(root, "build", "test", "joint_manifest")
null_data <- read.delim(file.path(fixture, "null.tsv"), check.names = FALSE)
harmonized_path <- paste0(output, ".harmonized.tsv")
scale_path <- paste0(output, ".scales.tsv")
prepared_path <- paste0(output, ".prepared.tsv")
write.table(null_data[setdiff(names(null_data), c("v_x", "v_m", "v_y"))],
            harmonized_path, sep = "\t", quote = FALSE, row.names = FALSE)
scale_table <- unique(null_data[c("role", "v_x", "v_m", "v_y")])
write.table(scale_table, scale_path, sep = "\t", quote = FALSE, row.names = FALSE)
prepare_status <- system2("Rscript", c(
    file.path(root, "research", "prepare_joint_graph_input.R"),
    harmonized_path, scale_path, prepared_path
))
if (!identical(prepare_status, 0L)) stop("joint input preparation failed")
manifest <- data.frame(
    protein = c("NULL_PROTEIN", "BAD_PROTEIN", "MEDIATOR_PROTEIN"),
    input = c(prepared_path,
              file.path(fixture, "underidentified.tsv"),
              file.path(fixture, "moderate_mediation.tsv")),
    ld = c(file.path(fixture, "null.ld.tsv"),
           file.path(fixture, "underidentified.ld.tsv"),
           file.path(fixture, "moderate_mediation.ld.tsv")),
    stringsAsFactors = FALSE
)
manifest_path <- paste0(output, ".manifest.tsv")
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
status <- suppressWarnings(system2(
    "Rscript", c(file.path(root, "research", "run_joint_graph_manifest.R"),
                 binary, manifest_path, output, "2")
))
if (!identical(status, 2L)) stop("incomplete manifest did not return status 2")
result <- read.delim(paste0(output, ".joint.tsv"), check.names = FALSE)
failures <- read.delim(paste0(output, ".failures.tsv"), check.names = FALSE)
if (nrow(result) != 2L || nrow(failures) != 1L) stop("unexpected manifest row counts")
if (!identical(result$protein, manifest$protein[c(1, 3)])) {
    stop("successful manifest order was not preserved")
}
if (failures$protein[[1]] != "BAD_PROTEIN") stop("failed protein was not isolated")
if (any(result$model_version != "JG-0.2.7")) stop("wrong manifest model version")
if (!(result$PP_two_path[[1]] < 0.01 && result$PP_two_path[[2]] > 0.80)) {
    stop("manifest results do not separate null and mediation fixtures")
}
if (any(result$family_complete) ||
    any(result$posterior_fdr_status != "UNAVAILABLE_INCOMPLETE_MANIFEST") ||
    any(!is.na(result$selected_bfdr05)) ||
    any(!is.na(result$posterior_rank)) ||
    any(!is.na(result$posterior_cumulative_fdr))) {
    stop("incomplete manifest exposed a family-wide selection")
}

complete_manifest <- manifest[c(1, 3), ]
complete_manifest_path <- paste0(output, ".complete.manifest.tsv")
complete_output <- paste0(output, ".complete")
write.table(complete_manifest, complete_manifest_path, sep = "\t", quote = FALSE,
            row.names = FALSE)
complete_status <- system2(
    "Rscript", c(file.path(root, "research", "run_joint_graph_manifest.R"),
                 binary, complete_manifest_path, complete_output, "2")
)
if (!identical(complete_status, 0L)) stop("complete joint manifest runner failed")
complete_result <- read.delim(paste0(complete_output, ".joint.tsv"),
                              check.names = FALSE)
if (!all(complete_result$family_complete) ||
    any(complete_result$posterior_fdr_status != "AVAILABLE_COMPLETE_MANIFEST") ||
    any(is.na(complete_result$selected_bfdr05)) ||
    any(!is.finite(complete_result$posterior_cumulative_fdr))) {
    stop("complete manifest did not expose an auditable family-wide selection")
}
cat("JG-0.2.7 manifest runner test passed.\n")
