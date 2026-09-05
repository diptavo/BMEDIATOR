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
status <- system2("Rscript", c(file.path(root, "research", "run_joint_graph_manifest.R"),
                               binary, manifest_path, output, "2"))
if (!identical(status, 0L)) stop("joint manifest runner failed")
result <- read.delim(paste0(output, ".joint.tsv"), check.names = FALSE)
failures <- read.delim(paste0(output, ".failures.tsv"), check.names = FALSE)
if (nrow(result) != 2L || nrow(failures) != 1L) stop("unexpected manifest row counts")
if (!identical(result$protein, manifest$protein[c(1, 3)])) {
    stop("successful manifest order was not preserved")
}
if (failures$protein[[1]] != "BAD_PROTEIN") stop("failed protein was not isolated")
if (any(result$model_version != "JG-0.2.3")) stop("wrong manifest model version")
if (!(result$PP_two_path[[1]] < 0.01 && result$PP_two_path[[2]] > 0.80)) {
    stop("manifest results do not separate null and mediation fixtures")
}
cat("JG-0.2.3 manifest runner test passed.\n")
