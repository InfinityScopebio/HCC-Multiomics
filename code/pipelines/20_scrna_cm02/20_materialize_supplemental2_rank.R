#!/usr/bin/env Rscript

# Materialize the accepted Supplemental Figure 2A historical rank survey.

suppressPackageStartupMessages({
  library(data.table)
})

script_path <- normalizePath(
  sub("^--file=", "", grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "config.R"))
config <- ensure_output_dirs(scrna_cm02_config(script_path))

manifest <- fread(config$supplemental2_rank_manifest_path)
required_columns <- c(
  "role", "source_path", "materialized_path", "bytes", "sha256",
  "article_use", "contract"
)
if (!all(required_columns %in% names(manifest))) {
  stop("Supplemental Figure 2A historical manifest schema is invalid.")
}

audit_rows <- vector("list", nrow(manifest))
for (row_index in seq_len(nrow(manifest))) {
  manifest_row <- manifest[row_index]
  source_path <- file.path(config$project_root, manifest_row$source_path)
  materialized_path <- file.path(
    config$project_root,
    manifest_row$materialized_path
  )
  if (!file.exists(source_path)) {
    stop("Missing historical S2A input: ", source_path)
  }
  if (
    sha256_file(source_path) != manifest_row$sha256 ||
      file.info(source_path)$size != manifest_row$bytes
  ) {
    stop("Historical S2A source checksum contract failed: ", source_path)
  }
  dir.create(dirname(materialized_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(materialized_path)) {
    copied <- file.copy(source_path, materialized_path, overwrite = FALSE)
    if (!copied) {
      stop("Failed to materialize S2A input: ", materialized_path)
    }
  }
  if (
    sha256_file(materialized_path) != manifest_row$sha256 ||
      file.info(materialized_path)$size != manifest_row$bytes
  ) {
    stop("Materialized S2A input is not checksum-exact: ", materialized_path)
  }
  audit_rows[[row_index]] <- data.table(
    role = manifest_row$role,
    source_path = manifest_row$source_path,
    materialized_path = manifest_row$materialized_path,
    bytes = manifest_row$bytes,
    sha256 = manifest_row$sha256,
    status = "PASS_CHECKSUM_EXACT"
  )
}

metrics <- fread(config$historical_downstream_paths$rank_metrics)
matrix_snapshot <- fread(
  config$historical_downstream_paths$rank_matrix,
  check.names = FALSE
)
if (
  nrow(metrics) != 27L ||
    !identical(sort(unique(metrics$method)), c("frobenius", "KL", "nsNMF")) ||
    !all(vapply(split(metrics$K, metrics$method), function(values) {
      identical(sort(as.integer(values)), 2:10)
    }, logical(1))) ||
    nrow(matrix_snapshot) != 21L ||
    ncol(matrix_snapshot) != 40L
) {
  stop("Supplemental Figure 2A historical table dimensions changed.")
}

audit <- rbindlist(audit_rows)
fwrite(
  audit,
  file.path(
    config$audit_root,
    "supplemental2_rank_survey_historical_boundary.tsv"
  ),
  sep = "\t"
)
message(
  "Supplemental Figure 2A historical boundary passed: ",
  nrow(audit),
  " checksum-exact inputs; no NMF rank survey was rerun."
)
