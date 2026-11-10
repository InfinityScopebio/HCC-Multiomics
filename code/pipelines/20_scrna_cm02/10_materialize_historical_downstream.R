#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
})

all_arguments <- commandArgs(trailingOnly = FALSE)
this_script <- normalizePath(
  sub("^--file=", "", all_arguments[grepl("^--file=", all_arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(this_script), "config.R"))
config <- ensure_output_dirs(scrna_cm02_config(this_script))

manifest <- fread(config$historical_downstream_manifest_path)
required_columns <- c(
  "role", "source_path", "materialized_path", "bytes", "sha256",
  "article_use", "contract"
)
if (!all(required_columns %in% names(manifest))) {
  stop("Historical downstream manifest schema is invalid.")
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
    stop("Missing historical article input: ", source_path)
  }
  source_hash <- sha256_file(source_path)
  source_bytes <- file.info(source_path)$size
  if (
    source_hash != manifest_row$sha256 ||
      source_bytes != manifest_row$bytes
  ) {
    stop("Historical source checksum contract failed: ", source_path)
  }
  if (normalizePath(source_path) != normalizePath(
    materialized_path,
    mustWork = FALSE
  )) {
    dir.create(
      dirname(materialized_path),
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!file.exists(materialized_path)) {
      copied <- file.copy(source_path, materialized_path, overwrite = FALSE)
      if (!copied) {
        stop("Failed to materialize historical input: ", materialized_path)
      }
    }
  }
  if (
    !file.exists(materialized_path) ||
      sha256_file(materialized_path) != manifest_row$sha256 ||
      file.info(materialized_path)$size != manifest_row$bytes
  ) {
    stop("Materialized historical input is not checksum-exact: ", materialized_path)
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

cm02_members <- config$cm02_anchor_states
nmf_model <- readRDS(config$historical_downstream_paths$nmf)
nmf_basis <- NMF::basis(nmf_model)
if (
  !identical(dim(nmf_basis), c(21L, 5L)) ||
    !"CM02" %in% colnames(nmf_basis) ||
    !all(cm02_members %in% rownames(nmf_basis))
) {
  stop("Historical CoVarNet K=5 model contract failed.")
}

tcga_proportions <- fread(
  config$historical_downstream_paths$tcga_proportions
)
bulk_frequency <- fread(config$historical_downstream_paths$bulk_frequency)
bulk_pseudotime <- fread(
  config$historical_downstream_paths$bulk_pseudotime
)
if (
  nrow(tcga_proportions) != 423L ||
    nrow(bulk_frequency) != 1104L ||
    nrow(bulk_pseudotime) != 1087L
) {
  stop("Historical downstream table dimensions changed.")
}
for (table_value in list(tcga_proportions, bulk_frequency)) {
  if (!all(c("sampleID", cm02_members) %in% names(table_value))) {
    stop("Historical frequency table is missing CM02 members.")
  }
  if (anyDuplicated(table_value$sampleID)) {
    stop("Historical frequency table contains duplicate sampleID values.")
  }
}
if (
  !all(c("sampleID", "Group", "Cohort", "Pseudotime") %in%
    names(bulk_pseudotime)) ||
    anyDuplicated(bulk_pseudotime$sampleID) ||
    length(intersect(
      bulk_pseudotime$sampleID,
      bulk_frequency$sampleID
    )) != 1087L
) {
  stop("Historical bulk trajectory join contract failed.")
}

audit <- rbindlist(audit_rows)
fwrite(
  audit,
  file.path(config$audit_root, "cm02_historical_downstream_boundary.tsv"),
  sep = "\t"
)
message(
  "Historical CM02 downstream boundary passed: ",
  nrow(audit),
  " checksum-exact inputs; no RCTD, UMAP, NMF, or deconvolution refit."
)
