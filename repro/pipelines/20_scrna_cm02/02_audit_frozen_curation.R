#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

arguments <- commandArgs(trailingOnly = FALSE)
this_script <- normalizePath(
  sub("^--file=", "", arguments[grepl("^--file=", arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(this_script), "config.R"))
config <- ensure_output_dirs(scrna_cm02_config(this_script))

manifest_path <- file.path(
  config$project_root,
  "repro",
  "manifests",
  "scrna_frozen_curation.tsv"
)
manifest <- fread(manifest_path, colClasses = "character")
required_manifest_columns <- c(
  "object", "filename", "role", "expected_sha256", "expected_features",
  "expected_cells", "expected_samples", "label_column"
)
if (!identical(colnames(manifest), required_manifest_columns)) {
  stop("Frozen-curation manifest schema is invalid.")
}
if (anyDuplicated(manifest$object) || anyDuplicated(manifest$filename)) {
  stop("Frozen-curation manifest contains duplicate objects or filenames.")
}

relative_to_project <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  root_prefix <- paste0(config$project_root, .Platform$file.sep)
  if (startsWith(normalized, root_prefix)) {
    return(substring(normalized, nchar(root_prefix) + 1L))
  }
  normalized
}

audit_rows <- vector("list", nrow(manifest))
atlas_metadata <- NULL
atlas_cells <- NULL
for (index in seq_len(nrow(manifest))) {
  specification <- manifest[index]
  path <- file.path(config$curation_root, specification$filename)
  assert_production_input(config, path)
  observed_sha256 <- sha256_file(path)
  object <- readRDS(path)
  is_seurat <- inherits(object, "Seurat")
  metadata <- if (is_seurat) object[[]] else data.frame()
  required_metadata <- c(
    "sample", "group", "source", "celltype", "sub_celltype"
  )
  metadata_columns_ok <- all(required_metadata %in% colnames(metadata))
  prohibited_hits <- if ("sample" %in% colnames(metadata)) {
    sum(as.character(metadata$sample) %in% config$prohibited_samples)
  } else {
    NA_integer_
  }
  observed_samples <- if ("sample" %in% colnames(metadata)) {
    uniqueN(as.character(metadata$sample))
  } else {
    NA_integer_
  }
  label_column <- specification$label_column
  observed_labels <- if (label_column %in% colnames(metadata)) {
    uniqueN(as.character(metadata[[label_column]]), na.rm = TRUE)
  } else {
    NA_integer_
  }
  expected_features <- as.integer(specification$expected_features)
  expected_cells <- as.integer(specification$expected_cells)
  expected_samples <- as.integer(specification$expected_samples)
  dimensions_ok <- is_seurat &&
    nrow(object) == expected_features &&
    ncol(object) == expected_cells
  current_cells <- if (is_seurat) colnames(object) else character()
  if (identical(specification$object, "atlas")) {
    atlas_cells <- current_cells
  }
  if (is.null(atlas_cells)) {
    stop("The unified atlas must be the first frozen-curation manifest row.")
  }
  cells_in_unified_atlas <- sum(current_cells %in% atlas_cells)
  cells_outside_unified_atlas <- sum(!current_cells %in% atlas_cells)
  samples_ok <- identical(as.integer(observed_samples), expected_samples)
  checksum_ok <- identical(observed_sha256, specification$expected_sha256)
  status <- if (
    is_seurat && metadata_columns_ok && dimensions_ok && samples_ok &&
      checksum_ok && identical(as.integer(prohibited_hits), 0L)
  ) {
    "PASS"
  } else {
    "FAIL"
  }
  audit_rows[[index]] <- data.table(
    object = specification$object,
    role = specification$role,
    path = relative_to_project(path),
    sha256 = observed_sha256,
    expected_sha256 = specification$expected_sha256,
    n_features = if (is_seurat) nrow(object) else NA_integer_,
    n_cells = if (is_seurat) ncol(object) else NA_integer_,
    cells_in_unified_atlas = cells_in_unified_atlas,
    cells_outside_unified_atlas = cells_outside_unified_atlas,
    n_samples = observed_samples,
    label_column = label_column,
    n_labels = observed_labels,
    has_umap = is_seurat && "umap" %in% Reductions(object),
    metadata_columns_ok = metadata_columns_ok,
    prohibited_cc_hits = prohibited_hits,
    curation_boundary = "frozen_expert_annotation",
    status = status
  )
  if (identical(specification$object, "atlas")) {
    atlas_metadata <- as.data.table(metadata, keep.rownames = "cell")
  }
  rm(object, metadata)
  invisible(gc())
}

audit <- rbindlist(audit_rows, fill = TRUE)
fwrite(audit, config$curation_audit_path, sep = "\t")
if (any(audit$status != "PASS")) {
  stop("One or more frozen-curation objects failed integrity validation.")
}

if (is.null(atlas_metadata)) {
  stop("Frozen-curation atlas metadata were not audited.")
}
atlas_total <- nrow(atlas_metadata)
count_rows <- list(
  data.table(
    aggregation = "total",
    sample = "ALL",
    group = "ALL",
    celltype = "ALL",
    sub_celltype = "ALL",
    n_cells = atlas_total
  ),
  atlas_metadata[
    ,
    .(n_cells = .N),
    by = .(sample = as.character(sample), group = as.character(group))
  ][
    ,
    `:=`(
      aggregation = "sample",
      celltype = "ALL",
      sub_celltype = "ALL"
    )
  ],
  atlas_metadata[
    ,
    .(n_cells = .N),
    by = .(group = as.character(group))
  ][
    ,
    `:=`(
      aggregation = "group",
      sample = "ALL",
      celltype = "ALL",
      sub_celltype = "ALL"
    )
  ],
  atlas_metadata[
    ,
    .(n_cells = .N),
    by = .(
      group = as.character(group),
      celltype = as.character(celltype)
    )
  ][
    ,
    `:=`(
      aggregation = "celltype",
      sample = "ALL",
      sub_celltype = "ALL"
    )
  ],
  atlas_metadata[
    ,
    .(n_cells = .N),
    by = .(
      group = as.character(group),
      celltype = as.character(celltype),
      sub_celltype = as.character(sub_celltype)
    )
  ][
    ,
    `:=`(
      aggregation = "sub_celltype",
      sample = "ALL"
    )
  ]
)
count_audit <- rbindlist(count_rows, use.names = TRUE)
count_audit[
  ,
  `:=`(
    denominator_cells = atlas_total,
    fraction_of_atlas = n_cells / atlas_total,
    input_provenance = "frozen_expert_annotation",
    status = "PASS"
  )
]
fwrite(
  count_audit,
  file.path(config$audit_root, "scrna_frozen_curation_counts.tsv"),
  sep = "\t"
)
message(
  "PASS: six frozen expert-curated scRNA objects passed checksum, dimension, ",
  "metadata, and CC1-CC3 exclusion checks."
)
