#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(Seurat)
  library(spacexr)
})

source(file.path(dirname(script_path <- normalizePath(
  sub(
    "^--file=",
    "",
    commandArgs(trailingOnly = FALSE)[
      grepl("^--file=", commandArgs(trailingOnly = FALSE))
    ]
  ),
  mustWork = TRUE
)), "common.R"))

seed <- configure_reproducibility()
root <- project_root()
raw_root <- assert_production_path(
  file.path(root, "repo", "data", "public", "raw", "HCC-ST")
)
atlas_path <- assert_production_path(
  file.path(
    root,
    "repo",
    "data",
    "curation",
    "scRNA",
    "combined_annotated_v5.rds"
  )
)
curation_audit_path <- assert_production_path(
  file.path(root, "repo", "audit", "scrna_frozen_curation_objects.tsv")
)
default_output <- file.path(spatial_result_root(), "rctd")
output_root <- normalizePath(
  parse_named_argument("output-root", default_output),
  mustWork = FALSE
)
section_argument <- parse_named_argument("sections", paste(all_sections(), collapse = ","))
sections <- strsplit(section_argument, ",", fixed = TRUE)[[1L]]
invalid_sections <- setdiff(sections, all_sections())
if (length(invalid_sections) > 0L) {
  stop("Unknown sections: ", paste(invalid_sections, collapse = ", "))
}
force <- flag_present("force")
rctd_workers <- suppressWarnings(as.integer(Sys.getenv("RCTD_WORKERS", "2")))
if (is.na(rctd_workers) || rctd_workers < 1L || rctd_workers > 2L) {
  stop("RCTD_WORKERS must be 1 or 2 for the audited 64-GB workstation.")
}

if (!file.exists(atlas_path) || !file.exists(curation_audit_path)) {
  stop(
    "Frozen-curated scRNA reference dependency is missing. Required: ",
    atlas_path,
    " and ",
    curation_audit_path
  )
}
curation_audit <- fread(curation_audit_path)
validate_scrna_curation_audit(curation_audit, atlas_path)

reference_dir <- file.path(output_root, "reference")
weight_dir <- file.path(output_root, "weights")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(weight_dir, recursive = TRUE, showWarnings = FALSE)
reference_path <- file.path(reference_dir, "rctd_reference.rds")
reference_counts_path <- file.path(reference_dir, "reference_cell_counts.tsv")
reference_label_mapping_path <- file.path(
  reference_dir,
  "reference_label_mapping.tsv"
)
reference_source_audit_path <- file.path(
  spatial_audit_root(),
  "spatial_rctd_reference_input_audit.tsv"
)
dir.create(spatial_audit_root(), recursive = TRUE, showWarnings = FALSE)

build_reference <- function() {
  atlas <- readRDS(atlas_path)
  if (!inherits(atlas, "Seurat")) {
    stop("Frozen-curated scRNA reference is not a Seurat object.")
  }
  required_metadata <- c("sub_celltype", "sample", "group", "source")
  missing_metadata <- setdiff(required_metadata, colnames(atlas[[]]))
  if (length(missing_metadata) > 0L) {
    stop(
      "Frozen-curated scRNA reference lacks metadata: ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  if (any(as.character(atlas$sample) %in% c("CC1", "CC2", "CC3"))) {
    stop("Frozen-curated scRNA reference includes prohibited CC1-CC3 samples.")
  }
  label_contract <- validate_rctd_label_contract(atlas$sub_celltype)
  label_contract$atlas_sha256 <- sha256_file(atlas_path)
  fwrite(
    as.data.table(label_contract),
    reference_label_mapping_path,
    sep = "\t"
  )
  DefaultAssay(atlas) <- "RNA"
  atlas <- tryCatch(JoinLayers(atlas), error = function(error) atlas)
  labels <- clean_rctd_label(atlas$sub_celltype)
  counts_all <- GetAssayData(atlas, assay = "RNA", layer = "counts")
  source_check <- audit_rctd_count_sources(
    counts_all,
    atlas$source
  )
  fwrite(
    as.data.table(source_check$audit),
    reference_source_audit_path,
    sep = "\t"
  )
  valid <- !is.na(labels) &
    labels != "" &
    labels != "rm" &
    source_check$keep
  label_counts <- table(labels[valid])
  valid_types <- names(label_counts)[label_counts >= 25L]
  keep <- valid & labels %in% valid_types
  labels <- droplevels(factor(labels[keep]))
  counts <- counts_all[
    ,
    keep,
    drop = FALSE
  ]
  names(labels) <- colnames(counts)
  if (any(abs(counts@x - round(counts@x)) > 1e-8)) {
    stop("RCTD reference still contains noninteger values after source audit.")
  }
  n_umi <- Matrix::colSums(counts)
  reference <- Reference(counts, labels, n_umi, require_int = TRUE)
  saveRDS(reference, reference_path)
  fwrite(
    data.table(
      celltype = levels(labels),
      n_cells = as.integer(table(labels)[levels(labels)])
    ),
    reference_counts_path,
    sep = "\t"
  )
  reference
}

if (file.exists(reference_path) && !force) {
  reference <- readRDS(reference_path)
  if (!file.exists(reference_label_mapping_path)) {
    atlas_for_contract <- readRDS(atlas_path)
    label_contract <- validate_rctd_label_contract(
      atlas_for_contract$sub_celltype
    )
    label_contract$atlas_sha256 <- sha256_file(atlas_path)
    fwrite(
      as.data.table(label_contract),
      reference_label_mapping_path,
      sep = "\t"
    )
  }
} else {
  reference <- build_reference()
}

required_reference_types <- unique(
  c(cm02_rctd_types(), paste0("Hep_C", 1:7))
)
reference_types <- levels(reference@cell_types)
missing_reference_types <- setdiff(required_reference_types, reference_types)
if (length(missing_reference_types) > 0L) {
  stop(
    "RCTD reference lacks required cell states: ",
    paste(missing_reference_types, collapse = ", ")
  )
}

run_section <- function(section) {
  section_index <- match(section, all_sections())
  if (is.na(section_index)) {
    stop("Unable to derive canonical seed index for ", section)
  }
  section_seed <- derive_r_seed(seed, section_index)
  output_rds <- file.path(weight_dir, paste0(section, "_weights.rds"))
  output_csv <- file.path(weight_dir, paste0(section, "_weights.csv.gz"))
  if (file.exists(output_rds) && file.exists(output_csv) && !force) {
    existing <- readRDS(output_rds)
    return(data.table(
      section = section,
      status = "cached_production",
      n_spots = nrow(existing),
      seed = section_seed,
      output_rds = output_rds,
      sha256 = sha256_file(output_rds)
    ))
  }
  set.seed(section_seed)
  section_dir <- assert_production_path(file.path(raw_root, section))
  matrix_path <- file.path(section_dir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(matrix_path)) {
    stop("Missing raw Visium matrix: ", matrix_path)
  }
  spatial <- Load10X_Spatial(
    section_dir,
    filename = "filtered_feature_bc_matrix.h5"
  )
  positions <- read_tissue_positions(section_dir)
  matched <- match(colnames(spatial), positions$barcode)
  if (anyNA(matched)) {
    stop(section, ": raw Visium barcodes are missing from tissue positions.")
  }
  coordinates <- positions[
    matched,
    c("array_row", "array_col"),
    drop = FALSE
  ]
  rownames(coordinates) <- colnames(spatial)
  counts <- GetAssayData(spatial, assay = "Spatial", layer = "counts")
  n_umi <- Matrix::colSums(counts)
  puck <- SpatialRNA(coordinates, counts, n_umi)
  rctd <- create.RCTD(
    puck,
    reference,
    max_cores = 1L,
    UMI_min = 100
  )
  rctd <- run.RCTD(rctd, doublet_mode = "full")
  weights <- as.matrix(rctd@results$weights)
  row_totals <- rowSums(weights)
  if (any(row_totals <= 0)) {
    stop(section, ": RCTD returned a non-positive row sum.")
  }
  weights <- sweep(weights, 1L, row_totals, "/")
  weights <- as.data.frame(weights, check.names = FALSE)
  weights$barcode <- rownames(weights)
  coordinate_table <- data.frame(
    barcode = rownames(coordinates),
    row = coordinates$array_row,
    col = coordinates$array_col,
    stringsAsFactors = FALSE
  )
  pixel_table <- data.frame(
    barcode = positions$barcode,
    imgrow = positions$pxl_row,
    imgcol = positions$pxl_col,
    stringsAsFactors = FALSE
  )
  weights <- merge(
    weights,
    coordinate_table,
    by = "barcode",
    sort = TRUE
  )
  weights <- merge(weights, pixel_table, by = "barcode", sort = TRUE)
  saveRDS(weights, output_rds)
  fwrite(as.data.table(weights), output_csv)
  data.table(
    section = section,
    status = "completed_raw_rctd",
    n_spots = nrow(weights),
    seed = section_seed,
    output_rds = output_rds,
    sha256 = sha256_file(output_rds)
  )
}

section_results <- if (rctd_workers == 1L || length(sections) == 1L) {
  lapply(sections, run_section)
} else {
  parallel::mclapply(
    sections,
    run_section,
    mc.cores = rctd_workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}
failed_sections <- vapply(section_results, inherits, logical(1), "try-error")
if (any(failed_sections)) {
  stop(
    "RCTD failed for: ",
    paste(sections[failed_sections], collapse = ", "),
    "\n",
    paste(as.character(section_results[failed_sections]), collapse = "\n")
  )
}
status <- rbindlist(section_results)
status[, workers := rctd_workers]
status[, run_variant := spatial_run_variant()]
status[, rctd_source := rctd_source()]
status[, article_equivalent := FALSE]
fwrite(status, file.path(output_root, "rctd_section_status.tsv"), sep = "\t")
cat(
  "RCTD completed for ",
  nrow(status),
  " sections using raw Visium and frozen-curated scRNA counts with ",
  rctd_workers,
  " worker(s).\n",
  sep = ""
)
