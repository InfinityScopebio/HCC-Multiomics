#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

script_path <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- sub(
    "^--file=",
    "",
    arguments[grepl("^--file=", arguments)]
  )
  if (length(file_argument) != 1L) {
    stop("Run this script with Rscript.")
  }
  normalizePath(file_argument, mustWork = TRUE)
}

project_root <- function() {
  normalizePath(
    file.path(dirname(script_path()), "..", "..", ".."),
    mustWork = TRUE
  )
}

spatial_run_variant <- function() {
  variant <- Sys.getenv(
    "LIUF_SPATIAL_RUN_VARIANT",
    unset = "article_equivalent"
  )
  allowed <- c("article_equivalent", "raw_rctd_diagnostic")
  if (!variant %in% allowed) {
    stop(
      "LIUF_SPATIAL_RUN_VARIANT must be one of ",
      paste(allowed, collapse = ", "),
      "; got ",
      variant
    )
  }
  variant
}

spatial_result_root <- function() {
  base <- file.path(project_root(), "repo", "results", "50_spatial_cci")
  if (identical(spatial_run_variant(), "raw_rctd_diagnostic")) {
    return(file.path(base, "raw_rctd_diagnostic"))
  }
  base
}

spatial_audit_root <- function() {
  if (identical(spatial_run_variant(), "raw_rctd_diagnostic")) {
    return(file.path(spatial_result_root(), "audit"))
  }
  file.path(project_root(), "repo", "audit")
}

rctd_source <- function() {
  if (identical(spatial_run_variant(), "raw_rctd_diagnostic")) {
    return("raw_recomputed")
  }
  "historical_frozen"
}

source(file.path(
  project_root(),
  "repro",
  "lib",
  "R",
  "deterministic_rng.R"
))

master_seed <- function() {
  value <- Sys.getenv("LIUF_MASTER_SEED", unset = "20260807")
  validate_r_seed(value, "LIUF_MASTER_SEED")
}

configure_reproducibility <- function(
  seed = master_seed(),
  parallel_workers = FALSE
) {
  seed <- validate_r_seed(seed)
  thread_count <- Sys.getenv("OMP_NUM_THREADS", unset = "1")
  parsed_threads <- suppressWarnings(as.integer(thread_count))
  if (is.na(parsed_threads) || parsed_threads <= 0L) {
    stop("OMP_NUM_THREADS must be a positive integer; got: ", thread_count)
  }
  Sys.setenv(
    OMP_NUM_THREADS = thread_count,
    OPENBLAS_NUM_THREADS = thread_count,
    MKL_NUM_THREADS = thread_count,
    BLIS_NUM_THREADS = thread_count,
    VECLIB_MAXIMUM_THREADS = thread_count,
    NUMEXPR_NUM_THREADS = thread_count,
    RCPP_PARALLEL_NUM_THREADS = thread_count
  )
  options(mc.cores = parsed_threads)
  if (isTRUE(parallel_workers)) {
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  } else {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  }
  set.seed(seed)
  invisible(seed)
}

assert_production_path <- function(path) {
  root <- project_root()
  resolved <- normalizePath(path, mustWork = FALSE)
  forbidden <- c(
    file.path(root, "repo", "data", "validation"),
    file.path(root, "repo", "data", "reference", "spatial_cci"),
    file.path(root, "hcc_st"),
    file.path(root, "dat"),
    file.path(root, "combined"),
    file.path(root, "plan"),
    file.path(
      root,
      "repo",
      "results",
      "50_spatial_cci",
      "invalidated_frozen_input_run_20260807"
    )
  )
  hits <- vapply(
    forbidden,
    function(prefix) {
      normalized_prefix <- normalizePath(prefix, mustWork = FALSE)
      resolved == normalized_prefix ||
        startsWith(resolved, paste0(normalized_prefix, .Platform$file.sep))
    },
    logical(1)
  )
  if (any(hits)) {
    stop("Forbidden production input: ", resolved)
  }
  resolved
}

all_sections <- function() {
  c(
    "HCC-1N", "HCC-1L", "HCC-1T",
    "HCC-2N", "HCC-2L", "HCC-2T",
    "HCC-3N", "HCC-3L", "HCC-3T",
    "HCC-4N", "HCC-4L", "HCC-4T",
    "HCC-5A", "HCC-5B", "HCC-5C", "HCC-5D"
  )
}

historical_spatial_manifest_path <- function() {
  file.path(
    project_root(),
    "repro",
    "manifests",
    "spatial_historical_intermediates.tsv"
  )
}

historical_spatial_intermediates <- function() {
  path <- historical_spatial_manifest_path()
  if (!file.exists(path)) {
    stop("Missing historical spatial intermediate manifest: ", path)
  }
  manifest <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expected <- c(
    "role", "section", "source_path", "materialized_path",
    "bytes", "sha256", "article_use"
  )
  if (!identical(colnames(manifest), expected)) {
    stop("Historical spatial intermediate manifest schema mismatch.")
  }
  manifest
}

historical_intermediate_path <- function(role, section = "ALL") {
  manifest <- historical_spatial_intermediates()
  record <- manifest[
    manifest$role == role & manifest$section == section,
    ,
    drop = FALSE
  ]
  if (nrow(record) != 1L) {
    stop(
      "Expected one historical intermediate for ",
      role,
      "/",
      section,
      "; found ",
      nrow(record),
      "."
    )
  }
  path <- assert_production_path(file.path(
    project_root(),
    record$materialized_path[[1L]]
  ))
  if (!file.exists(path)) {
    stop("Missing historical intermediate: ", path)
  }
  path
}

historical_rctd_weight_path <- function(section) {
  if (!section %in% all_sections()) {
    stop("Unknown spatial section: ", section)
  }
  historical_intermediate_path("rctd_weight", section)
}

active_rctd_weight_path <- function(section) {
  if (!section %in% all_sections()) {
    stop("Unknown spatial section: ", section)
  }
  if (identical(rctd_source(), "raw_recomputed")) {
    path <- file.path(
      spatial_result_root(),
      "rctd",
      "weights",
      paste0(section, "_weights.rds")
    )
    if (!file.exists(path)) {
      stop("Missing raw-recomputed RCTD weight: ", path)
    }
    return(assert_production_path(path))
  }
  historical_rctd_weight_path(section)
}

historical_covarnet_model_path <- function() {
  historical_intermediate_path("covarnet_model", "ALL")
}

historical_recovery_barcode_map <- function(recovery) {
  if (!is.data.frame(recovery) || !"section" %in% colnames(recovery)) {
    stop("CM recovery must be a data frame containing section.")
  }
  pieces <- lapply(all_sections(), function(section) {
    section_rows <- which(recovery$section == section)
    weights <- readRDS(historical_rctd_weight_path(section))
    expected_ids <- paste0(section, "_", seq_len(nrow(weights)))
    if (length(section_rows) != nrow(weights) ||
        !identical(rownames(recovery)[section_rows], expected_ids)) {
      stop(section, ": recovery rows do not match historical RCTD order.")
    }
    data.frame(
      spot_id = expected_ids,
      section = section,
      barcode = as.character(weights$barcode),
      stringsAsFactors = FALSE
    )
  })
  mapping <- do.call(rbind, pieces)
  rownames(mapping) <- NULL
  if (!identical(mapping$spot_id, rownames(recovery)) ||
      anyDuplicated(paste(mapping$section, mapping$barcode, sep = "\r"))) {
    stop("Historical recovery-to-barcode mapping is not one-to-one.")
  }
  mapping
}

attach_historical_recovery_barcodes <- function(recovery) {
  if ("barcode" %in% colnames(recovery)) {
    return(recovery)
  }
  mapping <- historical_recovery_barcode_map(recovery)
  output <- recovery
  output$barcode <- mapping$barcode
  leading <- intersect(
    c(
      "CM01", "CM02", "CM03", "CM04", "CM05",
      "section", "patient", "region", "barcode", "row", "col", "maxCM"
    ),
    colnames(output)
  )
  output[, c(leading, setdiff(colnames(output), leading)), drop = FALSE]
}

tumor_sections <- function() {
  c(
    "HCC-1T", "HCC-2T", "HCC-3T", "HCC-4T",
    "HCC-5A", "HCC-5B", "HCC-5C", "HCC-5D"
  )
}

cm02_rctd_types <- function() {
  c(
    "HSC", "Pericyte", "Treg", "IL1B_Mp",
    "SPP1_Mp", "TREM2_Mp", "CXCL10_Mp", "Hep_C4"
  )
}

cm02_covarnet_types <- function() {
  c(
    "HSC", "Pericyte", "Treg", "IL1B+M_",
    "SPP1+M_", "TREM2+M_", "CXCL10+M_", "Hep_C4"
  )
}

clean_rctd_label <- function(labels) {
  labels <- sub("M_$", "Mp", as.character(labels))
  gsub("[^A-Za-z0-9_]", "_", labels)
}

rctd_to_covarnet_label_map <- function() {
  c(
    CD4_Tn = "CD4+Tn",
    CD8_Tem = "CD8+Tem",
    CD8_Temra = "CD8+Temra",
    CD8_Trm = "CD8+Trm",
    CXCL10_Mp = "CXCL10+M_",
    IL1B_Mp = "IL1B+M_",
    SPP1_Mp = "SPP1+M_",
    TREM2_Mp = "TREM2+M_"
  )
}

restore_covarnet_labels <- function(labels) {
  labels <- as.character(labels)
  mapping <- rctd_to_covarnet_label_map()
  matched <- match(labels, names(mapping))
  replace <- !is.na(matched)
  labels[replace] <- unname(mapping[matched[replace]])
  labels
}

validate_rctd_label_contract <- function(labels) {
  original <- sort(unique(as.character(labels)))
  original <- original[
    !is.na(original) & original != "" & original != "rm"
  ]
  cleaned <- clean_rctd_label(original)
  restored <- restore_covarnet_labels(cleaned)
  collision <- duplicated(cleaned) | duplicated(cleaned, fromLast = TRUE)
  status <- ifelse(
    !collision & restored == original,
    "PASS",
    "FAIL"
  )
  audit <- data.frame(
    covarnet_label = original,
    rctd_label = cleaned,
    restored_covarnet_label = restored,
    clean_label_collision = collision,
    status = status,
    stringsAsFactors = FALSE
  )
  if (any(audit$status != "PASS")) {
    failures <- audit[audit$status != "PASS", , drop = FALSE]
    stop(
      "RCTD label cleaning is not injective and reversible: ",
      paste(
        paste0(
          failures$covarnet_label,
          "->",
          failures$rctd_label,
          "->",
          failures$restored_covarnet_label
        ),
        collapse = ", "
      )
    )
  }
  audit
}

validate_scrna_curation_audit <- function(table, atlas_path) {
  required_columns <- c(
    "object", "path", "sha256", "expected_sha256", "n_features",
    "n_cells", "n_samples", "prohibited_cc_hits", "curation_boundary",
    "status"
  )
  if (!all(required_columns %in% colnames(table))) {
    stop("Frozen-curated scRNA audit schema mismatch.")
  }
  expected_objects <- c(
    "atlas", "hepatocyte", "t_cell", "macrophage", "fibroblast", "b_cell"
  )
  if (nrow(table) != length(expected_objects) ||
      !setequal(as.character(table$object), expected_objects) ||
      anyDuplicated(as.character(table$object))) {
    stop("Frozen-curated scRNA audit does not contain the six declared objects.")
  }
  if (!all(as.character(table$status) == "PASS")) {
    stop("Frozen-curated scRNA audit contains non-PASS rows.")
  }
  if (!all(as.character(table$curation_boundary) ==
           "frozen_expert_annotation")) {
    stop("Frozen-curated scRNA audit has an unexpected boundary role.")
  }
  if (any(as.integer(table$prohibited_cc_hits) != 0L)) {
    stop("Frozen-curated scRNA audit contains prohibited CC1-CC3 cells.")
  }
  if (!all(as.character(table$sha256) ==
           as.character(table$expected_sha256))) {
    stop("Frozen-curated scRNA audit contains checksum mismatches.")
  }
  atlas_row <- table[as.character(table$object) == "atlas", , drop = FALSE]
  observed_sha256 <- sha256_file(atlas_path)
  if (!identical(observed_sha256, as.character(atlas_row$sha256[[1L]]))) {
    stop("Frozen-curated scRNA atlas changed after its audit was written.")
  }
  invisible(TRUE)
}

validate_scrna_count_audit <- function(table) {
  expected_columns <- c(
    "stage",
    "aggregation",
    "sample",
    "group",
    "source",
    "label",
    "n_cells",
    "denominator_cells",
    "retention_fraction",
    "input_provenance",
    "status"
  )
  if (!identical(colnames(table), expected_columns)) {
    stop(
      "Raw-derived scRNA count audit schema mismatch: ",
      paste(colnames(table), collapse = ", ")
    )
  }
  if (nrow(table) == 0L) {
    stop("Raw-derived scRNA count audit has no rows.")
  }
  text_columns <- c(
    "stage",
    "aggregation",
    "sample",
    "group",
    "source",
    "label",
    "input_provenance",
    "status"
  )
  audit_frame <- as.data.frame(table, stringsAsFactors = FALSE)
  blank_text <- vapply(
    audit_frame[, text_columns, drop = FALSE],
    function(values) any(is.na(values) | trimws(as.character(values)) == ""),
    logical(1)
  )
  if (any(blank_text)) {
    stop("Raw-derived scRNA count audit has blank required text fields.")
  }
  if (!all(as.character(table$status) == "PASS")) {
    stop("Raw-derived scRNA count audit contains non-PASS rows.")
  }
  key_columns <- c(
    "stage",
    "aggregation",
    "sample",
    "group",
    "source",
    "label"
  )
  if (anyDuplicated(audit_frame[, key_columns, drop = FALSE])) {
    stop("Raw-derived scRNA count audit has duplicate semantic keys.")
  }
  expected_stage_aggregations <- c(
    "raw_input\rtotal",
    "raw_input\rsample",
    "post_qc\rtotal",
    "post_qc\rsample",
    "final_atlas\rtotal",
    "final_atlas\rsample",
    "final_atlas\rgroup",
    "final_atlas\rcelltype",
    "final_atlas\rsub_celltype"
  )
  observed_stage_aggregations <- unique(
    paste(
      as.character(audit_frame$stage),
      as.character(audit_frame$aggregation),
      sep = "\r"
    )
  )
  if (!setequal(
    observed_stage_aggregations,
    expected_stage_aggregations
  )) {
    stop(
      "Raw-derived scRNA count audit stage/aggregation contract mismatch."
    )
  }
  prohibited <- intersect(
    unique(as.character(table$sample)),
    c("CC1", "CC2", "CC3")
  )
  if (length(prohibited) > 0L) {
    stop(
      "Raw-derived scRNA count audit contains prohibited samples: ",
      paste(prohibited, collapse = ", ")
    )
  }
  n_cells <- suppressWarnings(as.numeric(table$n_cells))
  denominators <- suppressWarnings(as.numeric(table$denominator_cells))
  retention <- suppressWarnings(as.numeric(table$retention_fraction))
  if (
    any(!is.finite(n_cells)) ||
      any(!is.finite(denominators)) ||
      any(!is.finite(retention))
  ) {
    stop("Raw-derived scRNA count audit has non-finite numeric values.")
  }
  if (
    any(n_cells < 0) ||
      any(denominators <= 0) ||
      any(n_cells > denominators) ||
      any(retention < 0 | retention > 1)
  ) {
    stop("Raw-derived scRNA count audit has out-of-range values.")
  }
  if (
    any(n_cells != round(n_cells)) ||
      any(denominators != round(denominators))
  ) {
    stop("Raw-derived scRNA count audit cell counts are not integers.")
  }
  if (any(abs(retention - n_cells / denominators) > 1e-12)) {
    stop(
      "Raw-derived scRNA count audit retention_fraction is inconsistent."
    )
  }
  total_rows <- (
    as.character(table$aggregation) == "total" &
      as.character(table$sample) == "ALL" &
      as.character(table$group) == "ALL" &
      as.character(table$source) == "ALL" &
      as.character(table$label) == "ALL"
  )
  total_stages <- as.character(table$stage[total_rows])
  expected_stages <- c("raw_input", "post_qc", "final_atlas")
  if (
    length(total_stages) != length(expected_stages) ||
      anyDuplicated(total_stages) ||
      !setequal(total_stages, expected_stages)
  ) {
    stop(
      "Raw-derived scRNA count audit lacks unique raw/post-QC/final totals."
    )
  }

  valid_sentinels <- list(
    total = (
      as.character(audit_frame$sample) == "ALL" &
        as.character(audit_frame$group) == "ALL" &
        as.character(audit_frame$source) == "ALL" &
        as.character(audit_frame$label) == "ALL"
    ),
    sample = (
      as.character(audit_frame$sample) != "ALL" &
        as.character(audit_frame$group) != "ALL" &
        as.character(audit_frame$source) != "ALL" &
        as.character(audit_frame$label) == "ALL"
    ),
    group = (
      as.character(audit_frame$sample) == "ALL" &
        as.character(audit_frame$group) != "ALL" &
        as.character(audit_frame$source) == "ALL" &
        as.character(audit_frame$label) == "ALL"
    ),
    celltype = (
      as.character(audit_frame$sample) == "ALL" &
        as.character(audit_frame$group) != "ALL" &
        as.character(audit_frame$source) == "ALL" &
        as.character(audit_frame$label) != "ALL"
    ),
    sub_celltype = (
      as.character(audit_frame$sample) == "ALL" &
        as.character(audit_frame$group) != "ALL" &
        as.character(audit_frame$source) == "ALL" &
        as.character(audit_frame$label) != "ALL"
    )
  )
  for (aggregation in names(valid_sentinels)) {
    selected <- as.character(audit_frame$aggregation) == aggregation
    if (
      any(selected) &&
        !all(valid_sentinels[[aggregation]][selected])
    ) {
      stop(
        "Raw-derived scRNA count audit has invalid sentinel fields for ",
        aggregation,
        "."
      )
    }
  }

  audit_frame$._n_cells <- round(n_cells)
  audit_frame$._denominator <- round(denominators)
  sample_rows <- audit_frame[
    as.character(audit_frame$aggregation) == "sample",
    ,
    drop = FALSE
  ]
  sample_rows$._key <- paste(
    as.character(sample_rows$sample),
    as.character(sample_rows$group),
    as.character(sample_rows$source),
    sep = "\r"
  )
  samples_by_stage <- lapply(expected_stages, function(stage) {
    sample_rows[
      as.character(sample_rows$stage) == stage,
      ,
      drop = FALSE
    ]
  })
  names(samples_by_stage) <- expected_stages
  raw_samples <- samples_by_stage[["raw_input"]]
  post_samples <- samples_by_stage[["post_qc"]]
  final_samples <- samples_by_stage[["final_atlas"]]
  if (
    nrow(raw_samples) == 0L ||
      !setequal(raw_samples$._key, post_samples$._key) ||
      !setequal(raw_samples$._key, final_samples$._key)
  ) {
    stop(
      "Raw-derived scRNA sample/group/source identities differ across stages."
    )
  }
  post_samples <- post_samples[
    match(raw_samples$._key, post_samples$._key),
    ,
    drop = FALSE
  ]
  final_samples <- final_samples[
    match(raw_samples$._key, final_samples$._key),
    ,
    drop = FALSE
  ]
  if (
    any(final_samples$._n_cells > post_samples$._n_cells) ||
      any(post_samples$._n_cells > raw_samples$._n_cells)
  ) {
    stop(
      "Raw-derived scRNA sample counts violate final <= post-QC <= raw."
    )
  }
  if (
    any(raw_samples$._denominator != raw_samples$._n_cells) ||
      any(post_samples$._denominator != raw_samples$._n_cells) ||
      any(final_samples$._denominator != post_samples$._n_cells)
  ) {
    stop(
      "Raw-derived scRNA sample denominators do not link adjacent stages."
    )
  }

  total_frame <- audit_frame[total_rows, , drop = FALSE]
  total_n <- setNames(total_frame$._n_cells, total_frame$stage)
  total_denominator <- setNames(
    total_frame$._denominator,
    total_frame$stage
  )
  raw_total <- total_n[["raw_input"]]
  post_total <- total_n[["post_qc"]]
  final_total <- total_n[["final_atlas"]]
  if (!(final_total <= post_total && post_total <= raw_total)) {
    stop(
      "Raw-derived scRNA totals violate final <= post-QC <= raw."
    )
  }
  if (
    total_denominator[["raw_input"]] != raw_total ||
      total_denominator[["post_qc"]] != raw_total ||
      total_denominator[["final_atlas"]] != post_total
  ) {
    stop(
      "Raw-derived scRNA total denominators do not link adjacent stages."
    )
  }
  for (stage in expected_stages) {
    observed_sum <- sum(samples_by_stage[[stage]]$._n_cells)
    if (observed_sum != total_n[[stage]]) {
      stop(
        "Raw-derived scRNA ",
        stage,
        " sample counts do not sum to its total."
      )
    }
  }

  final_rows <- audit_frame[
    as.character(audit_frame$stage) == "final_atlas",
    ,
    drop = FALSE
  ]
  final_aggregations <- c("group", "celltype", "sub_celltype")
  for (aggregation in final_aggregations) {
    rows <- final_rows[
      as.character(final_rows$aggregation) == aggregation,
      ,
      drop = FALSE
    ]
    if (sum(rows$._n_cells) != final_total) {
      stop(
        "Raw-derived scRNA final ",
        aggregation,
        " counts do not sum to the final total."
      )
    }
    if (any(rows$._denominator != final_total)) {
      stop(
        "Raw-derived scRNA final ",
        aggregation,
        " denominators do not equal the final total."
      )
    }
  }

  named_counts_equal <- function(left, right) {
    left <- sort(left)
    right <- sort(right)
    identical(names(left), names(right)) &&
      identical(as.numeric(left), as.numeric(right))
  }
  expected_group_counts <- tapply(
    final_samples$._n_cells,
    as.character(final_samples$group),
    sum
  )
  group_rows <- final_rows[
    as.character(final_rows$aggregation) == "group",
    ,
    drop = FALSE
  ]
  observed_group_counts <- setNames(
    group_rows$._n_cells,
    as.character(group_rows$group)
  )
  if (!named_counts_equal(
    expected_group_counts,
    observed_group_counts
  )) {
    stop(
      "Raw-derived scRNA final group counts do not match final samples."
    )
  }
  for (aggregation in c("celltype", "sub_celltype")) {
    rows <- final_rows[
      as.character(final_rows$aggregation) == aggregation,
      ,
      drop = FALSE
    ]
    observed <- tapply(
      rows$._n_cells,
      as.character(rows$group),
      sum
    )
    if (!named_counts_equal(expected_group_counts, observed)) {
      stop(
        "Raw-derived scRNA final ",
        aggregation,
        " counts do not conserve within group."
      )
    }
  }

  celltype_rows <- final_rows[
    as.character(final_rows$aggregation) == "celltype",
    ,
    drop = FALSE
  ]
  sub_celltype_rows <- final_rows[
    as.character(final_rows$aggregation) == "sub_celltype",
    ,
    drop = FALSE
  ]
  if (any(grepl(
    " / ",
    as.character(celltype_rows$label),
    fixed = TRUE
  ))) {
    stop(
      "Raw-derived scRNA celltype labels contain a sub-cell delimiter."
    )
  }
  sub_label_parts <- strsplit(
    as.character(sub_celltype_rows$label),
    " / ",
    fixed = TRUE
  )
  valid_sub_labels <- vapply(
    sub_label_parts,
    function(parts) {
      length(parts) == 2L &&
        all(nzchar(trimws(parts)))
    },
    logical(1)
  )
  if (!all(valid_sub_labels)) {
    stop(
      "Raw-derived scRNA sub_celltype labels must be ",
      "'celltype / sub_celltype'."
    )
  }
  parent_celltypes <- vapply(
    sub_label_parts,
    function(parts) trimws(parts[[1L]]),
    character(1)
  )
  celltype_key <- paste(
    as.character(celltype_rows$group),
    as.character(celltype_rows$label),
    sep = "\r"
  )
  sub_celltype_key <- paste(
    as.character(sub_celltype_rows$group),
    parent_celltypes,
    sep = "\r"
  )
  expected_parent_counts <- setNames(
    celltype_rows$._n_cells,
    celltype_key
  )
  observed_parent_counts <- tapply(
    sub_celltype_rows$._n_cells,
    sub_celltype_key,
    sum
  )
  if (!named_counts_equal(
    expected_parent_counts,
    observed_parent_counts
  )) {
    stop(
      "Raw-derived scRNA sub_celltype counts do not conserve within ",
      "group and parent celltype."
    )
  }
  invisible(TRUE)
}

audit_rctd_count_sources <- function(
  counts,
  sources,
  allowed_noninteger_sources = "GSE115469",
  tolerance = 1e-8
) {
  sources <- as.character(sources)
  if (length(sources) != ncol(counts)) {
    stop("Count matrix columns and source labels are not aligned.")
  }
  observed_sources <- sort(unique(sources))
  rows <- lapply(observed_sources, function(source) {
    indices <- which(sources == source)
    values <- counts[, indices, drop = FALSE]@x
    noninteger <- abs(values - round(values)) > tolerance
    data.frame(
      source = source,
      n_cells = length(indices),
      n_nonzero = length(values),
      n_noninteger = sum(noninteger),
      noninteger_fraction = if (length(values) == 0L) {
        0
      } else {
        mean(noninteger)
      },
      integer_count_source = !any(noninteger),
      used_in_rctd_reference = !any(noninteger),
      exclusion_reason = if (any(noninteger)) {
        "source matrix is noninteger normalized expression, not raw UMI counts"
      } else {
        ""
      },
      stringsAsFactors = FALSE
    )
  })
  audit <- do.call(rbind, rows)
  noninteger_sources <- audit$source[!audit$integer_count_source]
  unexpected <- setdiff(noninteger_sources, allowed_noninteger_sources)
  if (length(unexpected) > 0L) {
    stop(
      "Unexpected noninteger scRNA count sources: ",
      paste(unexpected, collapse = ", ")
    )
  }
  missing_expected <- setdiff(
    allowed_noninteger_sources,
    noninteger_sources
  )
  if (length(missing_expected) > 0L) {
    stop(
      "Expected normalized-only source was not detected as noninteger: ",
      paste(missing_expected, collapse = ", ")
    )
  }
  list(
    audit = audit,
    keep = !sources %in% noninteger_sources
  )
}

read_tissue_positions <- function(section_dir) {
  path <- file.path(section_dir, "spatial", "tissue_positions_list.csv")
  if (!file.exists(path)) {
    stop("Missing tissue positions: ", path)
  }
  positions <- read.csv(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(positions) < 6L) {
    stop("Tissue position table has fewer than six columns: ", path)
  }
  colnames(positions)[1:6] <- c(
    "barcode", "in_tissue", "array_row", "array_col",
    "pxl_row", "pxl_col"
  )
  positions
}

parse_named_argument <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  match <- arguments[startsWith(arguments, prefix)]
  if (length(match) == 0L) {
    return(default)
  }
  if (length(match) != 1L) {
    stop("Argument appears more than once: --", name)
  }
  sub(prefix, "", match, fixed = TRUE)
}

flag_present <- function(name) {
  paste0("--", name) %in% commandArgs(trailingOnly = TRUE)
}

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    c("-a", "256", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unable to hash ", path, ": ", paste(output, collapse = "\n"))
  }
  sub("\\s+.*$", "", output[[1L]])
}
