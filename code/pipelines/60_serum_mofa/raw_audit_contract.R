RAW_TRAINING_CHECK_IDS <- c(
  "raw_derived_training_grid",
  "fast1000_k_scan_counts",
  "fast1000_k_scan_numeric",
  "final_k15_native_scores",
  "final_k15_native_weights",
  "final_k15_native_variance",
  "final_k15_factor_recovery",
  "final_k15_file_hash",
  "final_k15_semantic_hash",
  "fast1000_repeatability",
  "slow2000_repeatability"
)

RAW_DOWNSTREAM_CHECK_IDS <- list(
  projection = c(
    "tcga_native_score_magnitudes",
    "native_pathway_membership_undirected",
    "native_pathway_padj",
    "local_only_execution"
  ),
  survival = c(
    "local_clinical_vs_axis_store",
    "allfactor_survival_statistics",
    "cutpoint_inference",
    "survival_scope"
  ),
  classifier = c(
    "top7_feature_selection",
    "top_n_model_comparison",
    "top7_mlr_loocv",
    "permutation_importance",
    "validation_design"
  ),
  supplemental = c(
    "dual_significant_factor_selection",
    "supplemental_raw_top_loadings",
    "supplemental_pathway_table",
    "supplemental_pathway_recomputation"
  )
)

RAW_EXACT_STATUSES <- c(
  "exact", "exact_subset", "numeric_exact", "file_exact", "semantic_exact",
  "numeric_and_semantic_exact", "fixed", "deterministic", "PASS",
  "post_selection_with_adjusted_audit", "informational"
)

RAW_DIAGNOSTIC_STATUSES <- c(
  "numeric_drift", "file_hash_differs", "semantic_hash_differs",
  "numeric_mismatch", "mismatch", "caveat", "numeric_close_rank_exact"
)

validate_audit_tables <- function(audits, id_columns = list()) {
  if (length(audits) == 0L || is.null(names(audits)) || any(!nzchar(names(audits)))) {
    stop("Audits must be a named, nonempty collection of tables.")
  }
  for (name in names(audits)) {
    table <- audits[[name]]
    if (!data.table::is.data.table(table) || nrow(table) == 0L) {
      stop("Audit is empty or is not a data.table: ", name)
    }
    if (!"status" %in% names(table)) {
      stop("Audit lacks a status column: ", name)
    }
    if (anyNA(table$status) || any(!nzchar(table$status))) {
      stop("Audit contains missing statuses: ", name)
    }
    id_column <- id_columns[[name]]
    if (is.null(id_column) && "check_id" %in% names(table)) {
      id_column <- "check_id"
    }
    if (!is.null(id_column)) {
      if (!id_column %in% names(table)) {
        stop("Audit lacks its identifier column: ", name, "/", id_column)
      }
      identifiers <- table[[id_column]]
      if (anyNA(identifiers) || any(!nzchar(identifiers)) ||
          anyDuplicated(identifiers)) {
        stop("Audit contains missing or duplicate identifiers: ", name)
      }
    }
    if ("observed" %in% names(table)) {
      observed_text <- as.character(table$observed)
      if (anyNA(observed_text)) {
        stop("Audit contains missing comparison metrics: ", name)
      }
      numeric_values <- suppressWarnings(as.numeric(observed_text))
      embedded_nonfinite <- grepl(
        "(^|[^A-Za-z])(Inf|NaN|NA)([^A-Za-z]|$)",
        observed_text,
        perl = TRUE
      )
      if (any(is.infinite(numeric_values) | is.nan(numeric_values), na.rm = TRUE) ||
          any(embedded_nonfinite)) {
        stop("Audit contains non-finite comparison metrics: ", name)
      }
    }
  }
  invisible(TRUE)
}

assert_expected_ids <- function(table, name, column, expected) {
  if (!column %in% names(table) || !setequal(table[[column]], expected)) {
    stop("Audit identifiers are incomplete or unexpected: ", name)
  }
  invisible(TRUE)
}

classify_raw_checks <- function(audits) {
  validate_audit_tables(audits)
  if (!setequal(audits$training$check_id, RAW_TRAINING_CHECK_IDS)) {
    stop("Raw training audit check IDs are incomplete or unexpected.")
  }
  for (name in names(RAW_DOWNSTREAM_CHECK_IDS)) {
    if (!setequal(
      audits[[name]]$check_id,
      RAW_DOWNSTREAM_CHECK_IDS[[name]]
    )) {
      stop("Downstream audit check IDs are incomplete or unexpected: ", name)
    }
  }
  classified <- data.table::rbindlist(lapply(names(audits), function(name) {
    table <- data.table::copy(audits[[name]])
    table[, audit := name]
    table[, check := if ("check_id" %in% names(table)) check_id else "runtime"]
    table[, verdict := ifelse(
      status %in% RAW_EXACT_STATUSES,
      "exact",
      ifelse(status %in% RAW_DIAGNOSTIC_STATUSES, "diagnostic_difference", "fatal")
    )]
    table
  }), fill = TRUE)

  if (any(classified$verdict == "fatal")) {
    fatal <- classified[verdict == "fatal"]
    stop(
      "Raw-derived audit contains fatal or unknown statuses: ",
      paste(fatal$audit, fatal$check, fatal$status, sep = "/", collapse = "; ")
    )
  }
  classified
}

assert_raw_input_status <- function(input_audit) {
  validate_audit_tables(
    list(input = input_audit),
    id_columns = list(input = "component")
  )
  assert_expected_ids(
    input_audit,
    "input",
    "component",
    c("Protein", "Phospho", "LCMS", "GCMS", "sample_metadata")
  )
  if (any(input_audit$status != "exact")) {
    stop("Raw input does not match the declared historical input contract.")
  }
  invisible(TRUE)
}

assert_raw_repeatability <- function(classified) {
  repeatability_rows <- classified[
    audit == "training" &
      check %in% c("fast1000_repeatability", "slow2000_repeatability")
  ]
  if (any(repeatability_rows$status != "numeric_and_semantic_exact")) {
    stop("Same-configuration repeatability failed; this is not expected article drift.")
  }
  invisible(TRUE)
}

article_difference_count <- function(classified) {
  sum(
    classified$verdict == "diagnostic_difference" &
      classified$check != "final_k15_file_hash"
  )
}

validate_source_hashes <- function(name, audits, model_sha256, input_sha256) {
  table <- audits[[name]]
  if (!all(c("source_model_sha256", "source_input_sha256") %in% names(table)) ||
      any(table$source_model_sha256 != model_sha256) ||
      any(table$source_input_sha256 != input_sha256)) {
    stop("Downstream audit is stale or mixed with another model/input: ", name)
  }
  invisible(TRUE)
}

assert_article_audit_contract <- function(audits) {
  validate_audit_tables(
    audits,
    id_columns = list(input = "component", boundary = "role")
  )
  assert_expected_ids(
    audits$input,
    "input",
    "component",
    c("Protein", "Phospho", "LCMS", "GCMS", "sample_metadata")
  )
  assert_expected_ids(
    audits$boundary,
    "boundary",
    "role",
    c("input", "k_scan_summary", "k15_model")
  )
  assert_expected_ids(
    audits$article_replay,
    "article_replay",
    "check_id",
    c(
      "historical_input_vs_raw_rebuild",
      "historical_k_scan_K15",
      "historical_k15_model_checksum",
      "historical_k15_dimensions"
    )
  )
  for (name in names(RAW_DOWNSTREAM_CHECK_IDS)) {
    assert_expected_ids(
      audits[[name]],
      name,
      "check_id",
      RAW_DOWNSTREAM_CHECK_IDS[[name]]
    )
  }
  if (any(audits$input$status != "exact")) {
    stop("Article input preprocessing is not numerically exact.")
  }
  if (any(audits$boundary$status != "PASS")) {
    stop("Historical-boundary contract did not pass.")
  }
  invalid <- data.table::rbindlist(lapply(names(audits), function(name) {
    table <- data.table::copy(audits[[name]])
    table[, audit := name]
    table[!status %in% RAW_EXACT_STATUSES]
  }), fill = TRUE)
  if (nrow(invalid) > 0L) {
    stop(
      "Article audit contains unknown or non-equivalent statuses: ",
      paste(invalid$audit, invalid$status, sep = "=", collapse = "; ")
    )
  }
  invisible(TRUE)
}
