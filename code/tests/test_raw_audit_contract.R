#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
test_dir <- dirname(normalizePath(file_arg[[1L]], mustWork = TRUE))
source(file.path(test_dir, "..", "pipelines", "60_serum_mofa", "raw_audit_contract.R"))

expect_error <- function(expression) {
  tryCatch({
    force(expression)
    FALSE
  }, error = function(error) TRUE)
}

classify <- function(audits) {
  classify_raw_checks(audits[setdiff(names(audits), "classifier_runtime")])
}

make_audits <- function() {
  list(
    training = data.table(
      check_id = RAW_TRAINING_CHECK_IDS,
      status = "numeric_and_semantic_exact"
    ),
    input = data.table(
      component = c("Protein", "Phospho", "LCMS", "GCMS", "sample_metadata"),
      status = "exact"
    ),
    projection = data.table(
      check_id = c(
        "tcga_native_score_magnitudes",
        "native_pathway_membership_undirected",
        "native_pathway_padj",
        "local_only_execution"
      ),
      status = "numeric_exact",
      observed = "0",
      source_model_sha256 = "model",
      source_input_sha256 = "input"
    ),
    survival = data.table(
      check_id = c(
        "local_clinical_vs_axis_store",
        "allfactor_survival_statistics",
        "cutpoint_inference",
        "survival_scope"
      ),
      status = "numeric_exact",
      observed = "0",
      source_model_sha256 = "model",
      source_input_sha256 = "input"
    ),
    classifier = data.table(
      check_id = c(
        "top7_feature_selection",
        "top_n_model_comparison",
        "top7_mlr_loocv",
        "permutation_importance",
        "validation_design"
      ),
      status = "exact",
      observed = "0"
    ),
    supplemental = data.table(
      check_id = c(
        "dual_significant_factor_selection",
        "supplemental_raw_top_loadings",
        "supplemental_pathway_table",
        "supplemental_pathway_recomputation"
      ),
      status = "exact",
      observed = "0",
      source_model_sha256 = "model",
      source_input_sha256 = "input"
    ),
    classifier_runtime = data.table(
      final_model_sha256 = "model",
      raw_input_sha256 = "input"
    ),
    classifier_seed = data.table(status = "fixed")
  )
}

valid <- make_audits()
classified <- classify(valid)
stopifnot(all(classified$verdict == "exact"))
assert_raw_input_status(valid$input)
assert_raw_repeatability(classified)
for (name in c("projection", "survival", "supplemental")) {
  validate_source_hashes(name, valid, "model", "input")
}

drift <- make_audits()
drift$training$status[3L] <- "numeric_drift"
classified <- classify(drift)
stopifnot(sum(classified$verdict == "diagnostic_difference") == 1L)

file_hash_difference <- make_audits()
file_hash_difference$training$status[8L] <- "file_hash_differs"
classified <- classify(file_hash_difference)
stopifnot(article_difference_count(classified) == 0L)

semantic_difference <- make_audits()
semantic_difference$training$status[9L] <- "semantic_hash_differs"
classified <- classify(semantic_difference)
stopifnot(article_difference_count(classified) == 1L)

numeric_close <- make_audits()
numeric_close$supplemental$status[4L] <- "numeric_close_rank_exact"
classified <- classify(numeric_close)
stopifnot(article_difference_count(classified) == 1L)

nonfinite <- make_audits()
nonfinite$projection$observed[1L] <- "max_abs_diff=Inf"
stopifnot(expect_error(classify(nonfinite)))

empty <- make_audits()
empty$projection <- empty$projection[0L]
stopifnot(expect_error(classify(empty)))

wrong_downstream_id <- make_audits()
wrong_downstream_id$projection$check_id[1L] <- "other"
stopifnot(expect_error(classify(wrong_downstream_id)))

incomplete <- make_audits()
incomplete$training$status[3L] <- "incomplete"
stopifnot(expect_error(classify(incomplete)))

unknown <- make_audits()
unknown$training$status[3L] <- "unknown_status"
stopifnot(expect_error(classify(unknown)))

repeatability_failure <- make_audits()
repeatability_failure$training$status[10L] <- "numeric_drift"
classified <- classify(repeatability_failure)
stopifnot(expect_error(assert_raw_repeatability(classified)))

bad_input <- make_audits()
bad_input$input$status <- "mismatch"
stopifnot(expect_error(assert_raw_input_status(bad_input$input)))

stale <- make_audits()
stale$projection$source_model_sha256 <- "another_model"
stopifnot(expect_error(
  validate_source_hashes("projection", stale, "model", "input")
))

duplicated_ids <- make_audits()
duplicated_ids$training <- rbind(duplicated_ids$training, duplicated_ids$training[1L])
stopifnot(expect_error(classify(duplicated_ids)))

missing_id <- make_audits()
missing_id$training <- missing_id$training[-1L]
stopifnot(expect_error(classify(missing_id)))

article <- make_audits()[c(
  "input", "projection", "survival", "classifier", "supplemental",
  "classifier_seed"
)]
article$input <- data.table(
  component = c("Protein", "Phospho", "LCMS", "GCMS", "sample_metadata"),
  status = "exact",
  observed = "0"
)
article$boundary <- data.table(
  role = c("input", "k_scan_summary", "k15_model"),
  status = "PASS"
)
article$article_replay <- data.table(
  check_id = c(
    "historical_input_vs_raw_rebuild",
    "historical_k_scan_K15",
    "historical_k15_model_checksum",
    "historical_k15_dimensions"
  ),
  status = c("exact", "exact", "file_exact", "exact"),
  observed = "0"
)
assert_article_audit_contract(article)

empty_article_input <- article
empty_article_input$input <- empty_article_input$input[0L]
stopifnot(expect_error(assert_article_audit_contract(empty_article_input)))

empty_article_boundary <- article
empty_article_boundary$boundary <- empty_article_boundary$boundary[0L]
stopifnot(expect_error(assert_article_audit_contract(empty_article_boundary)))

wrong_article_input <- article
wrong_article_input$input$status[1L] <- "mismatch"
stopifnot(expect_error(assert_article_audit_contract(wrong_article_input)))

unknown_article_status <- article
unknown_article_status$supplemental$status[4L] <- "numeric_close_rank_exact"
stopifnot(expect_error(assert_article_audit_contract(unknown_article_status)))

nonfinite_article <- article
nonfinite_article$survival$observed[2L] <- "max_abs_diff=NA"
stopifnot(expect_error(assert_article_audit_contract(nonfinite_article)))

source(file.path(test_dir, "..", "lib", "R", "runtime_flags.R"))
check_resume <- function(value, expected = NULL, fails = FALSE) {
  name <- "LIUF_TEST_RESUME"
  previous <- Sys.getenv(name, unset = NA_character_)
  on.exit({
    if (is.na(previous)) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, setNames(list(previous), name))
    }
  })
  do.call(Sys.setenv, setNames(list(value), name))
  if (fails) {
    stopifnot(expect_error(resume_enabled(name)))
  } else {
    stopifnot(identical(resume_enabled(name), expected))
  }
}
check_resume("0", FALSE)
check_resume("1", TRUE)
check_resume("true", fails = TRUE)
check_resume("", fails = TRUE)

spatial_script <- file.path(
  test_dir, "..", "pipelines", "50_spatial_cci", "20_recover_cm02.R"
)
spatial_lines <- readLines(spatial_script, warn = FALSE)
retained_assignment <- spatial_lines[
  grepl("^  retained_states <- rownames", spatial_lines)
]
stopifnot(length(retained_assignment) == 1L)
normalized_frequency <- matrix(
  0,
  nrow = 4L,
  ncol = 1L,
  dimnames = list(c("Hep_C4", "Treg", "HSC", "Pericyte"), "sample")
)
state_statistics <- data.table(
  subCluster = sort(rownames(normalized_frequency))
)
state_statistics[, dropped := subCluster == "HSC"]
eval(parse(text = retained_assignment))
stopifnot(identical(retained_states, c("Hep_C4", "Treg", "Pericyte")))
for (spatial_file in c("20_recover_cm02.R", "30_radial_from_raw.R")) {
  lines <- readLines(
    file.path(test_dir, "..", "pipelines", "50_spatial_cci", spatial_file),
    warn = FALSE
  )
  stopifnot(
    sum(grepl(
      "weight_path <- active_rctd_weight_path(section)",
      lines,
      fixed = TRUE
    )) == 1L,
    !any(grepl("^weight_root <-", lines))
  )
}
stopifnot(any(grepl(
  'article_equivalent = identical(rctd_source(), "historical_frozen")',
  spatial_lines,
  fixed = TRUE
)))

cat("Raw and article audit contract fixtures PASS\n")
