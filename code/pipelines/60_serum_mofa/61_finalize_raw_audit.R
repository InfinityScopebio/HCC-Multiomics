#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_dir <- dirname(
  normalizePath(
    sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
      grepl("^--file=", commandArgs(trailingOnly = FALSE))
    ]),
    mustWork = TRUE
  )
)
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "raw_audit_contract.R"))
cfg <- load_pipeline_config()
if (!isTRUE(cfg$is_raw_derived_analysis)) {
  stop("61_finalize_raw_audit.R requires MOFA_ANALYSIS_SOURCE=raw_derived.")
}
require_analysis_mofa_boundary(cfg)
require_historical_mofa_boundary(cfg)
require_local_files(c(cfg$raw_input_path, cfg$analysis_model_path))

input_sha256 <- sha256_file(cfg$raw_input_path)
model_sha256 <- sha256_file(cfg$analysis_model_path)

manifest_path <- file.path(cfg$result_dir, "mofa_training_manifest.tsv")
require_local_files(manifest_path)
manifest <- fread(manifest_path)
if (nrow(manifest) != 117L || any(manifest$status != "complete") ||
    any(manifest$input_sha256 != input_sha256)) {
  stop("Raw training manifest is incomplete or bound to another input.")
}

validate_production_model(
  model_path = cfg$final_model_path,
  metadata_path = cfg$final_model_metadata_path,
  input_sha256 = input_sha256,
  profile = "final_slow2000",
  k = 15L,
  seed = 1L,
  maxiter = 2000L,
  convergence_mode = "slow"
)
for (k in cfg$k_grid) {
  for (seed in cfg$k_scan_seeds) {
    if (!validate_production_model(
      model_path = fast_model_path(cfg, k, seed),
      metadata_path = fast_model_metadata_path(cfg, k, seed),
      input_sha256 = input_sha256,
      profile = "k_scan_fast1000",
      k = k,
      seed = seed,
      maxiter = 1000L,
      convergence_mode = "fast"
    )) {
      stop("Invalid or stale fast model provenance for K=", k, ", seed=", seed)
    }
  }
}
metadata_profiles <- c(
  fast1000_replay = "k_scan_fast1000_replay",
  slow2000_replay = "final_slow2000_replay"
)
for (profile in names(metadata_profiles)) {
  if (!validate_production_model(
    model_path = replay_model_path(cfg, profile, 15L, 1L),
    metadata_path = replay_model_metadata_path(cfg, profile, 15L, 1L),
    input_sha256 = input_sha256,
    profile = unname(metadata_profiles[[profile]]),
    k = 15L,
    seed = 1L,
    maxiter = if (profile == "fast1000_replay") 1000L else 2000L,
    convergence_mode = if (profile == "fast1000_replay") "fast" else "slow"
  )) {
    stop("Invalid or stale replay model provenance: ", profile)
  }
}

audit_paths <- c(
  input = file.path(cfg$audit_dir, "mofa_incdec_input_comparison.tsv"),
  training = file.path(cfg$audit_dir, "mofa_k_scan_variance_audit.tsv"),
  projection = file.path(cfg$audit_dir, "mofa_tcga_projection_audit.tsv"),
  survival = file.path(cfg$audit_dir, "mofa_tcga_survival_audit.tsv"),
  classifier = file.path(cfg$audit_dir, "mofa_factor2_classifier_audit.tsv"),
  classifier_runtime = file.path(cfg$result_dir, "factor2_classifier_runtime.tsv"),
  supplemental = file.path(cfg$audit_dir, "mofa_supplemental_audit.tsv"),
  classifier_seed = file.path(cfg$audit_dir, "mofa_classifier_seed_audit.tsv")
)
require_local_files(audit_paths)
audits <- lapply(audit_paths, fread)
classified <- classify_raw_checks(audits[
  setdiff(names(audits), "classifier_runtime")
])
assert_raw_input_status(audits$input)
assert_raw_repeatability(classified)
for (name in c("projection", "survival", "supplemental")) {
  validate_source_hashes(name, audits, model_sha256, input_sha256)
}
runtime <- audits$classifier_runtime
if (nrow(runtime) != 1L ||
    runtime$final_model_sha256[[1L]] != model_sha256 ||
    runtime$raw_input_sha256[[1L]] != input_sha256) {
  stop("Classifier runtime is bound to another model or input.")
}

production_files <- unique(c(
  list.files(cfg$result_dir, recursive = TRUE, full.names = TRUE),
  list.files(cfg$model_dir, recursive = TRUE, full.names = TRUE)
))
production_files <- production_files[
  file.info(production_files)$isdir %in% FALSE &
    !grepl("/audit/", production_files, fixed = TRUE) &
    !grepl("/logs/", production_files, fixed = TRUE) &
    !grepl("/receipts/", production_files, fixed = TRUE)
]
output_checksums <- data.table(
  path = normalizePath(production_files, mustWork = TRUE),
  bytes = as.numeric(file.info(production_files)$size),
  sha256 = vapply(production_files, sha256_file, character(1))
)
fwrite(output_checksums, file.path(cfg$audit_dir, "mofa_raw_output_checksums.tsv"), sep = "\t")
fwrite(classified, file.path(cfg$audit_dir, "mofa_raw_check_status.tsv"), sep = "\t")

n_differences <- article_difference_count(classified)
article_equivalence <- if (n_differences == 0L) {
  "ARTICLE_EQUIVALENCE_NUMERICALLY_EXACT"
} else {
  "ARTICLE_EQUIVALENCE_NOT_ESTABLISHED"
}
completion <- "complete"
verdict <- data.table(
  completion = completion,
  article_equivalence = article_equivalence,
  execution_mode = ifelse(
    identical(Sys.getenv("MOFA_RUN_RAW_DIAGNOSTIC", "0"), "1"),
    "raw_training_and_downstream",
    "audit_of_existing_raw_outputs"
  ),
  n_checks = nrow(classified),
  n_exact = sum(classified$verdict == "exact"),
  n_diagnostic_difference = n_differences,
  n_fatal = 0L,
  raw_input_sha256 = input_sha256,
  raw_model_sha256 = model_sha256,
  historical_input_sha256 = sha256_file(cfg$article_input_path),
  historical_k_scan_sha256 = sha256_file(cfg$article_k_scan_summary_path),
  historical_model_sha256 = sha256_file(cfg$article_model_path),
  hdf5_byte_equality = audits$training$status[
    audits$training$check_id == "final_k15_file_hash"
  ]
)
fwrite(verdict, file.path(cfg$audit_dir, "mofa_raw_verdict.tsv"), sep = "\t")

report <- c(
  "# Raw-derived MOFA audit",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Completion and equivalence",
  "",
  paste0("- Completion: ", completion, "."),
  paste0("- Article equivalence: ", article_equivalence, "."),
  paste0("- Execution mode: ", verdict$execution_mode, "."),
  paste0("- Checks: ", nrow(classified), "; exact: ", verdict$n_exact,
         "; diagnostic differences: ", n_differences, "; fatal: 0."),
  "- Diagnostic differences are not an overall PASS and are not silently replaced by article values.",
  "",
  "## Provenance",
  "",
  paste0("- Raw input SHA256: `", input_sha256, "`."),
  paste0("- Raw model SHA256: `", model_sha256, "`."),
  paste0("- Historical input SHA256: `", verdict$historical_input_sha256, "`."),
  paste0("- Historical K-scan SHA256: `", verdict$historical_k_scan_sha256, "`."),
  paste0("- Historical model SHA256: `", verdict$historical_model_sha256, "`."),
  "- Complete K=2:20 by six-seed grid, final, and replay metadata were revalidated against the current input.",
  "- Downstream audit and classifier runtime records are bound to the same input and model hashes.",
  "- Historical boundary hashes were revalidated without materializing historical outputs.",
  "",
  "## Interpretation",
  "",
  "This is a raw-derived computation/provenance audit. Article numerical equivalence is established only when every designated numerical check is exact. Matched absolute correlation is diagnostic only; HDF5 byte equality is reported separately."
)
writeLines(report, file.path(cfg$audit_dir, "MOFA_RAW_AUDIT.md"))
writeLines(capture.output(sessionInfo()), file.path(cfg$audit_dir, "mofa_raw_session_info.txt"))

cat(
  "Raw-derived MOFA audit complete; article equivalence=", article_equivalence,
  "; diagnostic differences=", n_differences,
  "; production files=", nrow(output_checksums), ".\n",
  sep = ""
)
