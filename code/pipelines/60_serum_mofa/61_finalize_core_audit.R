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
cfg <- load_pipeline_config()
require_historical_mofa_boundary(cfg)

audit_paths <- c(
  boundary = file.path(cfg$audit_dir, "mofa_historical_boundary_validation.tsv"),
  article_replay = file.path(cfg$audit_dir, "mofa_article_boundary_replay.tsv"),
  input = file.path(cfg$audit_dir, "mofa_incdec_input_comparison.tsv"),
  projection = file.path(cfg$audit_dir, "mofa_tcga_projection_audit.tsv"),
  survival = file.path(cfg$audit_dir, "mofa_tcga_survival_audit.tsv"),
  classifier = file.path(cfg$audit_dir, "mofa_factor2_classifier_audit.tsv"),
  supplemental = file.path(cfg$audit_dir, "mofa_supplemental_audit.tsv")
)
classifier_seed_path <- file.path(
  cfg$audit_dir,
  "mofa_classifier_seed_audit.tsv"
)
require_local_files(c(audit_paths, classifier_seed_path))
audits <- lapply(audit_paths, fread)

failure_pattern <- "fail|mismatch|error|missing|not_run|incomplete"
failures <- rbindlist(lapply(names(audits), function(audit_name) {
  table <- copy(audits[[audit_name]])
  if (!"status" %in% names(table)) {
    stop("Audit table lacks status column: ", audit_name)
  }
  table[, audit := audit_name]
  table[grepl(failure_pattern, status, ignore.case = TRUE)]
}), fill = TRUE)
if (nrow(failures) > 0L) {
  stop(
    "Core audit contains failures:\n",
    paste(failures$audit, failures$status, sep = "=", collapse = "\n")
  )
}
if (!all(audits$boundary$status == "PASS")) {
  stop("Historical-boundary contract did not pass.")
}

seed_audit <- rbind(
  data.table(
    pipeline = "60_serum_mofa",
    step = c(
      "raw_M0_preprocessing",
      "historical_k_scan_boundary",
      "historical_K15_boundary"
    ),
    script_path = c(
      "code/pipelines/60_serum_mofa/10_rebuild_incdec_input.R",
      "code/pipelines/60_serum_mofa/18_article_k_scan_and_variance.R",
      "code/pipelines/60_serum_mofa/18_article_k_scan_and_variance.R"
    ),
    language = "R",
    algorithm = c(
      "impute::impute.knn and strict-monotonic filter",
      "checksum-locked historical K=2:20 aggregate loading",
      "checksum-locked historical native K15 HDF5 loading"
    ),
    seed = c(cfg$master_seed, NA_integer_, NA_integer_),
    seed_source = c(
      sprintf("LIUF_MASTER_SEED=%d", cfg$master_seed),
      "not applicable; no stochastic fitting is executed",
      "not applicable; no stochastic fitting is executed"
    ),
    rng_kind = c("Mersenne-Twister/Inversion/Rejection", "none", "none"),
    status = c("fixed", "deterministic", "deterministic"),
    notes = c(
      "Raw preprocessing is checked element-wise against the historical input.",
      "The archived stochastic aggregate is an article-equivalence boundary.",
      "The archived stochastic model is an article-equivalence boundary."
    )
  ),
  fread(classifier_seed_path),
  fill = TRUE
)
setorder(seed_audit, step, algorithm)
fwrite(seed_audit, file.path(cfg$audit_dir, "mofa_seed_audit.tsv"), sep = "\t")

lookup_status <- function(audit_name, check_id) {
  table <- audits[[audit_name]]
  requested_id <- check_id
  value <- table[get("check_id") == requested_id, status]
  if (length(value) != 1L) {
    stop("Expected one audit row for ", audit_name, "/", check_id)
  }
  value
}

lineage <- data.table(
  artifact = c(
    "strict_monotonic_input",
    "archived_K_scan",
    "historical_K15_model",
    "TCGA_native_projection",
    "TCGA_survival",
    "Factor2_classifier",
    "Factor2_permutation_B1000",
    "Supplemental_survival_loadings",
    "Supplemental_pathway_tables"
  ),
  article_input = c(
    "raw M0 serum matrices and sample sheet",
    "checksum-locked historical K=2:20 summary",
    "checksum-locked historical native K15 HDF5",
    "historical native K15 HDF5 + local TCGA expression",
    "recomputed native TCGA projection + local clinical data",
    "historical native Factor2 loadings + rebuilt GCMS matrix",
    "recomputed Factor2 MLR with historical serial seed 42",
    "recomputed native survival and loading tables",
    "recomputed directional enrichment tables"
  ),
  validation_target_used_as_computation_input = FALSE,
  article_equivalence_role = c(
    "raw preprocessing verified exactly",
    "irrecoverable historical stochastic aggregate boundary",
    "irrecoverable historical stochastic model boundary",
    rep("downstream recomputation", 6L)
  )
)
fwrite(lineage, file.path(cfg$audit_dir, "mofa_production_lineage.tsv"), sep = "\t")

status <- data.table(
  analysis = c(
    "serum_K_scan",
    "K15_variance",
    "TCGA_native_projection_survival",
    "Factor2_pathways",
    "Factor2_model_comparison",
    "Factor2_top7_MLR",
    "Factor2_permutation",
    "Supplemental_survival_loadings",
    "Supplemental_pathway_tables"
  ),
  rebuild_status = c(
    lookup_status("article_replay", "historical_k_scan_K15"),
    lookup_status("article_replay", "historical_k15_model_checksum"),
    paste(
      lookup_status("projection", "tcga_native_score_magnitudes"),
      lookup_status("survival", "allfactor_survival_statistics"),
      sep = ";"
    ),
    lookup_status("projection", "native_pathway_membership_undirected"),
    lookup_status("classifier", "top_n_model_comparison"),
    lookup_status("classifier", "top7_mlr_loocv"),
    lookup_status("classifier", "permutation_importance"),
    lookup_status("survival", "allfactor_survival_statistics"),
    lookup_status("supplemental", "supplemental_pathway_table")
  ),
  caveat = c(
    "Archived K-scan summary is a checksum-locked stochastic boundary.",
    "Archived K15 HDF5 is a checksum-locked stochastic boundary.",
    "Cox is univariable; no progression sign is applied.",
    "Native 20% loading tails and within-call BH reproduce the historical method.",
    "Feature ranking uses all 54 samples; comparison is exploratory.",
    "Same-cohort non-nested LOOCV is not external validation.",
    "Complete serial B=1000 run preserves historical seed 42.",
    "Best-cutoff log-rank p-value is post-selection.",
    "Directional enrichment is reported as numerical tables only."
  )
)
fwrite(status, file.path(cfg$audit_dir, "mofa_reproduction_status.tsv"), sep = "\t")

production_files <- list.files(
  cfg$result_dir,
  recursive = TRUE,
  full.names = TRUE
)
production_files <- production_files[
  file.info(production_files)$isdir %in% FALSE &
    !grepl("/logs/", production_files, fixed = TRUE)
]
output_checksums <- data.table(
  path = normalizePath(production_files, mustWork = TRUE),
  bytes = as.numeric(file.info(production_files)$size),
  sha256 = vapply(production_files, sha256_file, character(1))
)
fwrite(output_checksums, file.path(cfg$audit_dir, "mofa_output_checksums.tsv"), sep = "\t")

report <- c(
  "# Serum MOFA core computation audit",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Verdict",
  "",
  "PASS. Historical K-scan and native K15 inputs passed checksum validation, followed by local downstream numerical recomputation.",
  "",
  "## Historical boundary",
  "",
  paste0("- Input SHA256: `", sha256_file(cfg$article_input_path), "`."),
  paste0("- K-scan SHA256: `", sha256_file(cfg$article_k_scan_summary_path), "`."),
  paste0("- K15 HDF5 SHA256: `", sha256_file(cfg$article_model_path), "`."),
  "- The raw M0 preprocessing rebuild is checked element-wise against the historical input.",
  "",
  "## Downstream validation",
  "",
  paste0("- TCGA projection: ", lookup_status("projection", "tcga_native_score_magnitudes"), "."),
  paste0("- TCGA survival: ", lookup_status("survival", "allfactor_survival_statistics"), "."),
  paste0("- Factor2 classifier: ", lookup_status("classifier", "top_n_model_comparison"), "; permutation ", lookup_status("classifier", "permutation_importance"), "."),
  "",
  "## Interpretation boundary",
  "",
  "The raw-derived MOFA retraining branch is an optional diagnostic and is not substituted for article values without a separate numerical equivalence assessment. Native factor signs are preserved for downstream calculations."
)
writeLines(report, file.path(cfg$audit_dir, "MOFA_CORE_AUDIT.md"))
writeLines(capture.output(sessionInfo()), file.path(cfg$audit_dir, "mofa_session_info.txt"))

cat(
  "Serum MOFA core computation audit PASS; checksummed",
  nrow(output_checksums),
  "production files.\n"
)
