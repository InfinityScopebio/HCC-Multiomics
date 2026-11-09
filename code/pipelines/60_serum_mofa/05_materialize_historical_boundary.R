#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MOFA2)
})

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

expected_columns <- c(
  "role",
  "source_path",
  "materialized_path",
  "bytes",
  "sha256",
  "article_use",
  "contract"
)
manifest <- fread(cfg$historical_boundary_manifest, encoding = "UTF-8")
if (!identical(names(manifest), expected_columns)) {
  stop("Historical MOFA boundary manifest schema has changed.")
}
if (!identical(sort(manifest$role), sort(c("input", "k_scan_summary", "k15_model")))) {
  stop("Historical MOFA boundary manifest must contain exactly three roles.")
}

materialize_record <- function(record) {
  source_path <- file.path(cfg$project_root, record$source_path)
  target_path <- file.path(cfg$project_root, record$materialized_path)
  if (!file.exists(target_path)) {
    require_local_files(source_path)
    if (sha256_file(source_path) != record$sha256) {
      stop("Historical source checksum mismatch: ", source_path)
    }
    dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
    temporary_path <- paste0(target_path, ".tmp")
    on.exit(unlink(temporary_path), add = TRUE)
    if (!file.copy(source_path, temporary_path, overwrite = TRUE)) {
      stop("Failed to materialize historical boundary: ", target_path)
    }
    if (!file.rename(temporary_path, target_path)) {
      stop("Failed to atomically install historical boundary: ", target_path)
    }
  }
  observed_bytes <- as.numeric(file.info(target_path)$size)
  observed_sha256 <- sha256_file(target_path)
  if (observed_bytes != as.numeric(record$bytes) ||
      observed_sha256 != record$sha256) {
    stop(
      "Materialized historical boundary differs from its manifest: ",
      target_path
    )
  }
  list(
    source_path = source_path,
    target_path = target_path,
    observed_bytes = observed_bytes,
    observed_sha256 = observed_sha256
  )
}

records <- lapply(seq_len(nrow(manifest)), function(index) {
  materialize_record(manifest[index])
})
names(records) <- manifest$role

historical_input <- readRDS(records$input$target_path)
if (!is.list(historical_input$views) || length(historical_input$views) != 4L) {
  stop("Historical MOFA input must contain four views.")
}
expected_features <- c(Protein = 49L, Phospho = 17L, LCMS = 766L, GCMS = 85L)
observed_features <- vapply(historical_input$views, nrow, integer(1))
observed_samples <- unique(vapply(historical_input$views, ncol, integer(1)))
if (!identical(observed_features[names(expected_features)], expected_features) ||
    !identical(observed_samples, 54L) ||
    length(historical_input$common) != 54L) {
  stop("Historical MOFA input dimensions have changed.")
}

k_scan <- fread(records$k_scan_summary$target_path)
required_k_scan_columns <- c(
  "K",
  "n_factors",
  "n_active",
  "n_robust",
  "n_strict_robust",
  "total_VE",
  "mean_repro_active"
)
if (!identical(names(k_scan), required_k_scan_columns) ||
    !identical(k_scan$K, 2:20) ||
    nrow(k_scan[K == 15L & abs(total_VE - 220.38) < 1e-12]) != 1L) {
  stop("Historical K-scan semantic contract has changed.")
}

model <- load_model(records$k15_model$target_path, verbose = FALSE)
scores <- get_factors(model)[[1]]
weights <- get_weights(model)
if (!identical(dim(scores), c(54L, 15L)) ||
    length(weights) != 4L ||
    !setequal(rownames(scores), historical_input$common) ||
    !identical(colnames(scores), paste0("Factor", 1:15))) {
  stop("Historical K15 model semantic contract has changed.")
}

receipt <- copy(manifest)
receipt[, resolved_source_path := vapply(records, function(item) {
  normalizePath(item$source_path, mustWork = FALSE)
}, character(1))]
receipt[, resolved_materialized_path := vapply(records, function(item) {
  normalizePath(item$target_path, mustWork = TRUE)
}, character(1))]
receipt[, observed_bytes := vapply(records, `[[`, numeric(1), "observed_bytes")]
receipt[, observed_sha256 := vapply(records, `[[`, character(1), "observed_sha256")]
receipt[, checksum_pass := observed_bytes == bytes & observed_sha256 == sha256]
receipt[, semantic_pass := TRUE]
receipt[, status := ifelse(checksum_pass & semantic_pass, "PASS", "FAIL")]
fwrite(
  receipt,
  file.path(cfg$audit_dir, "mofa_historical_boundary_validation.tsv"),
  sep = "\t"
)

cat(
  "Validated historical MOFA article boundary:",
  "54 samples, four views, K=2:20 summary, and native K15 model.\n"
)
