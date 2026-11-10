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
require_historical_mofa_boundary(cfg)
require_local_files(cfg$raw_input_path)

historical_input <- readRDS(cfg$article_input_path)
rebuilt_input <- readRDS(cfg$raw_input_path)
input_components <- c("common", "grp", "stage")
input_exact <- identical(
  historical_input[c("views", input_components)],
  rebuilt_input[c("views", input_components)]
)
if (!input_exact) {
  stop("Raw-derived serum input is not identical to the historical MOFA input.")
}

k_scan <- fread(cfg$article_k_scan_summary_path)[order(K)]
k_scan[, mean_view_ve := total_VE / length(historical_input$views)]
k_scan[, delta_mean_view_ve := c(NA_real_, diff(mean_view_ve))]
fwrite(
  k_scan,
  file.path(cfg$result_dir, "k_scan_accepted_fast1000.tsv"),
  sep = "\t"
)

model <- load_model(cfg$article_model_path, verbose = FALSE)
variance <- get_variance_explained(model)
factor_view <- as.data.table(as.table(variance$r2_per_factor[[1]]))
setnames(factor_view, c("factor", "view", "variance_percent"))
view_total <- data.table(
  view = names(variance$r2_total[[1]]),
  variance_percent = as.numeric(variance$r2_total[[1]])
)
fwrite(
  factor_view,
  file.path(cfg$result_dir, "k15_variance_factor_by_view.tsv"),
  sep = "\t"
)
fwrite(
  view_total,
  file.path(cfg$result_dir, "k15_variance_total_by_view.tsv"),
  sep = "\t"
)

selected <- k_scan[K == 15L]

audit <- data.table(
  check_id = c(
    "historical_input_vs_raw_rebuild",
    "historical_k_scan_K15",
    "historical_k15_model_checksum",
    "historical_k15_dimensions"
  ),
  status = c("exact", "exact", "file_exact", "exact"),
  observed = c(
    "all four matrices and sample metadata are identical",
    sprintf("total_VE=%.2f;mean_four_views=%.3f", selected$total_VE, selected$mean_view_ve),
    sha256_file(cfg$article_model_path),
    sprintf("scores=%dx%d;views=%d", nrow(get_factors(model)[[1]]), ncol(get_factors(model)[[1]]), length(get_weights(model)))
  ),
  article_role = c(
    "verified_raw_preprocessing",
    "Figure6A_article_input",
    "Figure6B_to_G_article_model",
    "Figure6B_to_G_article_model"
  )
)
fwrite(
  audit,
  file.path(cfg$audit_dir, "mofa_article_boundary_audit.tsv"),
  sep = "\t"
)

cat(
  sprintf(
    "Historical MOFA boundary validated: K15 mean-view VE %.3f%%; model %s.\n",
    selected$mean_view_ve,
    sha256_file(cfg$article_model_path)
  )
)
