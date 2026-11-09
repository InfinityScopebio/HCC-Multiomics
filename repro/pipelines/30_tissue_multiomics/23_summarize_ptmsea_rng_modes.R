#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
project_root <- normalizePath(
  file.path(dirname(normalizePath(script_path)), "..", "..", ".."),
  mustWork = TRUE
)
audit_root <- file.path(project_root, "repo", "audit", "ptmsea_rng_modes")
historical_path <- file.path(
  project_root,
  "f4_phospho_unbiased",
  "PTMSEA_official_v2.0.0_all.csv"
)
production_path <- file.path(
  project_root,
  "repo",
  "results",
  "30_tissue_multiomics",
  "ptmsea",
  "PTMSEA_official_v2.0.0_all.csv"
)
dependency_paths <- c(
  gmt_sha256 = file.path(
    project_root,
    "repo",
    "data",
    "reference",
    "PTMsigDB",
    "ptm.sig.db.all.uniprot.human.v2.0.0.gmt"
  ),
  protein_sha256 = file.path(
    project_root,
    "repo",
    "results",
    "00_preprocess",
    "multiomics",
    "dat",
    "tissue_protein.rds"
  ),
  phospho_sha256 = file.path(
    project_root,
    "repo",
    "results",
    "00_preprocess",
    "multiomics",
    "dat",
    "tissue_phospho.rds"
  ),
  differential_sha256 = file.path(
    project_root,
    "repo",
    "results",
    "30_tissue_multiomics",
    "ptmsea",
    "phosphosite_paired_limma.csv"
  ),
  worker_script_sha256 = file.path(
    project_root,
    "repro",
    "pipelines",
    "30_tissue_multiomics",
    "22_probe_ptmsea_rng_mode.R"
  ),
  production_script_sha256 = file.path(
    project_root,
    "repro",
    "pipelines",
    "30_tissue_multiomics",
    "20_ptmsea.R"
  ),
  orchestrator_script_sha256 = file.path(
    project_root,
    "repro",
    "pipelines",
    "30_tissue_multiomics",
    "22_run_ptmsea_rng_probe.R"
  ),
  summarizer_script_sha256 = normalizePath(script_path, mustWork = TRUE),
  run_all_script_sha256 = file.path(
    project_root,
    "repro",
    "pipelines",
    "30_tissue_multiomics",
    "run_all.sh"
  ),
  project_config_sha256 = file.path(
    project_root,
    "repro",
    "config",
    "project.toml"
  ),
  input_manifest_sha256 = file.path(
    project_root,
    "repro",
    "manifests",
    "input_sources.tsv"
  )
)
required <- c(historical_path, production_path, dependency_paths)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing comparison target(s): ", paste(missing, collapse = ", "))
}

modes <- c(
  "historical_default",
  "historical_serial",
  "explicit_seed_default",
  "explicit_seed_serial"
)
probe_paths <- unlist(lapply(modes, function(mode) {
  file.path(audit_root, paste0(mode, "_repeat", 1:2, ".csv"))
}))
metadata_paths <- sub("\\.csv$", ".metadata.tsv", probe_paths)
missing <- c(probe_paths, metadata_paths)[
  !file.exists(c(probe_paths, metadata_paths))
]
if (length(missing) > 0L) {
  stop("Missing RNG probe result(s): ", paste(missing, collapse = ", "))
}

read_result <- function(path) {
  result <- fread(path)
  setorder(result, ID)
  result
}

sha256_file <- function(path) {
  unname(tools::sha256sum(path))
}
expected_dependency_hashes <- vapply(
  dependency_paths,
  sha256_file,
  character(1)
)

historical <- read_result(historical_path)
production <- read_result(production_path)
numeric_columns <- c(
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue"
)

compare_results <- function(left, right) {
  if (!identical(left$ID, right$ID)) {
    stop("PTM-SEA ID order/set differs.")
  }
  deltas <- vapply(
    numeric_columns,
    function(column) {
      max(abs(left[[column]] - right[[column]]), na.rm = TRUE)
    },
    numeric(1)
  )
  list(
    exact_numeric = all(deltas == 0),
    max_numeric_delta = max(deltas),
    max_NES_delta = deltas[["NES"]],
    max_p_adjust_delta = deltas[["p.adjust"]],
    significant_set_exact = identical(
      left$p.adjust < 0.05,
      right$p.adjust < 0.05
    ),
    significant_count_left = sum(left$p.adjust < 0.05, na.rm = TRUE),
    significant_count_right = sum(right$p.adjust < 0.05, na.rm = TRUE)
  )
}

rows <- rbindlist(lapply(modes, function(mode) {
  first_path <- file.path(
    audit_root,
    paste0(mode, "_repeat1.csv")
  )
  second_path <- file.path(
    audit_root,
    paste0(mode, "_repeat2.csv")
  )
  first <- read_result(first_path)
  second <- read_result(second_path)
  first_metadata <- fread(sub("\\.csv$", ".metadata.tsv", first_path))
  second_metadata <- fread(sub("\\.csv$", ".metadata.tsv", second_path))
  if (nrow(first_metadata) != 1L || nrow(second_metadata) != 1L) {
    stop(mode, ": probe metadata must contain exactly one row per repeat.")
  }
  required_metadata <- c(
    "mode",
    "replicate",
    "external_seed",
    "GSEA_seed_argument",
    "registered_backend",
    "gmt_sha256",
    "protein_sha256",
    "phospho_sha256",
    "differential_sha256",
    "worker_script_sha256",
    "production_script_sha256",
    "orchestrator_script_sha256",
    "summarizer_script_sha256",
    "run_all_script_sha256",
    "project_config_sha256",
    "input_manifest_sha256",
    "result_sha256",
    "clusterProfiler_version",
    "fgsea_version",
    "BiocParallel_version",
    "R_version"
  )
  if (!all(required_metadata %in% names(first_metadata))) {
    stop(mode, ": probe metadata is incomplete.")
  }
  if (!all(required_metadata %in% names(second_metadata))) {
    stop(mode, ": second-repeat probe metadata is incomplete.")
  }
  metadata_match_columns <- setdiff(required_metadata, c("replicate", "result_sha256"))
  if (!identical(
    first_metadata[, ..metadata_match_columns],
    second_metadata[, ..metadata_match_columns]
  )) {
    stop(mode, ": probe metadata differs across clean-process repeats.")
  }
  if (
    first_metadata$mode != mode ||
      second_metadata$mode != mode ||
      first_metadata$replicate != 1L ||
      second_metadata$replicate != 2L ||
      first_metadata$external_seed != 42L ||
      second_metadata$external_seed != 42L
  ) {
    stop(mode, ": probe identity or external seed is invalid.")
  }
  expected_internal_seed <- startsWith(mode, "explicit_seed_")
  if (
    !identical(first_metadata$GSEA_seed_argument, expected_internal_seed) ||
      !identical(second_metadata$GSEA_seed_argument, expected_internal_seed)
  ) {
    stop(mode, ": GSEA seed argument does not match mode.")
  }
  if (
    endsWith(mode, "_serial") &&
      (
        first_metadata$registered_backend != "SerialParam" ||
          second_metadata$registered_backend != "SerialParam"
      )
  ) {
    stop(mode, ": serial mode did not register SerialParam.")
  }
  if (
    first_metadata$result_sha256 != sha256_file(first_path) ||
      second_metadata$result_sha256 != sha256_file(second_path)
  ) {
    stop(mode, ": result checksum does not match metadata.")
  }
  for (hash_column in names(expected_dependency_hashes)) {
    expected_hash <- expected_dependency_hashes[[hash_column]]
    if (
      first_metadata[[hash_column]] != expected_hash ||
        second_metadata[[hash_column]] != expected_hash
    ) {
      stop(
        mode,
        ": ",
        hash_column,
        " does not match the current dependency."
      )
    }
  }
  repeat_comparison <- compare_results(first, second)
  historical_comparison <- compare_results(first, historical)
  production_comparison <- compare_results(first, production)
  data.table(
    mode = mode,
    repeat_file_sha256_exact =
      first_metadata$result_sha256 == second_metadata$result_sha256,
    repeat_exact_numeric = repeat_comparison$exact_numeric,
    repeat_max_numeric_delta = repeat_comparison$max_numeric_delta,
    repeat_significant_set_exact = repeat_comparison$significant_set_exact,
    historical_exact_numeric = historical_comparison$exact_numeric,
    historical_max_NES_delta = historical_comparison$max_NES_delta,
    historical_max_p_adjust_delta =
      historical_comparison$max_p_adjust_delta,
    historical_significant_set_exact =
      historical_comparison$significant_set_exact,
    probe_significant_count =
      historical_comparison$significant_count_left,
    historical_significant_count =
      historical_comparison$significant_count_right,
    production_exact_numeric = production_comparison$exact_numeric,
    production_max_NES_delta = production_comparison$max_NES_delta,
    production_max_p_adjust_delta =
      production_comparison$max_p_adjust_delta,
    production_significant_set_exact =
      production_comparison$significant_set_exact,
    production_significant_count =
      production_comparison$significant_count_right
  )
}))

historical_modes <- rows[
  mode %in% c("historical_default", "historical_serial")
]
if (
  any(!rows$repeat_exact_numeric) ||
    any(!rows$repeat_file_sha256_exact) ||
    any(!rows$repeat_significant_set_exact) ||
    any(!historical_modes$historical_significant_set_exact) ||
    any(!historical_modes$production_significant_set_exact) ||
    any(historical_modes$historical_max_NES_delta > 0.001) ||
    any(historical_modes$historical_max_p_adjust_delta > 0.005)
) {
  stop("PTM-SEA RNG-mode audit failed its repeatability/concordance contract.")
}

fwrite(
  rows,
  file.path(
    project_root,
    "repo",
    "audit",
    "tissue_ptmsea_rng_mode_comparison.tsv"
  ),
  sep = "\t"
)
print(rows)
