#!/usr/bin/env Rscript

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
pipeline_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(
  file.path(pipeline_dir, "..", "..", ".."),
  mustWork = TRUE
)
audit_root <- file.path(project_root, "repo", "audit", "ptmsea_rng_modes")
log_root <- file.path(
  project_root,
  "repo",
  "results",
  "30_tissue_multiomics",
  "logs",
  "ptmsea_rng_modes"
)
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
comparison_path <- file.path(
  project_root,
  "repo",
  "audit",
  "tissue_ptmsea_rng_mode_comparison.tsv"
)
unlink(comparison_path, force = TRUE)

modes <- c(
  "historical_default",
  "historical_serial",
  "explicit_seed_default",
  "explicit_seed_serial"
)
expected_paths <- unlist(lapply(modes, function(mode) {
  unlist(lapply(1:2, function(replicate_id) {
    stem <- sprintf("%s_repeat%d", mode, replicate_id)
    file.path(audit_root, paste0(stem, c(".csv", ".metadata.tsv")))
  }))
}))
unlink(expected_paths, force = TRUE)

rscript <- Sys.getenv("LIUF_RSCRIPT", unset = "")
if (!nzchar(rscript)) {
  rscript <- file.path(R.home("bin"), "Rscript")
}
resolved_rscript <- Sys.which(rscript)
if (!nzchar(resolved_rscript)) {
  if (!file.exists(rscript)) {
    stop("Rscript runtime is unavailable: ", rscript)
  }
  resolved_rscript <- normalizePath(rscript, mustWork = TRUE)
}
worker <- file.path(pipeline_dir, "22_probe_ptmsea_rng_mode.R")

for (mode in modes) {
  for (replicate_id in 1:2) {
    stem <- sprintf("%s_repeat%d", mode, replicate_id)
    log_path <- file.path(log_root, paste0(stem, ".log"))
    status <- system2(
      resolved_rscript,
      c(
        worker,
        paste0("--mode=", mode),
        paste0("--replicate=", replicate_id)
      ),
      stdout = log_path,
      stderr = log_path
    )
    if (!identical(status, 0L)) {
      stop("PTM-SEA RNG probe failed for ", stem, "; see ", log_path)
    }
  }
}

missing <- expected_paths[!file.exists(expected_paths)]
if (length(missing) > 0L) {
  stop("PTM-SEA RNG probe omitted expected output(s): ", paste(missing, collapse = ", "))
}
cat("Completed 4 PTM-SEA RNG modes x 2 clean-process repeats.\n")
