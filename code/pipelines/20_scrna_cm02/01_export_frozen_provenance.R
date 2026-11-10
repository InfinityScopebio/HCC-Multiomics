#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

arguments <- commandArgs(trailingOnly = FALSE)
this_script <- normalizePath(
  sub("^--file=", "", arguments[grepl("^--file=", arguments)]),
  mustWork = TRUE
)
mode_arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(mode_arguments) >= 1L) mode_arguments[[1L]] else "audit"
if (!mode %in% c("audit", "curated", "full")) {
  stop("Mode must be one of: audit, curated, full.")
}
source(file.path(dirname(this_script), "config.R"))
config <- ensure_output_dirs(scrna_cm02_config(this_script))
set_recorded_seed(config, "runtime", "R_session")
write_seed_registry(config)

runtime_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "data.table", "NMF", "harmony",
  "psych", "igraph", "RANN", "future", "foreach"
)
runtime <- data.table(
  component = c("R", runtime_packages),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(
      runtime_packages,
      function(package_name) {
        if (requireNamespace(package_name, quietly = TRUE)) {
          as.character(packageVersion(package_name))
        } else {
          "not_installed"
        }
      },
      character(1)
    )
  )
)
fwrite(runtime, file.path(config$audit_root, "scrna_runtime_R.tsv"), sep = "\t")

curation_manifest <- file.path(
  config$project_root,
  "config",
  "manifests",
  "scrna_frozen_curation.tsv"
)
raw_and_rule_inputs <- if (mode == "curated") {
  character()
} else {
  c(
    unlist(config$raw_paths, use.names = FALSE),
    list.files(
      config$pipeline_dir,
      pattern = "^annotation_.*_raw_derived_overlay\\.csv$",
      full.names = TRUE
    ),
    file.path(
      config$annotation_root,
      c(
        "annotation.csv", "annotation_t.csv", "annotation_mp.csv",
        "annotation_fib.csv", "annotation_b.csv", "annotation_hep.csv"
      )
    )
  )
}
inventory <- data.table(
  path = c(
    raw_and_rule_inputs,
    curation_manifest,
    unlist(config$curated_atlas_paths, use.names = FALSE)
  ),
  input_role = c(
    rep("raw_reconstruction_or_explicit_rule", length(raw_and_rule_inputs)),
    "frozen_curation_manifest",
    rep(
      "frozen_expert_curation_production",
      length(config$curated_atlas_paths)
    )
  )
)
inventory[
  ,
  exists := file.exists(path)
]
inventory[
  ,
  production_input_allowed := vapply(
    path,
    function(current_path) {
      if (!file.exists(current_path)) {
        return(FALSE)
      }
      tryCatch(
        {
          assert_production_input(config, current_path)
          TRUE
        },
        error = function(error) FALSE
      )
    },
    logical(1)
  )
]
inventory[
  ,
  status := ifelse(exists & production_input_allowed, "PASS", "FAIL")
]
fwrite(
  inventory,
  file.path(config$audit_root, "scrna_production_input_paths.tsv"),
  sep = "\t"
)
if (any(inventory$status != "PASS")) {
  stop("Missing or prohibited production inputs.")
}

message(
  "PASS: runtime, seed plan, and checksum-locked frozen expert-curation ",
  "production inputs audited in ", mode, " mode; ",
  if (mode == "curated") {
    "raw scRNA paths were intentionally outside this production dependency; "
  } else {
    "raw-audit inputs were also checked; "
  },
  "strict seed audit remains unpublished until script 09 completes."
)
