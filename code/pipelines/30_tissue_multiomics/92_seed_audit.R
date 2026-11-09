#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_path <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

seed_audit <- data.table(
  pipeline = "tissue_multiomics",
  step = c(
    "pipeline_default",
    "seed_derivation_helper",
    "PTMSEA_null_sampling",
    "PTMSEA_rng_mode_worker",
    "PTMSEA_rng_mode_orchestrator",
    "PTMSEA_rng_mode_summarizer",
    "cross_layer_bootstrap"
  ),
  script_path = c(
    "code/pipelines/30_tissue_multiomics/common.R",
    "code/lib/R/deterministic_rng.R",
    "code/pipelines/30_tissue_multiomics/20_ptmsea.R",
    "code/pipelines/30_tissue_multiomics/22_probe_ptmsea_rng_mode.R",
    "code/pipelines/30_tissue_multiomics/22_run_ptmsea_rng_probe.R",
    "code/pipelines/30_tissue_multiomics/23_summarize_ptmsea_rng_modes.R",
    "code/pipelines/30_tissue_multiomics/30_cross_layer.R"
  ),
  language = "R",
  algorithm = c(
    "R base RNG fallback",
    "deterministic modulo-2^31 seed derivation",
    "PTM-SEA null sampling",
    "PTM-SEA RNG-mode repeatability worker",
    "clean-process PTM-SEA probe orchestration",
    "PTM-SEA RNG-mode comparison and fail-closed validation",
    "nonparametric bootstrap"
  ),
  seed = c(
    master_seed,
    master_seed,
    42L,
    42L,
    42L,
    42L,
    42L
  ),
  seed_source = c(
    "project_default",
    "project_default",
    "historical",
    "historical",
    "historical",
    "historical",
    "historical"
  ),
  rng_kind = c(
    "Mersenne-Twister;Inversion;Rejection",
    "none",
    rep("Mersenne-Twister;Inversion;Rejection", 5L)
  ),
  status = c("fixed", "deterministic", rep("fixed", 5L)),
  notes = c(
    "Used before any step-specific seed is set.",
    paste(
      "Shared audited helper wraps deterministic offsets into the valid",
      "non-negative R seed domain."
    ),
    paste(
      "Preserved to reproduce the accepted PTM-SEA result;",
      "external set.seed(42), clusterProfiler seed=FALSE,",
      "BiocParallel SerialParam, and an exact in-run repeat."
    ),
    paste(
      "Audit-only clean-process comparison of historical/default and",
      "explicit/internal GSEA seed modes; every mode is repeated twice."
    ),
    paste(
      "Deletes stale probe and comparison files, then launches each repeat",
      "in an independent R process."
    ),
    paste(
      "Validates result and dependency checksums, clean-process equality,",
      "historical concordance, and the production result contract."
    ),
    "Preserved to reproduce the accepted bootstrap estimates."
  )
)

fwrite(
  seed_audit,
  file.path(audit_root, "tissue_seed_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)
