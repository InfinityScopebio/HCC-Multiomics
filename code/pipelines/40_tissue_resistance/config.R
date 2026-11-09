script_file_path <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- sub(
    "^--file=",
    "",
    arguments[grepl("^--file=", arguments)]
  )
  if (length(file_argument) != 1L) {
    stop("Run this file with Rscript.")
  }
  normalizePath(file_argument, mustWork = TRUE)
}

tissue_resistance_config <- function(script_path = script_file_path()) {
  pipeline_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
  project_root <- normalizePath(
    file.path(pipeline_dir, "..", "..", ".."),
    mustWork = TRUE
  )
  result_root <- file.path(
    project_root,
    "results",
    "40_tissue_resistance"
  )
  curation_root <- file.path(
    project_root,
    "data",
    "curation",
    "resistance"
  )
  depmap_root <- file.path(
    project_root,
    "data",
    "public",
    "raw",
    "DepMap"
  )
  list(
    project_root = project_root,
    pipeline_dir = pipeline_dir,
    result_root = result_root,
    audit_root = file.path(project_root, "audit"),
    curation_root = curation_root,
    manifest_path = file.path(
      project_root,
      "config",
      "manifests",
      "depmap_prism_historical_boundary.tsv"
    ),
    inputs = list(
      model = file.path(depmap_root, "Model.csv"),
      expression = file.path(
        depmap_root,
        "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
      ),
      prism_auc = file.path(
        depmap_root,
        "Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv"
      ),
      signatures = file.path(curation_root, "HepC_signatures_main.gmt")
    ),
    references = list(
      score_table = file.path(
        project_root,
        "resistance_landscape",
        "R3d_cell_line_scores.csv"
      ),
      tki_profile = file.path(
        project_root,
        "resistance_landscape",
        "R3d_profile_liver_biliary_TKI.csv"
      ),
      chemo_profile = file.path(
        project_root,
        "resistance_landscape",
        "R3d_profile_liver_biliary_TACE_HAIC_chemo.csv"
      ),
      panel = file.path(
        project_root,
        "resistance_landscape",
        "Fig_R3d_PanelA.pdf"
      )
    ),
    clusters = paste0("Hep_C", seq_len(7L)),
    composites = list(
      TKI = c("SORAFENIB", "LENVATINIB", "REGORAFENIB", "CABOZANTINIB"),
      TACE_HAIC_chemo = c(
        "DOXORUBICIN", "EPIRUBICIN", "IDARUBICIN", "CISPLATIN",
        "OXALIPLATIN", "FLUOROURACIL", "CAPECITABINE", "GEMCITABINE",
        "RALTITREXED"
      )
    )
  )
}

ensure_tissue_resistance_dirs <- function(config) {
  dir.create(config$result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$audit_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$curation_root, recursive = TRUE, showWarnings = FALSE)
  invisible(config)
}

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    c("-a", "256", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Failed to hash file: ", path, "\n", paste(output, collapse = "\n"))
  }
  sub("\\s+.*$", "", output[[1L]])
}
