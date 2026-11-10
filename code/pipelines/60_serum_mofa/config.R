resolve_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg) != 1L) {
    stop("Run pipeline scripts with Rscript.")
  }

  script_path <- normalizePath(file_arg[[1]], mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = TRUE)
}

source(file.path(
  resolve_project_root(),
  "code",
  "lib",
  "R",
  "deterministic_rng.R"
))

load_pipeline_config <- function() {
  project_root <- resolve_project_root()
  canonical_result_dir <- file.path(
    project_root,
    "results",
    "60_serum_mofa"
  )
  analysis_source <- Sys.getenv("MOFA_ANALYSIS_SOURCE", "article")
  if (!analysis_source %in% c("article", "raw_derived")) {
    stop("MOFA_ANALYSIS_SOURCE must be article or raw_derived.")
  }
  result_dir <- if (analysis_source == "raw_derived") {
    file.path(canonical_result_dir, "raw_to_figure")
  } else {
    canonical_result_dir
  }
  audit_dir <- if (analysis_source == "raw_derived") {
    file.path(result_dir, "audit")
  } else {
    file.path(project_root, "audit")
  }
  model_dir <- file.path(result_dir, "models")
  fast_model_dir <- file.path(model_dir, "k_scan_fast1000")
  final_model_dir <- file.path(model_dir, "final_slow2000")
  replay_model_dir <- file.path(model_dir, "repeatability")
  validation_dir <- file.path(
    project_root,
    "data",
    "validation",
    "frozen_targets",
    "MOFA"
  )
  validation_step49 <- file.path(validation_dir, "CURRENT_step49")
  historical_boundary_dir <- file.path(
    project_root,
    "data",
    "curation",
    "serum_mofa"
  )
  historical_boundary_manifest <- file.path(
    project_root,
    "config",
    "manifests",
    "serum_mofa_historical_boundary.tsv"
  )

  for (path in c(
    result_dir,
    audit_dir,
    model_dir,
    fast_model_dir,
    final_model_dir,
    replay_model_dir,
    file.path(result_dir, "logs")
  )) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  list(
    project_root = project_root,
    pipeline_dir = file.path(project_root, "code", "pipelines", "60_serum_mofa"),
    result_dir = result_dir,
    canonical_result_dir = canonical_result_dir,
    analysis_source = analysis_source,
    is_raw_derived_analysis = analysis_source == "raw_derived",
    audit_dir = audit_dir,
    model_dir = model_dir,
    fast_model_dir = fast_model_dir,
    final_model_dir = final_model_dir,
    replay_model_dir = replay_model_dir,
    log_dir = file.path(result_dir, "logs"),
    raw_input_path = file.path(
      result_dir,
      "incdec_input_rebuilt.rds"
    ),
    article_input_path = file.path(historical_boundary_dir, "incdec_input.rds"),
    article_k_scan_summary_path = file.path(
      historical_boundary_dir,
      "incdec_Kselect_summary.csv"
    ),
    article_model_path = file.path(
      historical_boundary_dir,
      "incdec_mofa_K15.hdf5"
    ),
    historical_boundary_dir = historical_boundary_dir,
    historical_boundary_manifest = historical_boundary_manifest,
    final_model_path = file.path(
      final_model_dir,
      "incdec_mofa_K15_seed1_slow2000.hdf5"
    ),
    final_model_metadata_path = file.path(
      final_model_dir,
      "incdec_mofa_K15_seed1_slow2000.metadata.tsv"
    ),
    raw_k_scan_summary_path = file.path(
      result_dir,
      "k_scan_fast1000_summary_raw_derived.tsv"
    ),
    analysis_model_path = if (analysis_source == "raw_derived") {
      file.path(final_model_dir, "incdec_mofa_K15_seed1_slow2000.hdf5")
    } else {
      file.path(historical_boundary_dir, "incdec_mofa_K15.hdf5")
    },
    analysis_k_scan_summary_path = if (analysis_source == "raw_derived") {
      file.path(
        result_dir,
        "k_scan_fast1000_summary_raw_derived.tsv"
      )
    } else {
      file.path(historical_boundary_dir, "incdec_Kselect_summary.csv")
    },
    m0_dir = file.path(
      project_root,
      "results",
      "00_preprocess",
      "multiomics",
      "dat"
    ),
    tcga_raw_dir = file.path(project_root, "data", "public", "raw", "TCGA-LIHC"),
    local_kegg_compound_pathway = file.path(
      project_root,
      "data",
  "reference",
  "KEGG",
      "kegg_cpd_pathway.tsv"
    ),
    python_path = Sys.getenv(
      "MOFA_PYTHON",
      "/Users/pyramidsnail/miniconda3/bin/python"
    ),
    master_seed = validate_r_seed(
      Sys.getenv("LIUF_MASTER_SEED", "20260807"),
      "LIUF_MASTER_SEED"
    ),
    k_grid = 2:20,
    k_scan_seeds = 1:6,
    validation_dir = validation_dir,
    validation = list(
      input = file.path(validation_step49, "incdec_input.rds"),
      final_model = file.path(validation_step49, "incdec_mofa_K15.hdf5"),
      fast_k_scan_summary = file.path(
        validation_dir,
        "05_血清MOFA",
        "incdec_Kselect_summary.csv"
      ),
      slow_k_scan_cache = file.path(
        validation_step49,
        "incdec_Kselect_runs_slow2000.rds"
      ),
      kegg_compound_pathway = file.path(
        validation_step49,
        "kegg_cpd_pathway.tsv"
      ),
      axis_store = file.path(validation_step49, "incdec_K15_all_axes.rds"),
      pathways = file.path(
        validation_step49,
        "incdec_K15_allfactors_pathways.csv"
      ),
      survival_stats = file.path(
        validation_dir,
        "tcga_compartment",
        "Fig_incdec_K15_allfactors_TCGA_survival_stats.csv"
      ),
      model_comparison = file.path(
        validation_step49,
        "Fig1_model_comparison.csv"
      ),
      performance = file.path(
        validation_step49,
        "Fig2_performance_metrics.csv"
      ),
      permutation_b1000 = file.path(
        validation_step49,
        "MLR7_permutation_importance_B1000.csv"
      ),
      supplemental_factor_selection = file.path(
        validation_step49,
        "incdec_K15_dual_significant_survival_factors.csv"
      ),
      supplemental_loadings = file.path(
        validation_step49,
        "incdec_K15_dual_significant_raw_top_loadings.csv"
      ),
      supplemental_pathways = file.path(
        validation_step49,
        "incdec_K15_selected_factor_pathways_full.csv"
      ),
      supplemental_ranked_profiles = file.path(
        validation_step49,
        "incdec_K15_selected_factor_ranked_profiles.csv"
      )
    )
  )
}

require_local_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      "Required local inputs are missing; this pipeline never downloads data: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(paths)
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

parse_named_arguments <- function(arguments) {
  parsed <- list()
  for (argument in arguments) {
    if (!grepl("^--[^=]+=", argument)) {
      stop("Arguments must use --name=value syntax: ", argument)
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    parsed[[pieces[[1]]]] <- paste(pieces[-1], collapse = "=")
  }
  parsed
}

fast_model_path <- function(cfg, k, seed) {
  file.path(
    cfg$fast_model_dir,
    sprintf("incdec_mofa_K%02d_seed%02d_fast1000.hdf5", k, seed)
  )
}

fast_model_metadata_path <- function(cfg, k, seed) {
  sub(
    "\\.hdf5$",
    ".metadata.tsv",
    fast_model_path(cfg, k, seed)
  )
}

replay_model_path <- function(cfg, profile, k, seed) {
  file.path(
    cfg$replay_model_dir,
    sprintf("incdec_mofa_%s_K%02d_seed%02d.hdf5", profile, k, seed)
  )
}

replay_model_metadata_path <- function(cfg, profile, k, seed) {
  sub(
    "\\.hdf5$",
    ".metadata.tsv",
    replay_model_path(cfg, profile, k, seed)
  )
}

semantic_model_hash <- function(model) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for semantic model hashes.")
  }
  variance <- MOFA2::get_variance_explained(model)
  payload <- list(
    factors = MOFA2::get_factors(model),
    weights = MOFA2::get_weights(model),
    r2_per_factor = variance$r2_per_factor,
    r2_total = variance$r2_total,
    samples = rownames(MOFA2::get_factors(model)[[1]]),
    features = lapply(MOFA2::get_weights(model), rownames),
    data_options = model@data_options,
    model_options = model@model_options,
    training_options = model@training_options
  )
  digest::digest(payload, algo = "sha256", serialize = TRUE)
}

read_model_metadata <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  data.table::fread(path)
}

validate_production_model <- function(
    model_path,
    metadata_path,
    input_sha256,
    profile,
    k,
    seed,
    maxiter,
    convergence_mode) {
  if (!file.exists(model_path) || !file.exists(metadata_path)) {
    return(FALSE)
  }
  metadata <- read_model_metadata(metadata_path)
  if (is.null(metadata) || nrow(metadata) != 1L) {
    return(FALSE)
  }
  identical(metadata$status, "complete") &&
    identical(metadata$input_sha256, input_sha256) &&
    identical(metadata$profile, profile) &&
    identical(as.integer(metadata$K), as.integer(k)) &&
    identical(as.integer(metadata$seed), as.integer(seed)) &&
    identical(as.integer(metadata$maxiter), as.integer(maxiter)) &&
    identical(metadata$convergence_mode, convergence_mode) &&
    identical(
      as.integer(metadata$master_seed),
      as.integer(Sys.getenv("LIUF_MASTER_SEED", "20260807"))
    ) &&
    identical(metadata$file_sha256, sha256_file(model_path))
}

require_production_model <- function(cfg) {
  require_local_files(c(
    cfg$raw_input_path,
    cfg$final_model_path,
    cfg$final_model_metadata_path
  ))
  valid <- validate_production_model(
    model_path = cfg$final_model_path,
    metadata_path = cfg$final_model_metadata_path,
    input_sha256 = sha256_file(cfg$raw_input_path),
    profile = "final_slow2000",
    k = 15L,
    seed = 1L,
    maxiter = 2000L,
    convergence_mode = "slow"
  )
  if (!valid) {
    stop(
      "The production K15 model is absent or has invalid raw-derived provenance. ",
      "Frozen validation models are never used as fallback."
    )
  }
  invisible(cfg$final_model_path)
}

require_historical_mofa_boundary <- function(cfg) {
  expected_paths <- c(
    input = cfg$article_input_path,
    k_scan_summary = cfg$article_k_scan_summary_path,
    k15_model = cfg$article_model_path
  )
  require_local_files(c(cfg$historical_boundary_manifest, expected_paths))
  manifest <- data.table::fread(
    cfg$historical_boundary_manifest,
    encoding = "UTF-8"
  )
  if (!identical(sort(manifest$role), sort(names(expected_paths)))) {
    stop("Historical MOFA boundary manifest roles have changed.")
  }
  for (role in names(expected_paths)) {
    role_name <- role
    record <- manifest[manifest$role == role_name]
    if (nrow(record) != 1L) {
      stop("Historical MOFA boundary requires one row for role: ", role)
    }
    manifest_path <- normalizePath(
      file.path(cfg$project_root, record$materialized_path),
      mustWork = TRUE
    )
    configured_path <- normalizePath(expected_paths[[role]], mustWork = TRUE)
    if (manifest_path != configured_path ||
        as.numeric(file.info(configured_path)$size) != as.numeric(record$bytes) ||
        sha256_file(configured_path) != record$sha256) {
      stop("Historical MOFA boundary contract failed for role: ", role)
    }
    if (startsWith(configured_path, normalizePath(cfg$validation_dir))) {
      stop("Article inputs must not resolve inside validation/frozen_targets.")
    }
  }
  invisible(expected_paths)
}

require_analysis_mofa_boundary <- function(cfg) {
  if (isTRUE(cfg$is_raw_derived_analysis)) {
    require_production_model(cfg)
    require_local_files(c(
      cfg$raw_input_path,
      cfg$analysis_model_path,
      cfg$analysis_k_scan_summary_path
    ))
    return(invisible(c(
      input = cfg$raw_input_path,
      k_scan_summary = cfg$analysis_k_scan_summary_path,
      k15_model = cfg$analysis_model_path
    )))
  }
  require_historical_mofa_boundary(cfg)
}
