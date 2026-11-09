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

scrna_cm02_config <- function(script_path = script_file_path()) {
  pipeline_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
  project_root <- normalizePath(
    file.path(pipeline_dir, "..", "..", ".."),
    mustWork = TRUE
  )
  result_root <- file.path(
    project_root,
    "results",
    "20_scrna_cm02"
  )
  audit_root <- file.path(project_root, "audit")
  curation_root <- file.path(
    project_root,
    "data",
    "curation",
    "scRNA"
  )
  downstream_curation_root <- file.path(
    curation_root,
    "cm02_downstream"
  )
  downstream_source <- Sys.getenv("FIG3_DOWNSTREAM_SOURCE", "historical")
  if (!downstream_source %in% c("historical", "raw_derived")) {
    stop("FIG3_DOWNSTREAM_SOURCE must be historical or raw_derived.")
  }
  raw_to_figure3_root <- file.path(result_root, "raw_to_figure3")
  raw_music_root <- file.path(raw_to_figure3_root, "music")
  raw_trajectory_root <- file.path(raw_to_figure3_root, "trajectory")

  master_seed <- suppressWarnings(
    as.integer(Sys.getenv("LIUF_MASTER_SEED", unset = "20260807"))
  )
  if (is.na(master_seed)) {
    stop("LIUF_MASTER_SEED must be an integer.")
  }

  historical_downstream_paths <- list(
    nmf = file.path(downstream_curation_root, "NMF_K.rds"),
    tcga_proportions = file.path(downstream_curation_root, "prop_TCGA.csv"),
    bulk_frequency = file.path(downstream_curation_root, "combined_freq.csv"),
    bulk_pseudotime = file.path(
      downstream_curation_root,
      "umap_CM02_pseudotime.csv"
    ),
    tcga_clinical = file.path(
      project_root,
      "data",
      "public",
      "raw",
      "TCGA-LIHC",
      "LIHC_clinicalMatrix.tsv"
    )
  )
  raw_downstream_paths <- list(
    nmf = file.path(result_root, "covarnet", "NMF_K5_frozen_curated.rds"),
    tcga_proportions = file.path(raw_music_root, "prop_TCGA.csv"),
    bulk_frequency = file.path(raw_trajectory_root, "combined_freq.csv"),
    bulk_pseudotime = file.path(
      raw_trajectory_root,
      "umap_CM02_pseudotime.csv"
    ),
    tcga_clinical = historical_downstream_paths$tcga_clinical
  )

  list(
    master_seed = master_seed,
    project_root = project_root,
    pipeline_dir = pipeline_dir,
    result_root = result_root,
    audit_root = audit_root,
    seed_audit_path = file.path(audit_root, "scrna_seed_audit.tsv"),
    curation_audit_path = file.path(
      audit_root,
      "scrna_frozen_curation_objects.tsv"
    ),
    project_config_path = file.path(
      project_root,
      "config",
      "project.toml"
    ),
    materialized_input_manifest_path = file.path(
      project_root,
      "manifests",
      "input_inventory_materialized.tsv"
    ),
    validation_root = file.path(
      project_root,
      "data",
      "validation",
      "frozen_targets",
      "scRNA"
    ),
    historical_downstream_manifest_path = file.path(
      project_root,
      "config",
      "manifests",
      "cm02_historical_downstream.tsv"
    ),
    supplemental2_rank_manifest_path = file.path(
      project_root,
      "config",
      "manifests",
      "supplemental2_rank_survey_historical.tsv"
    ),
    downstream_curation_root = downstream_curation_root,
    downstream_source = downstream_source,
    raw_to_figure3_root = raw_to_figure3_root,
    raw_music_root = raw_music_root,
    raw_trajectory_root = raw_trajectory_root,
    figure3_output_root = if (downstream_source == "raw_derived") {
      file.path(raw_to_figure3_root, "figure3_downstream")
    } else {
      file.path(result_root, "figure3_downstream")
    },
    historical_downstream_paths = c(
      historical_downstream_paths,
      list(
      rank_metrics = file.path(
        downstream_curation_root,
        "rank_survey_filtered_metrics.csv"
      ),
      rank_matrix = file.path(
        downstream_curation_root,
        "matnorm_filtered_21x39.csv"
      )
      )
    ),
    raw_downstream_paths = raw_downstream_paths,
    analysis_downstream_paths = if (downstream_source == "raw_derived") {
      raw_downstream_paths
    } else {
      historical_downstream_paths
    },
    internal_samples = c("HCC1", "HCC2", "HCC3", "NC1", "NC2", "NC3"),
    internal_suffix = c(
      HCC1 = 4L,
      HCC2 = 5L,
      HCC3 = 6L,
      NC1 = 7L,
      NC2 = 8L,
      NC3 = 9L
    ),
    prohibited_samples = c("CC1", "CC2", "CC3"),
    group_levels = c("Healthy", "Cirrhosis", "Adjacent", "HCC"),
    raw_paths = list(
      internal_root = file.path(
        project_root,
        "data",
        "private",
        "raw",
        "data",
        "scRNAData"
      ),
      gse115469 = file.path(
        project_root,
        "data",
        "public",
        "raw",
        "GSE115469",
        "GSE115469_Data.csv.gz"
      ),
      gse136103 = file.path(
        project_root,
        "data",
        "public",
        "raw",
        "GSE136103",
        "GSE136103_extracted"
      ),
      gse149614_counts = file.path(
        project_root,
        "data",
        "public",
        "raw",
        "GSE149614",
        "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
      ),
      gse149614_metadata = file.path(
        project_root,
        "data",
        "public",
        "raw",
        "GSE149614",
        "GSE149614_HCC.metadata.updated.txt.gz"
      )
    ),
    annotation_root = file.path(
      project_root,
      "data",
      "reference",
      "annotations"
    ),
    curation_root = curation_root,
    curated_atlas_paths = list(
      final = file.path(curation_root, "combined_annotated_v5.rds"),
      hepatocyte = file.path(curation_root, "hep_annotated.rds"),
      t_cell = file.path(curation_root, "t_annotated.rds"),
      macrophage = file.path(curation_root, "mp_annotated.rds"),
      fibroblast = file.path(curation_root, "fib_annotated.rds"),
      b_cell = file.path(curation_root, "b_annotated.rds")
    ),
    atlas_paths = list(
      qc = file.path(result_root, "atlas_qc_raw_derived.rds"),
      uncurated = file.path(
        result_root,
        "atlas_uncurated_raw_derived.rds"
      ),
      integrated = file.path(result_root, "atlas_major_raw_derived.rds"),
      final = file.path(result_root, "atlas_final_raw_derived.rds"),
      epithelial = file.path(result_root, "epithelial_raw_derived.rds"),
      hepatocyte_initial = file.path(
        result_root,
        "hep_initial_raw_derived.rds"
      ),
      hepatocyte = file.path(result_root, "hep_raw_derived.rds"),
      t_cell = file.path(result_root, "t_raw_derived.rds"),
      macrophage = file.path(result_root, "mp_raw_derived.rds"),
      fibroblast = file.path(result_root, "fib_raw_derived.rds"),
      b_cell = file.path(result_root, "b_raw_derived.rds")
    ),
    covarnet_paths = list(
      frequency = file.path(
        result_root,
        "covarnet",
        "freq_inputs_frozen_curated.rds"
      ),
      nmf = file.path(
        result_root,
        "covarnet",
        "NMF_K5_frozen_curated.rds"
      ),
      network = file.path(
        result_root,
        "covarnet",
        "cm_network_frozen_curated.rds"
      )
    ),
    cm02_anchor_states = c(
      "Treg", "Pericyte", "IL1B+M_", "SPP1+M_", "HSC",
      "TREM2+M_", "CXCL10+M_", "Hep_C4"
    )
  )
}

seed_registry <- function(config) {
  make_row <- function(
      stage,
      step_name,
      algorithm,
      seed,
      library,
      source,
      reason,
      rng_kind = "Mersenne-Twister/Inversion/Rejection") {
    data.frame(
      stage = stage,
      step_name = step_name,
      algorithm = algorithm,
      seed = as.integer(seed),
      rng_kind = rng_kind,
      library = library,
      source = source,
      reason = reason,
      stringsAsFactors = FALSE
    )
  }
  rpca_stage <- function(stage) {
    rbind(
      make_row(
        stage, "RunPCA", "irlba PCA", 42L, "Seurat", "historical",
        "Seurat historical default seed.use=42"
      ),
      make_row(
        stage, "RPCAIntegration", "Seurat RPCA", 42L, "Seurat",
        "historical", "Historical RPCA run inherited Seurat seed 42"
      ),
      make_row(
        stage, "FindNeighbors", "Annoy KNN", config$master_seed,
        "Seurat", "default",
        "No seed argument; master seed set immediately before KNN"
      ),
      make_row(
        stage, "FindClusters", "Louvain modularity", 0L, "Seurat",
        "historical", "Seurat historical default random.seed=0"
      ),
      make_row(
        stage, "RunUMAP", "uwot UMAP", 42L, "uwot", "historical",
        "Seurat historical default seed.use=42"
      )
    )
  }
  harmony_stage <- function(stage) {
    rbind(
      make_row(
        stage, "RunPCA", "irlba PCA", 42L, "Seurat", "historical",
        "Seurat historical default seed.use=42"
      ),
      make_row(
        stage, "HarmonyIntegration", "Harmony", config$master_seed,
        "harmony", "default",
        "Historical call had no seed argument; master seed fixed before call"
      ),
      make_row(
        stage, "FindNeighbors", "Annoy KNN", config$master_seed,
        "Seurat", "default",
        "No seed argument; master seed set immediately before KNN"
      ),
      make_row(
        stage, "FindClusters", "Louvain modularity", 0L, "Seurat",
        "historical", "Seurat historical default random.seed=0"
      ),
      make_row(
        stage, "RunUMAP", "uwot UMAP", 42L, "uwot", "historical",
        "Seurat historical default seed.use=42"
      )
    )
  }
  rbind(
    make_row(
      "runtime", "R_session", "R base RNG", config$master_seed, "base",
      "default", "Project-wide master seed"
    ),
    make_row(
      "raw_ingest", "matrix_construction", "deterministic sparse import",
      config$master_seed, "Matrix/Seurat", "default",
      "No sampling; master seed set before object construction"
    ),
    rpca_stage("atlas"),
    rpca_stage("epithelial"),
    rpca_stage("t_cell"),
    harmony_stage("macrophage"),
    harmony_stage("fibroblast"),
    harmony_stage("b_cell"),
    rpca_stage("hepatocyte_initial"),
    rpca_stage("hepatocyte_clean"),
    make_row(
      "covarnet", "nsNMF_K5", "NMF nsNMF", 77L, "NMF", "historical",
      "Historical CoVarNet final K=5 nsNMF seed"
    ),
    make_row(
      "covarnet", "foreach_backend", "foreach doSEQ", 77L, "foreach",
      "historical", "Force sequential NMF runs; eliminate parallel RNG"
    ),
    make_row(
      "python", "python_runtime", "random+numpy", config$master_seed,
      "Python random/numpy", "default",
      "Project-wide master seed and PYTHONHASHSEED",
      rng_kind = "Python random MT19937/numpy MT19937/PYTHONHASHSEED"
    )
  )
}

set_recorded_seed <- function(config, stage, step_name) {
  registry <- seed_registry(config)
  selected <- registry[
    registry$stage == stage & registry$step_name == step_name,
    ,
    drop = FALSE
  ]
  if (nrow(selected) != 1L) {
    stop("Seed registry lookup is not unique: ", stage, "/", step_name)
  }
  if (
    selected$rng_kind[[1L]] ==
      "Mersenne-Twister/Inversion/Rejection"
  ) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  } else {
    stop(
      "Unsupported R RNG kind for ",
      stage,
      "/",
      step_name,
      ": ",
      selected$rng_kind[[1L]]
    )
  }
  set.seed(as.integer(selected$seed[[1L]]))
  invisible(as.integer(selected$seed[[1L]]))
}

seed_audit_dependency_table <- function(config) {
  dependency_paths <- c(
    list.files(
      config$pipeline_dir,
      pattern = "\\.(R|py|sh|csv)$",
      full.names = TRUE
    ),
    config$project_config_path,
    config$materialized_input_manifest_path
  )
  dependency_paths <- sort(unique(normalizePath(
    dependency_paths,
    mustWork = TRUE
  )))
  dependency_info <- file.info(dependency_paths)
  data.table::data.table(
    path = sub(
      paste0("^", config$project_root, .Platform$file.sep),
      "",
      dependency_paths
    ),
    sha256 = vapply(dependency_paths, sha256_file, character(1)),
    source_mtime_epoch = as.numeric(dependency_info$mtime),
    source_mtime_utc = format(
      dependency_info$mtime,
      "%Y-%m-%dT%H:%M:%OS6Z",
      tz = "UTC"
    )
  )
}

write_seed_registry <- function(config) {
  dependency_audit <- seed_audit_dependency_table(config)
  dependency_summary <- paste(
    paste0(
      basename(dependency_audit$path),
      "=",
      substr(dependency_audit$sha256, 1L, 12L)
    ),
    collapse = ";"
  )
  detailed_path <- file.path(
    config$audit_root,
    "scrna_seed_audit_detailed.tsv"
  )
  data.table::fwrite(
    seed_registry(config),
    detailed_path,
    sep = "\t"
  )
  detailed <- data.table::as.data.table(seed_registry(config))
  detailed[
    ,
    script_path := data.table::fcase(
      stage == "runtime",
      "code/pipelines/20_scrna_cm02/01_export_frozen_provenance.R",
      stage == "raw_ingest",
      "code/pipelines/20_scrna_cm02/03_rebuild_atlas_from_local_raw.R",
      stage %in% c(
        "atlas"
      ),
      "code/pipelines/20_scrna_cm02/03_rebuild_atlas_from_local_raw.R",
      stage %in% c(
        "epithelial", "t_cell", "macrophage", "fibroblast", "b_cell",
        "hepatocyte_initial", "hepatocyte_clean"
      ),
      "code/pipelines/20_scrna_cm02/04_rebuild_atlas_panel_data.R",
      stage == "covarnet",
      "code/pipelines/20_scrna_cm02/05_rebuild_covarnet_cm02.R",
      default = "code/pipelines/20_scrna_cm02/02_audit_local_raw_coverage.py"
    )
  ]
  aggregate <- detailed[
    ,
    .(
      pipeline = "20_scrna_cm02",
      step = paste(stage, step_name, sep = "/"),
      script_path = script_path,
      language = ifelse(grepl("\\.py$", script_path), "Python", "R"),
      algorithm = algorithm,
      seed = seed,
      seed_source = source,
      rng_kind = rng_kind,
      planned_status = "fixed",
      publication_gate = "none",
      notes = reason
    )
  ]
  aggregate[
    grepl(
      "/(FindNeighbors|FindClusters|RunUMAP|nsNMF_K5)$",
      step
    ),
    `:=`(
      planned_status = "fixed_or_blocked_after_execution",
      publication_gate = "09 exact/semantic duplicate completed",
      notes = paste(
        notes,
        "Strict status is intentionally unpublished until duplicate execution.",
        sep = " | "
      )
    )
  ]
  deterministic_scripts <- data.table::data.table(
    pipeline = "20_scrna_cm02",
    step = c(
      "raw_input_hash_and_cc_assertion",
      "validation_comparison",
      "panel_audit_reporting",
      "output_validation",
      "pipeline_environment"
      ,
      "seed_repeatability"
      ,
      "rng_dispatch"
    ),
    script_path = c(
      "code/pipelines/20_scrna_cm02/02_audit_local_raw_coverage.py",
      "code/pipelines/20_scrna_cm02/06_audit_cm02_downstream.R",
      "code/pipelines/20_scrna_cm02/07_build_panel_audit.py",
      "code/pipelines/20_scrna_cm02/08_validate_outputs.py",
      "code/pipelines/20_scrna_cm02/run_all.sh"
      ,
      "code/pipelines/20_scrna_cm02/09_test_seed_repeatability.R"
      ,
      "code/pipelines/20_scrna_cm02/config.R"
    ),
    language = c("Python", "R", "Python", "Python", "Shell", "R", "R"),
    algorithm = c(
      "SHA256 inventory and directory assertions",
      "deterministic raw-derived versus frozen comparison",
      "deterministic audit table rendering",
      "deterministic invariant checks",
      "LIUF_MASTER_SEED environment propagation"
      ,
      "Annoy/Louvain/UMAP/NMF same-input duplicate execution"
      ,
      "Per-stage RNGkind and set.seed dispatcher"
    ),
    seed = c(
      rep(config$master_seed, 5L),
      config$master_seed,
      config$master_seed
    ),
    seed_source = rep("default", 7L),
    rng_kind = c(
      "Python random/numpy/PYTHONHASHSEED",
      "none",
      "Python random/numpy/PYTHONHASHSEED",
      "Python random/numpy/PYTHONHASHSEED",
      "environment"
      ,
      "Mersenne-Twister/Inversion/Rejection"
      ,
      "Mersenne-Twister/Inversion/Rejection"
    ),
    planned_status = c(
      "deterministic", "deterministic", "deterministic", "deterministic",
      "fixed", "fixed_or_blocked_after_execution", "deterministic"
    ),
    publication_gate = c(
      rep("none", 5L),
      "all 09 repeatability rows completed",
      "none"
    ),
    notes = c(
      "No sampling; RNG initialized to prevent hidden library defaults",
      "No random call",
      "No sampling; RNG initialized to prevent hidden library defaults",
      "No sampling; RNG initialized to prevent hidden library defaults",
      paste0(
        "Default seed 20260807; caller may override LIUF_MASTER_SEED. ",
        "LIUF_THREADS controls OMP/OPENBLAS/MKL/VECLIB/NUMEXPR/BLIS/",
        "RCPP_PARALLEL; default 1"
      )
      ,
      "Strict audit is withheld until all duplicate executions finish"
      ,
      "Uses historical/default registry entry before every R random call"
    )
  )
  aggregate <- unique(rbind(aggregate, deterministic_scripts))
  aggregate[
    step == "rng_dispatch",
    notes := paste(
      notes,
      paste0(
        "Registry regenerated from source/config/materialized-manifest SHA256: ",
        dependency_summary
      ),
      sep = " | "
    )
  ]
  data.table::setcolorder(
    aggregate,
    c(
      "pipeline", "step", "script_path", "language", "algorithm", "seed",
      "seed_source", "rng_kind", "planned_status", "publication_gate",
      "notes"
    )
  )
  plan_path <- file.path(
    config$audit_root,
    "scrna_seed_audit_publication_plan.tsv"
  )
  data.table::fwrite(aggregate, plan_path, sep = "\t")
  strict_dependency_path <- file.path(
    config$audit_root,
    "scrna_seed_audit_dependencies.tsv"
  )
  if (file.exists(config$seed_audit_path)) {
    unlink(config$seed_audit_path)
  }
  if (file.exists(strict_dependency_path)) {
    unlink(strict_dependency_path)
  }
  invisible(plan_path)
}

publish_strict_seed_audit <- function(config, repeatability) {
  format_repeat_delta <- function(record) {
    paste0(
      "max_abs_diff=",
      ifelse(
        is.na(record$max_abs_diff),
        "NA",
        format(record$max_abs_diff, digits = 17L, scientific = TRUE)
      ),
      "; agreement=",
      ifelse(
        is.na(record$agreement),
        "NA",
        format(record$agreement, digits = 17L, scientific = TRUE)
      ),
      "; ",
      record$semantic_metric,
      "=",
      format(record$semantic_value, digits = 17L, scientific = TRUE)
    )
  }
  stages <- c(
    "atlas", "epithelial", "t_cell", "macrophage", "fibroblast",
    "b_cell", "hepatocyte_initial", "hepatocyte_clean"
  )
  expected_repeat_keys <- c(
    as.vector(outer(
      stages,
      c("FindNeighbors", "FindClusters", "RunUMAP"),
      paste,
      sep = "/"
    )),
    "covarnet/nsNMF_K5"
  )
  repeatability <- data.table::as.data.table(repeatability)
  repeatability[
    ,
    step := data.table::fcase(
      grepl("FindNeighbors", algorithm),
      paste(stage, "FindNeighbors", sep = "/"),
      grepl("FindClusters", algorithm),
      paste(stage, "FindClusters", sep = "/"),
      grepl("RunUMAP", algorithm),
      paste(stage, "RunUMAP", sep = "/"),
      grepl("NMF/nsNMF", algorithm),
      "covarnet/nsNMF_K5",
      default = NA_character_
    )
  ]
  if (
    anyNA(repeatability$step) ||
      anyDuplicated(repeatability$step) ||
      !setequal(repeatability$step, expected_repeat_keys)
  ) {
    stop("Strict seed audit withheld: repeatability suite is incomplete.")
  }
  if (
    anyNA(repeatability$exact_match) ||
      any(!repeatability$status %in% c("fixed", "blocked")) ||
      any(repeatability$status == "fixed" & !repeatability$exact_match) ||
      any(repeatability$status == "blocked" & repeatability$exact_match)
  ) {
    stop("Strict seed audit withheld: invalid executed repeat statuses.")
  }

  plan_path <- file.path(
    config$audit_root,
    "scrna_seed_audit_publication_plan.tsv"
  )
  if (!file.exists(plan_path)) {
    stop("Strict seed audit withheld: publication plan is missing.")
  }
  strict <- data.table::fread(plan_path)
  strict[, status := planned_status]
  strict[
    status == "fixed_or_blocked_after_execution",
    status := NA_character_
  ]
  for (index in seq_len(nrow(repeatability))) {
    record <- repeatability[index]
    delta <- format_repeat_delta(record)
    execution_note <- if (record$status == "blocked") {
      paste0(
        "residual non-determinism after executed duplicate: ",
        "numeric drift / output drift; ",
        delta
      )
    } else {
      paste0("Executed duplicate exact; ", delta)
    }
    strict[
      step == record$step,
      `:=`(
        status = record$status,
        notes = paste(
          notes,
          execution_note,
          sep = " | "
        )
      )
    ]
  }
  blocked_summary <- if (any(repeatability$status == "blocked")) {
    paste(
      vapply(
        which(repeatability$status == "blocked"),
        function(index) {
          record <- repeatability[index]
          paste0(record$step, " [", format_repeat_delta(record), "]")
        },
        character(1)
      ),
      collapse = "; "
    )
  } else {
    "all executed duplicate comparisons were exact"
  }
  strict[
    step == "seed_repeatability",
    `:=`(
      status = ifelse(
        any(repeatability$status == "blocked"),
        "blocked",
        "fixed"
      ),
      notes = paste(
        notes,
        if (any(repeatability$status == "blocked")) {
          paste0(
            "residual non-determinism after executed repeats: ",
            "numeric drift / output drift; ",
            blocked_summary
          )
        } else {
          blocked_summary
        },
        sep = " | "
      )
    )
  ]
  if (
    anyNA(strict$status) ||
      any(!strict$status %in% c("fixed", "deterministic", "blocked"))
  ) {
    stop("Strict seed audit withheld: unresolved publication statuses.")
  }
  strict[
    ,
    c("planned_status", "publication_gate") := NULL
  ]
  data.table::setcolorder(
    strict,
    c(
      "pipeline", "step", "script_path", "language", "algorithm", "seed",
      "seed_source", "rng_kind", "status", "notes"
    )
  )

  dependency_audit <- seed_audit_dependency_table(config)
  data.table::fwrite(strict, config$seed_audit_path, sep = "\t")
  seed_audit_mtime <- file.info(config$seed_audit_path)$mtime
  dependency_mtime <- as.POSIXct(
    dependency_audit$source_mtime_epoch,
    origin = "1970-01-01",
    tz = "UTC"
  )
  if (seed_audit_mtime <= max(dependency_mtime)) {
    unlink(config$seed_audit_path)
    stop("Strict seed audit withheld: it is not newer than dependencies.")
  }
  dependency_audit[
    ,
    `:=`(
      seed_audit_generated_utc = format(
        seed_audit_mtime,
        "%Y-%m-%dT%H:%M:%OS6Z",
        tz = "UTC"
      ),
      seed_audit_newer_than_source = seed_audit_mtime > dependency_mtime
    )
  ]
  data.table::fwrite(
    dependency_audit,
    file.path(config$audit_root, "scrna_seed_audit_dependencies.tsv"),
    sep = "\t"
  )
  invisible(config$seed_audit_path)
}

assert_no_cc <- function(values, label, prohibited = c("CC1", "CC2", "CC3")) {
  observed <- unique(as.character(values))
  hits <- intersect(observed, prohibited)
  if (length(hits) > 0L) {
    stop(
      label,
      " contains prohibited CC samples: ",
      paste(hits, collapse = ", ")
    )
  }
  invisible(TRUE)
}

assert_production_input <- function(config, path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  validation <- normalizePath(config$validation_root, mustWork = TRUE)
  if (startsWith(normalized, paste0(validation, .Platform$file.sep))) {
    stop("Production input points into validation-only root: ", normalized)
  }
  invisible(normalized)
}

validate_frozen_curation_audit <- function(config) {
  if (!file.exists(config$curation_audit_path)) {
    stop(
      "Frozen-curation audit is missing: ",
      config$curation_audit_path,
      ". Run 02_audit_frozen_curation.R first."
    )
  }
  audit <- data.table::fread(config$curation_audit_path)
  required_columns <- c("object", "path", "sha256", "status")
  if (!all(required_columns %in% colnames(audit))) {
    stop("Frozen-curation audit schema is invalid.")
  }
  if (nrow(audit) != length(config$curated_atlas_paths)) {
    stop("Frozen-curation audit does not contain all declared objects.")
  }
  if (any(audit$status != "PASS")) {
    stop("Frozen-curation audit contains non-PASS rows.")
  }
  for (path in unlist(config$curated_atlas_paths, use.names = FALSE)) {
    relative_path <- sub(
      paste0("^", config$project_root, .Platform$file.sep),
      "",
      normalizePath(path, mustWork = TRUE)
    )
    matched <- audit[audit$path == relative_path, ]
    if (nrow(matched) != 1L || sha256_file(path) != matched$sha256[[1L]]) {
      stop("Frozen-curation checksum audit is stale for: ", path)
    }
  }
  invisible(TRUE)
}

read_annotation <- function(config, filename) {
  path <- file.path(config$annotation_root, filename)
  assert_production_input(config, path)
  annotation <- data.table::fread(
    path,
    header = FALSE,
    col.names = c("cluster", "label", "evidence")
  )
  annotation[, cluster := as.character(cluster)]
  annotation
}

ensure_output_dirs <- function(config) {
  dir.create(config$result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$audit_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(
    config$downstream_curation_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
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
