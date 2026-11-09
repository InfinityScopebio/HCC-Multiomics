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

require_production_model(cfg)
require_local_files(c(
  cfg$validation$fast_k_scan_summary,
  cfg$validation$final_model
))
input <- readRDS(cfg$raw_input_path)
group_levels <- c("Healthy", "HepFibrosis", "Cirrhosis", "HCC")
n_views <- length(input$views)
if (n_views != 4L) {
  stop("Expected four serum views, found ", n_views)
}

greedy_match <- function(reference_scores, target_scores) {
  correlations <- abs(suppressWarnings(cor(reference_scores, target_scores)))
  correlations[is.na(correlations)] <- 0
  best <- numeric(ncol(reference_scores))
  for (iteration in seq_len(ncol(reference_scores))) {
    maximum <- which(
      correlations == max(correlations),
      arr.ind = TRUE
    )[1L, ]
    best[maximum[1L]] <- correlations[maximum[1L], maximum[2L]]
    correlations[maximum[1L], ] <- -1
    correlations[, maximum[2L]] <- -1
  }
  best
}

strict_monotonic <- function(values) {
  means <- tapply(values, input$grp, mean)[group_levels]
  if (any(is.na(means))) {
    return(FALSE)
  }
  all(diff(means) > 0) || all(diff(means) < 0)
}

load_scan_run <- function(k, seed) {
  model_path <- fast_model_path(cfg, k, seed)
  metadata_path <- fast_model_metadata_path(cfg, k, seed)
  valid <- validate_production_model(
    model_path = model_path,
    metadata_path = metadata_path,
    input_sha256 = sha256_file(cfg$raw_input_path),
    profile = "k_scan_fast1000",
    k = k,
    seed = seed,
    maxiter = 1000L,
    convergence_mode = "fast"
  )
  if (!valid) {
    stop(
      "Missing or invalid raw-derived K-scan model: K=",
      k,
      " seed=",
      seed
    )
  }
  model <- load_model(model_path, verbose = FALSE)
  variance <- get_variance_explained(model)
  scores <- get_factors(model)[[1]]
  if (!all(input$common %in% rownames(scores))) {
    stop("K=", k, " seed=", seed, " lacks input samples.")
  }
  list(
    Z = scores[input$common, , drop = FALSE],
    VE = variance$r2_per_factor[[1]],
    VE_total = variance$r2_total[[1]],
    metadata = fread(metadata_path)
  )
}

runs <- setNames(vector("list", length(cfg$k_grid)), as.character(cfg$k_grid))
for (k in cfg$k_grid) {
  runs[[as.character(k)]] <- setNames(
    lapply(cfg$k_scan_seeds, function(seed) load_scan_run(k, seed)),
    as.character(cfg$k_scan_seeds)
  )
}

detail <- data.table()
summary <- data.table()
for (k in cfg$k_grid) {
  run_set <- runs[[as.character(k)]]
  reference <- run_set[["1"]]
  matches <- sapply(
    run_set[names(run_set) != "1"],
    function(target) greedy_match(reference$Z, target$Z)
  )
  if (is.null(dim(matches))) {
    matches <- matrix(matches, nrow = ncol(reference$Z))
  }
  reproducibility <- apply(matches, 1, median)
  variance_sum <- rowSums(reference$VE)
  variance_mean <- rowMeans(reference$VE)
  mean_view_variance <- vapply(
    run_set,
    function(run) mean(run$VE_total),
    numeric(1)
  )
  factor_attributed_variance <- vapply(
    run_set,
    function(run) sum(run$VE) / ncol(run$VE),
    numeric(1)
  )
  factor_names <- colnames(reference$Z)
  detail <- rbind(
    detail,
    data.table(
      K = k,
      factor = factor_names,
      VE_total = round(variance_sum, 3),
      VE_sum_views_pp = round(variance_sum, 3),
      VE_mean_views = round(variance_mean, 3),
      reproducibility = round(reproducibility, 3),
      active = variance_sum >= 1,
      robust = variance_sum >= 1 & reproducibility >= 0.8,
      strict_monotonic = apply(reference$Z, 2, strict_monotonic),
      rho_stage = round(
        apply(
          reference$Z,
          2,
          cor,
          y = input$stage,
          method = "spearman"
        ),
        3
      )
    )
  )
  summary <- rbind(
    summary,
    data.table(
      K = k,
      n_factors = ncol(reference$Z),
      n_active = sum(variance_sum >= 1),
      n_robust = sum(variance_sum >= 1 & reproducibility >= 0.8),
      n_strict_robust = sum(
        variance_sum >= 1 &
          reproducibility >= 0.8 &
          apply(reference$Z, 2, strict_monotonic)
      ),
      total_VE = round(sum(variance_sum), 2),
      mean_repro_active = round(
        mean(reproducibility[variance_sum >= 1]),
        3
      ),
      mean_view_VE_seed1 = mean(reference$VE_total),
      median_mean_view_VE = median(mean_view_variance),
      min_mean_view_VE = min(mean_view_variance),
      max_mean_view_VE = max(mean_view_variance),
      mean_factor_attributed_VE_seed1 = sum(reference$VE) / ncol(reference$VE),
      median_factor_attributed_VE = median(factor_attributed_variance),
      min_factor_attributed_VE = min(factor_attributed_variance),
      max_factor_attributed_VE = max(factor_attributed_variance),
      n_successful_seeds = length(run_set)
    )
  )
}

saveRDS(
  list(
    provenance = list(
      source = "raw_M0_retrained",
      input_path = normalizePath(cfg$raw_input_path, mustWork = TRUE),
      input_sha256 = sha256_file(cfg$raw_input_path),
      K = cfg$k_grid,
      seeds = cfg$k_scan_seeds,
      convergence_mode = "fast",
      maxiter = 1000L,
      scale_views = TRUE,
      spikeslab_weights = TRUE,
      ard_weights = TRUE
    ),
    runs = lapply(runs, function(run_set) {
      lapply(run_set, function(run) run[c("Z", "VE", "VE_total")])
    })
  ),
  file.path(cfg$result_dir, "k_scan_fast1000_runs_raw_derived.rds")
)
fwrite(
  detail,
  file.path(cfg$result_dir, "k_scan_fast1000_factor_detail_raw_derived.tsv"),
  sep = "\t"
)
fwrite(
  summary,
  file.path(cfg$result_dir, "k_scan_fast1000_summary_raw_derived.tsv"),
  sep = "\t"
)
fwrite(
  summary,
  file.path(cfg$result_dir, "k_scan_accepted_fast1000.tsv"),
  sep = "\t"
)

manifest_paths <- c(
  unlist(lapply(cfg$k_grid, function(k) {
    vapply(
      cfg$k_scan_seeds,
      function(seed) fast_model_metadata_path(cfg, k, seed),
      character(1)
    )
  })),
  cfg$final_model_metadata_path,
  replay_model_metadata_path(cfg, "fast1000_replay", 15L, 1L),
  replay_model_metadata_path(cfg, "slow2000_replay", 15L, 1L)
)
manifest <- rbindlist(lapply(manifest_paths, fread), fill = TRUE)
setorder(manifest, profile, K, seed)
fwrite(
  manifest,
  file.path(cfg$result_dir, "mofa_training_manifest.tsv"),
  sep = "\t"
)

validation_summary <- fread(cfg$validation$fast_k_scan_summary)[order(K)]
validation_columns <- c(
  "K",
  "n_factors",
  "n_active",
  "n_robust",
  "n_strict_robust",
  "total_VE",
  "mean_repro_active"
)
if (!all(validation_columns %in% names(validation_summary))) {
  stop("Validation K-scan target lacks required columns.")
}
comparison <- merge(
  validation_summary[, ..validation_columns],
  summary[, ..validation_columns],
  by = "K",
  suffixes = c("_target", "_raw_derived"),
  all = TRUE,
  sort = TRUE
)
for (column in setdiff(validation_columns, "K")) {
  comparison[
    ,
    paste0(column, "_abs_diff") := abs(
      get(paste0(column, "_target")) -
        get(paste0(column, "_raw_derived"))
    )
  ]
}
fwrite(
  comparison,
  file.path(cfg$result_dir, "k_scan_fast1000_vs_validation.tsv"),
  sep = "\t"
)

final_model <- load_model(cfg$final_model_path, verbose = FALSE)
final_variance <- get_variance_explained(final_model)
factor_view <- as.data.table(as.table(final_variance$r2_per_factor[[1]]))
setnames(factor_view, c("factor", "view", "variance_percent"))
view_total <- data.table(
  view = names(final_variance$r2_total[[1]]),
  variance_percent = as.numeric(final_variance$r2_total[[1]])
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

validation_model <- load_model(cfg$validation$final_model, verbose = FALSE)
native_scores_diff <- max(abs(
  get_factors(final_model)[[1]] -
    get_factors(validation_model)[[1]]
))
native_weights_diff <- max(vapply(names(get_weights(final_model)), function(view) {
  max(abs(get_weights(final_model)[[view]] - get_weights(validation_model)[[view]]))
}, numeric(1)))
native_variance_diff <- max(abs(
  final_variance$r2_per_factor[[1]] -
    get_variance_explained(validation_model)$r2_per_factor[[1]]
))
matched_score_correlations <- greedy_match(
  get_factors(validation_model)[[1]],
  get_factors(final_model)[[1]]
)

fast_primary_model <- load_model(fast_model_path(cfg, 15L, 1L), verbose = FALSE)
fast_replay_model <- load_model(
  replay_model_path(cfg, "fast1000_replay", 15L, 1L),
  verbose = FALSE
)
final_replay_model <- load_model(
  replay_model_path(cfg, "slow2000_replay", 15L, 1L),
  verbose = FALSE
)

compare_models <- function(first, second) {
  first_variance <- get_variance_explained(first)
  second_variance <- get_variance_explained(second)
  list(
    score_max_abs_diff = max(abs(
      get_factors(first)[[1]] -
        get_factors(second)[[1]]
    )),
    weight_max_abs_diff = max(vapply(names(get_weights(first)), function(view) {
      max(abs(get_weights(first)[[view]] - get_weights(second)[[view]]))
    }, numeric(1))),
    variance_max_abs_diff = max(abs(
      first_variance$r2_per_factor[[1]] -
        second_variance$r2_per_factor[[1]]
    )),
    semantic_hash_equal = identical(
      semantic_model_hash(first),
      semantic_model_hash(second)
    )
  )
}
fast_repeatability <- compare_models(fast_primary_model, fast_replay_model)
final_repeatability <- compare_models(final_model, final_replay_model)

numeric_columns <- setdiff(validation_columns, c("K"))
k_scan_numeric_max_diff <- max(
  unlist(
    comparison[
      ,
      paste0(numeric_columns, "_abs_diff"),
      with = FALSE
    ]
  ),
  na.rm = TRUE
)
k_scan_count_exact <- all(
  comparison$n_factors_target == comparison$n_factors_raw_derived &
    comparison$n_active_target == comparison$n_active_raw_derived &
    comparison$n_robust_target == comparison$n_robust_raw_derived &
    comparison$n_strict_robust_target ==
      comparison$n_strict_robust_raw_derived
)

audit <- data.table(
  check_id = c(
    "raw_derived_training_grid",
    "fast1000_k_scan_counts",
    "fast1000_k_scan_numeric",
    "final_k15_native_scores",
    "final_k15_native_weights",
    "final_k15_native_variance",
    "final_k15_factor_recovery",
    "final_k15_file_hash",
    "final_k15_semantic_hash",
    "fast1000_repeatability",
    "slow2000_repeatability"
  ),
  status = c(
    ifelse(nrow(manifest[profile == "k_scan_fast1000"]) == 114L, "exact", "incomplete"),
    ifelse(k_scan_count_exact, "exact", "mismatch"),
    ifelse(k_scan_numeric_max_diff <= 1e-10, "numeric_exact", "numeric_drift"),
    ifelse(native_scores_diff <= 1e-10, "numeric_exact", "numeric_drift"),
    ifelse(native_weights_diff <= 1e-10, "numeric_exact", "numeric_drift"),
    ifelse(native_variance_diff <= 1e-10, "numeric_exact", "numeric_drift"),
    ifelse(min(matched_score_correlations) >= 1 - 1e-10, "exact", "numeric_drift"),
    ifelse(
      sha256_file(cfg$final_model_path) == sha256_file(cfg$validation$final_model),
      "file_exact",
      "file_hash_differs"
    ),
    ifelse(
      semantic_model_hash(final_model) == semantic_model_hash(validation_model),
      "semantic_exact",
      "semantic_hash_differs"
    ),
    ifelse(
      max(
        fast_repeatability$score_max_abs_diff,
        fast_repeatability$weight_max_abs_diff,
        fast_repeatability$variance_max_abs_diff
      ) <= 1e-12 && fast_repeatability$semantic_hash_equal,
      "numeric_and_semantic_exact",
      "numeric_drift"
    ),
    ifelse(
      max(
        final_repeatability$score_max_abs_diff,
        final_repeatability$weight_max_abs_diff,
        final_repeatability$variance_max_abs_diff
      ) <= 1e-12 && final_repeatability$semantic_hash_equal,
      "numeric_and_semantic_exact",
      "numeric_drift"
    )
  ),
  observed = c(
    sprintf(
      "%d fast models plus %d final model",
      nrow(manifest[profile == "k_scan_fast1000"]),
      nrow(manifest[profile == "final_slow2000"])
    ),
    sprintf("count_rows_exact=%s", k_scan_count_exact),
    sprintf("max_abs_diff=%.12g", k_scan_numeric_max_diff),
    sprintf("max_abs_diff=%.12g", native_scores_diff),
    sprintf("max_abs_diff=%.12g", native_weights_diff),
    sprintf("max_abs_diff=%.12g", native_variance_diff),
    sprintf("minimum matched |cor|=%.12g", min(matched_score_correlations)),
    paste(
      sha256_file(cfg$final_model_path),
      sha256_file(cfg$validation$final_model),
      sep = " vs "
    ),
    paste(
      semantic_model_hash(final_model),
      semantic_model_hash(validation_model),
      sep = " vs "
    ),
    sprintf(
      "scores=%.12g;weights=%.12g;variance=%.12g;semantic=%s",
      fast_repeatability$score_max_abs_diff,
      fast_repeatability$weight_max_abs_diff,
      fast_repeatability$variance_max_abs_diff,
      fast_repeatability$semantic_hash_equal
    ),
    sprintf(
      "scores=%.12g;weights=%.12g;variance=%.12g;semantic=%s",
      final_repeatability$score_max_abs_diff,
      final_repeatability$weight_max_abs_diff,
      final_repeatability$variance_max_abs_diff,
      final_repeatability$semantic_hash_equal
    )
  ),
  note = c(
    "All K=2:20 x seeds=1:6 models were trained from raw-derived input.",
    "Frozen summary is validation-only and never supplies production values.",
    "Numeric drift is reported and never replaced by target values.",
    "Native factor order and sign are compared without target-based alignment.",
    "Native factor order and sign are compared without target-based alignment.",
    "Variance is recomputed from the newly trained final HDF5.",
    "Greedy matching is audit-only; downstream uses native raw-derived factors.",
    "HDF5 byte hashes may differ despite equivalent numeric model content.",
    "Semantic hash covers factors, weights, variance and model options.",
    "Independent replay uses the same raw input/config/seed and is validation-only.",
    "Independent replay uses the same raw input/config/seed and is validation-only."
  )
)
fwrite(
  audit,
  file.path(cfg$audit_dir, "mofa_k_scan_variance_audit.tsv"),
  sep = "\t"
)

model_seed_source <- function(profile) {
  if (grepl("replay", profile)) {
    return("repro repeatability requirement")
  }
  if (profile == "k_scan_fast1000") {
    return(
      "historical run_l5_step49e_kselect.R: K=2:20, seeds=1:6"
    )
  }
  "historical run_l5_step49g_K15.R: final K15 seed=1"
}
model_seed_reason <- function(profile) {
  if (grepl("replay", profile)) {
    return("independent deterministic replay of the corresponding production fit")
  }
  "preserve the explicitly documented historical model seed"
}

model_seed_rows <- rbindlist(lapply(seq_len(nrow(manifest)), function(index) {
  record <- manifest[index]
  worker_script <- if (record$profile == "k_scan_fast1000") {
    "code/pipelines/60_serum_mofa/16_train_fast_worker.sh"
  } else {
    "code/pipelines/60_serum_mofa/run_all.sh"
  }
  shared <- list(
    step = sprintf("%s_K%d_seed%d", record$profile, record$K, record$seed),
    seed_source = model_seed_source(record$profile),
    notes = sprintf(
      paste0(
        "K=%d;convergence_mode=%s;maxiter=%d;scale_views=TRUE;",
        "spikeslab_weights=TRUE;ard_weights=TRUE;math_threads=1;",
        "reason=%s"
      ),
      record$K,
      record$convergence_mode,
      record$maxiter,
      model_seed_reason(record$profile)
    )
  )
  rbind(
    data.table(
      pipeline = "60_serum_mofa",
      step = shared$step,
      script_path = "code/pipelines/60_serum_mofa/15_train_mofa_model.R",
      language = "R/Python",
      algorithm = "MOFA2/mofapy2 variational Bayes",
      seed = record$seed,
      seed_source = shared$seed_source,
      rng_kind = "MOFA training_options$seed",
      status = "fixed",
      notes = shared$notes
    ),
    data.table(
      pipeline = "60_serum_mofa",
      step = paste0(shared$step, "_R"),
      script_path = "code/pipelines/60_serum_mofa/15_train_mofa_model.R",
      language = "R",
      algorithm = "per-process R RNG",
      seed = record$seed,
      seed_source = shared$seed_source,
      rng_kind = "Mersenne-Twister/Inversion/Rejection",
      status = "fixed",
      notes = paste(shared$notes, "RNG fixed before model construction")
    ),
    data.table(
      pipeline = "60_serum_mofa",
      step = paste0(shared$step, "_Python"),
      script_path = "code/pipelines/60_serum_mofa/15_train_mofa_model.R",
      language = "Python",
      algorithm = "Python random and hash seed",
      seed = record$python_random_seed,
      seed_source = shared$seed_source,
      rng_kind = "random.seed/PYTHONHASHSEED",
      status = "fixed",
      notes = paste0(
        shared$notes,
        ";PYTHONHASHSEED=",
        record$python_hash_seed
      )
    ),
    data.table(
      pipeline = "60_serum_mofa",
      step = paste0(shared$step, "_NumPy"),
      script_path = "code/pipelines/60_serum_mofa/15_train_mofa_model.R",
      language = "Python",
      algorithm = "NumPy RNG used by mofapy2",
      seed = record$numpy_seed,
      seed_source = shared$seed_source,
      rng_kind = "NumPy RandomState/MT19937",
      status = "fixed",
      notes = shared$notes
    ),
    data.table(
      pipeline = "60_serum_mofa",
      step = paste0(shared$step, "_sklearn_namespace"),
      script_path = "code/pipelines/60_serum_mofa/15_train_mofa_model.R",
      language = "Python",
      algorithm = "sklearn.utils.check_random_state",
      seed = record$sklearn_random_state,
      seed_source = if (record$sklearn_random_state == record$seed) {
        shared$seed_source
      } else {
        "project master seed in completed fit; sklearn estimator not used by MOFA"
      },
      rng_kind = "scikit-learn RandomState/MT19937",
      status = "fixed",
      notes = paste0(
        shared$notes,
        ";namespace_only_no_sklearn_estimator=",
        record$sklearn_random_state
      )
    ),
    data.table(
      pipeline = "60_serum_mofa",
      step = paste0(shared$step, "_worker"),
      script_path = worker_script,
      language = "Bash/R",
      algorithm = "independent parallel training worker",
      seed = record$worker_seed,
      seed_source = shared$seed_source,
      rng_kind = "independent process-local RNG",
      status = "fixed",
      notes = paste0(
        shared$notes,
        ";worker_seed=",
        record$worker_seed,
        ";scheduler_order_does_not_change_process_rng"
      )
    )
  )
}))

seed_audit <- rbind(
  data.table(
    pipeline = "60_serum_mofa",
    step = "k_scan_summary_and_variance",
    script_path = "code/pipelines/60_serum_mofa/20_k_scan_and_variance.R",
    language = "R",
    algorithm = "MOFA model loading, factor matching, and variance aggregation",
    seed = NA_integer_,
    seed_source = "not applicable; this stage consumes no RNG",
    rng_kind = "none",
    status = "deterministic",
    notes = paste(
      "Reads only raw-derived production HDF5 files;",
      "frozen HDF5 and summaries are loaded afterward as validation targets"
    )
  ),
  data.table(
    pipeline = "60_serum_mofa",
    step = "raw_M0_preprocessing",
    script_path = "code/pipelines/60_serum_mofa/10_rebuild_incdec_input.R",
    language = "R",
    algorithm = "impute::impute.knn and strict-monotonic filter",
    seed = cfg$master_seed,
    seed_source = sprintf(
      "LIUF_MASTER_SEED=%d (default is 20260807)",
      cfg$master_seed
    ),
    rng_kind = "Mersenne-Twister/Inversion/Rejection",
    status = "fixed",
    notes = paste(
      "RNGkind=Mersenne-Twister/Inversion/Rejection;",
      "input=results/00_preprocess/multiomics/dat"
    )
  ),
  data.table(
    pipeline = "60_serum_mofa",
    step = "global_master_seed",
    script_path = "code/pipelines/60_serum_mofa/run_all.sh",
    language = "Bash",
    algorithm = "pipeline master seed namespace",
    seed = cfg$master_seed,
    seed_source = sprintf(
      "user requirement; LIUF_MASTER_SEED=%d (default is 20260807)",
      cfg$master_seed
    ),
    rng_kind = "environment seed",
    status = "fixed",
    notes = sprintf("LIUF_MASTER_SEED=%d", cfg$master_seed)
  ),
  data.table(
    pipeline = "60_serum_mofa",
    step = "k_scan_parallel_scheduler",
    script_path = "code/pipelines/60_serum_mofa/run_all.sh",
    language = "Bash",
    algorithm = "bounded Bash independent-process worker pool",
    seed = NA_integer_,
    seed_source = "not applicable; scheduler consumes no RNG",
    rng_kind = "none",
    status = "deterministic",
    notes = paste(
      "The bounded worker pool only schedules independent processes;",
      "each worker receives K and seed as arguments and exports",
      "MOFA_WORKER_SEED/PYTHONHASHSEED/SKLEARN_RANDOM_STATE before R starts"
    )
  ),
  model_seed_rows
)
fwrite(
  seed_audit,
  file.path(cfg$audit_dir, "mofa_seed_audit.tsv"),
  sep = "\t"
)


cat(
  sprintf(
    paste0(
      "Summarized 114 raw-derived fast/1000 models and one raw-derived ",
      "slow/2000 K15 model; validation max differences: K-scan %.12g, ",
      "scores %.12g, weights %.12g.\n"
    ),
    k_scan_numeric_max_diff,
    native_scores_diff,
    native_weights_diff
  )
)
