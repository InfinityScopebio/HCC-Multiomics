#!/usr/bin/env Rscript

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
arguments <- parse_named_arguments(commandArgs(trailingOnly = TRUE))

required_arguments <- c("profile", "k", "seed")
missing_arguments <- setdiff(required_arguments, names(arguments))
if (length(missing_arguments) > 0L) {
  stop("Missing arguments: ", paste(missing_arguments, collapse = ", "))
}

profile <- arguments$profile
k <- as.integer(arguments$k)
seed <- as.integer(arguments$seed)
if (!is.finite(k) || !is.finite(seed)) {
  stop("K and seed must be finite integers.")
}
worker_seed_environment <- Sys.getenv("MOFA_WORKER_SEED", unset = "")
if (
  nzchar(worker_seed_environment) &&
    !identical(as.integer(worker_seed_environment), seed)
) {
  stop(
    "Worker seed mismatch: --seed=",
    seed,
    " but MOFA_WORKER_SEED=",
    worker_seed_environment
  )
}

profile_options <- switch(
  profile,
  fast1000 = list(
    convergence_mode = "fast",
    maxiter = 1000L,
    model_path = fast_model_path(cfg, k, seed),
    metadata_path = fast_model_metadata_path(cfg, k, seed),
    profile_name = "k_scan_fast1000"
  ),
  slow2000 = {
    if (k != 15L || seed != 1L) {
      stop("The final slow2000 profile is fixed to K=15, seed=1.")
    }
    list(
      convergence_mode = "slow",
      maxiter = 2000L,
      model_path = cfg$final_model_path,
      metadata_path = cfg$final_model_metadata_path,
      profile_name = "final_slow2000"
    )
  },
  fast1000_replay = {
    if (k != 15L || seed != 1L) {
      stop("The fast1000 replay profile is fixed to K=15, seed=1.")
    }
    list(
      convergence_mode = "fast",
      maxiter = 1000L,
      model_path = replay_model_path(cfg, "fast1000_replay", k, seed),
      metadata_path = replay_model_metadata_path(
        cfg,
        "fast1000_replay",
        k,
        seed
      ),
      profile_name = "k_scan_fast1000_replay"
    )
  },
  slow2000_replay = {
    if (k != 15L || seed != 1L) {
      stop("The slow2000 replay profile is fixed to K=15, seed=1.")
    }
    list(
      convergence_mode = "slow",
      maxiter = 2000L,
      model_path = replay_model_path(cfg, "slow2000_replay", k, seed),
      metadata_path = replay_model_metadata_path(
        cfg,
        "slow2000_replay",
        k,
        seed
      ),
      profile_name = "final_slow2000_replay"
    )
  },
  stop("Unknown profile: ", profile)
)

require_local_files(c(cfg$raw_input_path, cfg$python_path))
input_sha256 <- sha256_file(cfg$raw_input_path)
force_retrain <- identical(Sys.getenv("MOFA_FORCE_RETRAIN", "0"), "1")
existing_valid <- validate_production_model(
  model_path = profile_options$model_path,
  metadata_path = profile_options$metadata_path,
  input_sha256 = input_sha256,
  profile = profile_options$profile_name,
  k = k,
  seed = seed,
  maxiter = profile_options$maxiter,
  convergence_mode = profile_options$convergence_mode
)
if (existing_valid && !force_retrain) {
  cat(
    sprintf(
      "[reuse raw-derived] %s K=%d seed=%d %s\n",
      profile,
      k,
      seed,
      profile_options$model_path
    )
  )
  quit(save = "no", status = 0L)
}
if (
  !force_retrain &&
    (file.exists(profile_options$model_path) ||
      file.exists(profile_options$metadata_path))
) {
  stop(
    "Existing model or metadata has invalid provenance. Set MOFA_FORCE_RETRAIN=1 ",
    "to replace it; frozen models are never used as fallback."
  )
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
expected_rng_kind <- c("Mersenne-Twister", "Inversion", "Rejection")
if (!identical(RNGkind(), expected_rng_kind)) {
  stop("Failed to establish the historical R RNG kind.")
}
Sys.setenv(
  PYTHONHASHSEED = as.character(seed),
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)
Sys.setenv(
  LIUF_MASTER_SEED = as.character(cfg$master_seed),
  MOFA_WORKER_SEED = as.character(seed),
  SKLEARN_RANDOM_STATE = as.character(seed)
)

suppressPackageStartupMessages({
  library(data.table)
  library(reticulate)
})
reticulate::use_python(cfg$python_path, required = TRUE)
numpy <- reticulate::import("numpy", delay_load = FALSE)
numpy$random$seed(as.integer(seed))
python_random <- reticulate::import("random", delay_load = FALSE)
python_random$seed(as.integer(seed))
sklearn <- reticulate::import("sklearn", delay_load = FALSE)
sklearn$utils$check_random_state(as.integer(seed))
suppressPackageStartupMessages(library(MOFA2))

input <- readRDS(cfg$raw_input_path)
views <- input$views
sample_metadata <- data.frame(
  sample = input$common,
  group = "group1",
  stage = input$stage,
  dx = as.character(input$grp)
)

started_at <- Sys.time()
model <- create_mofa(views)
samples_metadata(model) <- sample_metadata
data_options <- get_default_data_options(model)
data_options$scale_views <- TRUE
model_options <- get_default_model_options(model)
model_options$num_factors <- k
model_options$spikeslab_weights <- TRUE
model_options$ard_weights <- TRUE
training_options <- get_default_training_options(model)
training_options$seed <- seed
training_options$convergence_mode <- profile_options$convergence_mode
training_options$maxiter <- profile_options$maxiter
model <- prepare_mofa(
  model,
  data_options = data_options,
  model_options = model_options,
  training_options = training_options
)

temporary_model <- paste0(
  profile_options$model_path,
  ".tmp.",
  Sys.getpid(),
  ".hdf5"
)
on.exit(unlink(temporary_model), add = TRUE)
model <- run_mofa(
  model,
  outfile = temporary_model,
  use_basilisk = FALSE,
  save_data = FALSE
)
if (!file.rename(temporary_model, profile_options$model_path)) {
  stop("Failed to atomically install trained model: ", profile_options$model_path)
}
finished_at <- Sys.time()

variance <- get_variance_explained(model)
elbo <- as.numeric(model@training_stats$elbo)
finite_elbo <- elbo[is.finite(elbo)]
metadata <- data.table(
  status = "complete",
  profile = profile_options$profile_name,
  K = k,
  seed = seed,
  convergence_mode = profile_options$convergence_mode,
  maxiter = profile_options$maxiter,
  scale_views = TRUE,
  spikeslab_weights = TRUE,
  ard_weights = TRUE,
  use_basilisk = FALSE,
  save_data = FALSE,
  r_rng_kind = paste(RNGkind(), collapse = ";"),
  python_hash_seed = seed,
  python_random_seed = seed,
  numpy_seed = seed,
  sklearn_random_state = seed,
  master_seed = cfg$master_seed,
  worker_seed = seed,
  worker_seed_environment = as.integer(
    Sys.getenv("MOFA_WORKER_SEED", as.character(seed))
  ),
  math_threads = 1L,
  input_path = normalizePath(cfg$raw_input_path, mustWork = TRUE),
  input_sha256 = input_sha256,
  model_path = normalizePath(profile_options$model_path, mustWork = TRUE),
  file_sha256 = sha256_file(profile_options$model_path),
  semantic_sha256 = semantic_model_hash(model),
  n_samples = nrow(get_factors(model)[[1]]),
  n_factors = ncol(get_factors(model)[[1]]),
  mean_total_variance_percent = mean(variance$r2_total[[1]]),
  sum_factor_attributed_variance_percent = (
    sum(variance$r2_per_factor[[1]]) /
      ncol(variance$r2_per_factor[[1]])
  ),
  final_finite_elbo = if (length(finite_elbo) > 0L) tail(finite_elbo, 1L) else NA_real_,
  started_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
  finished_at = format(finished_at, "%Y-%m-%dT%H:%M:%S%z"),
  elapsed_seconds = as.numeric(
    difftime(finished_at, started_at, units = "secs")
  ),
  r_version = as.character(getRversion()),
  mofa2_version = as.character(packageVersion("MOFA2")),
  python_version = as.character(reticulate::py_config()$version),
  mofapy2_version = as.character(
    reticulate::import("mofapy2")$`__version__`
  ),
  numpy_version = as.character(numpy$`__version__`)
)
fwrite(metadata, profile_options$metadata_path, sep = "\t")

cat(
  sprintf(
    "[trained] %s K=%d seed=%d elapsed=%.2fs file=%s semantic=%s\n",
    profile,
    k,
    seed,
    metadata$elapsed_seconds,
    metadata$file_sha256,
    metadata$semantic_sha256
  )
)
