#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(e1071)
  library(glmnet)
  library(MASS)
  library(nnet)
  library(MOFA2)
  library(pROC)
  library(randomForest)
  library(tidyr)
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
classifier_started_at <- Sys.time()
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(cfg$master_seed)

rebuilt_input_path <- cfg$raw_input_path
require_analysis_mofa_boundary(cfg)
require_local_files(c(
  rebuilt_input_path,
  cfg$analysis_model_path,
  cfg$validation$model_comparison,
  cfg$validation$performance,
  cfg$validation$permutation_b1000
))

input <- readRDS(rebuilt_input_path)
gcms <- input$views$GCMS
feature_names <- rownames(gcms)
all_features <- scale(t(gcms), center = TRUE, scale = TRUE)
rownames(all_features) <- input$common
colnames(all_features) <- feature_names
truth <- factor(
  ifelse(
    input$grp %in% c("HepFibrosis", "Cirrhosis"),
    "FibCir",
    as.character(input$grp)
  ),
  levels = c("Healthy", "FibCir", "HCC")
)

production_model <- load_model(cfg$analysis_model_path, verbose = FALSE)
weights <- get_weights(production_model)$GCMS[, "Factor2"]
weights <- weights[feature_names]
weights <- weights[!is.na(weights)]
feature_rank <- names(sort(-abs(weights)))
selected_features <- feature_rank[1:7]
levels_truth <- levels(truth)

balanced_accuracy <- function(predicted) {
  mean(vapply(
    levels_truth,
    function(class_name) {
      sum(
        truth == class_name & predicted == class_name,
        na.rm = TRUE
      ) / sum(truth == class_name)
    },
    numeric(1)
  ))
}

lda_accuracy <- function(x) {
  predicted <- vapply(seq_len(nrow(x)), function(index) {
    set.seed(derive_r_seed(cfg$master_seed, index))
    fit <- lda(x = x[-index, , drop = FALSE], grouping = truth[-index])
    as.character(stats::predict(fit, x[index, , drop = FALSE])$class)
  }, character(1))
  balanced_accuracy(predicted)
}

random_forest_accuracy <- function(x) {
  set.seed(1)
  fit <- randomForest(x = x, y = truth)
  balanced_accuracy(fit$predicted)
}

svm_accuracy <- function(x) {
  predicted <- vapply(seq_len(nrow(x)), function(index) {
    set.seed(derive_r_seed(cfg$master_seed, index))
    fit <- svm(
      x = x[-index, , drop = FALSE],
      y = truth[-index],
      kernel = "radial"
    )
    as.character(stats::predict(fit, x[index, , drop = FALSE]))
  }, character(1))
  balanced_accuracy(predicted)
}

mlr_accuracy <- function(x) {
  predicted <- vapply(seq_len(nrow(x)), function(index) {
    set.seed(derive_r_seed(cfg$master_seed, index))
    fit <- tryCatch(
      multinom(
        truth[-index] ~ .,
        data = as.data.frame(x[-index, , drop = FALSE]),
        trace = FALSE,
        maxit = 300
      ),
      error = function(error) NULL
    )
    if (is.null(fit)) {
      return(NA_character_)
    }
    as.character(stats::predict(fit, as.data.frame(x[index, , drop = FALSE])))
  }, character(1))
  balanced_accuracy(predicted)
}

glmnet_accuracy <- function(x, alpha) {
  set.seed(1)
  predicted <- character(nrow(x))
  for (index in seq_len(nrow(x))) {
    fit <- tryCatch(
      cv.glmnet(
        x[-index, , drop = FALSE],
        truth[-index],
        family = "multinomial",
        alpha = alpha,
        nfolds = 5
      ),
      error = function(error) NULL
    )
    if (is.null(fit)) {
      predicted[index] <- NA_character_
      next
    }
    new_x <- matrix(x[index, ], nrow = 1)
    colnames(new_x) <- colnames(x)
    prediction <- tryCatch(
      stats::predict(
        fit,
        newx = new_x,
        s = "lambda.1se",
        type = "class"
      ),
      error = function(error) NULL
    )
    predicted[index] <- if (is.null(prediction)) {
      NA_character_
    } else {
      as.character(prediction)[1]
    }
  }
  balanced_accuracy(predicted)
}

top_n <- 3:12
comparison <- data.table(N = top_n)
comparison[, LDA := vapply(top_n, function(n) {
  lda_accuracy(all_features[, feature_rank[1:n], drop = FALSE])
}, numeric(1))]
comparison[, RF := vapply(top_n, function(n) {
  random_forest_accuracy(all_features[, feature_rank[1:n], drop = FALSE])
}, numeric(1))]
comparison[, SVM := vapply(top_n, function(n) {
  svm_accuracy(all_features[, feature_rank[1:n], drop = FALSE])
}, numeric(1))]
comparison[, MLR := vapply(top_n, function(n) {
  mlr_accuracy(all_features[, feature_rank[1:n], drop = FALSE])
}, numeric(1))]
comparison[, LASSO := vapply(top_n, function(n) {
  glmnet_accuracy(all_features[, feature_rank[1:n], drop = FALSE], 1)
}, numeric(1))]
comparison[, Ridge := vapply(top_n, function(n) {
  glmnet_accuracy(all_features[, feature_rank[1:n], drop = FALSE], 0)
}, numeric(1))]
comparison[, ENet := vapply(top_n, function(n) {
  glmnet_accuracy(all_features[, feature_rank[1:n], drop = FALSE], 0.5)
}, numeric(1))]
fwrite(
  comparison,
  file.path(cfg$result_dir, "factor2_model_comparison.tsv"),
  sep = "\t"
)

selected_matrix <- all_features[, selected_features, drop = FALSE]
safe_matrix <- selected_matrix
colnames(safe_matrix) <- paste0("V", seq_len(ncol(safe_matrix)))
set.seed(cfg$master_seed)
loocv <- lapply(seq_len(nrow(safe_matrix)), function(index) {
  set.seed(derive_r_seed(cfg$master_seed, index))
  fit <- multinom(
    truth[-index] ~ .,
    data = as.data.frame(safe_matrix[-index, , drop = FALSE]),
    trace = FALSE,
    maxit = 300
  )
  list(
    class = as.character(
      stats::predict(
        fit,
        as.data.frame(safe_matrix[index, , drop = FALSE]),
        type = "class"
      )
    ),
    probabilities = as.numeric(
      stats::predict(
        fit,
        as.data.frame(safe_matrix[index, , drop = FALSE]),
        type = "probs"
      )
    )
  )
})
predicted <- vapply(loocv, function(result) result$class, character(1))
probability_matrix <- do.call(
  rbind,
  lapply(loocv, function(result) result$probabilities)
)
colnames(probability_matrix) <- levels_truth

confusion <- table(Actual = truth, Predicted = predicted)
confusion_normalized <- prop.table(confusion, margin = 1) * 100
confusion_table <- as.data.table(as.data.frame(confusion_normalized))
fwrite(
  confusion_table,
  file.path(cfg$result_dir, "factor2_mlr_confusion_row_percent.tsv"),
  sep = "\t"
)

metric_rows <- lapply(levels_truth, function(class_name) {
  binary_truth <- as.integer(truth == class_name)
  roc_object <- roc(
    binary_truth,
    probability_matrix[, class_name],
    quiet = TRUE,
    levels = c(0, 1),
    direction = "<"
  )
  true_positive <- sum(truth == class_name & predicted == class_name)
  false_positive <- sum(truth != class_name & predicted == class_name)
  false_negative <- sum(truth == class_name & predicted != class_name)
  precision <- true_positive / (true_positive + false_positive)
  recall <- true_positive / (true_positive + false_negative)
  f1 <- 2 * precision * recall / (precision + recall)
  data.table(
    class = class_name,
    AUC = as.numeric(auc(roc_object)),
    precision = precision,
    recall = recall,
    F1 = f1
  )
})
performance <- rbindlist(metric_rows)
fwrite(
  performance,
  file.path(cfg$result_dir, "factor2_mlr_performance.tsv"),
  sep = "\t"
)

predictions <- data.table(
  sample = rownames(safe_matrix),
  truth = as.character(truth),
  predicted = predicted
)
predictions <- cbind(predictions, as.data.table(probability_matrix))
fwrite(
  predictions,
  file.path(cfg$result_dir, "factor2_mlr_loocv_predictions.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    rank = seq_along(selected_features),
    feature = selected_features,
    factor2_loading = weights[selected_features]
  ),
  file.path(cfg$result_dir, "factor2_top7_features.tsv"),
  sep = "\t"
)

permutation_count <- as.integer(Sys.getenv("MOFA_PERM_B", "1000"))
if (!identical(permutation_count, 1000L)) {
  stop(
    "The production audit requires the complete historical B=1000 permutation run."
  )
}
permutation_summary <- NULL
if (is.finite(permutation_count) && permutation_count > 0L) {
  permutation_predict <- function(x) {
    predicted_perm <- vapply(seq_len(nrow(x)), function(index) {
      fit <- tryCatch(
        multinom(
          truth[-index] ~ .,
          data = as.data.frame(x[-index, , drop = FALSE]),
          trace = FALSE,
          maxit = 400
        ),
        error = function(error) NULL,
        warning = function(warning) NULL
      )
      if (is.null(fit)) {
        return(NA_character_)
      }
      as.character(
        stats::predict(fit, as.data.frame(x[index, , drop = FALSE]))
      )
    }, character(1))
    factor(predicted_perm, levels = levels_truth)
  }

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(42)
  baseline_accuracy <- balanced_accuracy(permutation_predict(selected_matrix))
  permutation_rows <- lapply(selected_features, function(feature_name) {
      drops <- vapply(seq_len(permutation_count), function(iteration) {
        permuted <- selected_matrix
        permuted[, feature_name] <- sample(permuted[, feature_name])
        baseline_accuracy - balanced_accuracy(permutation_predict(permuted))
      }, numeric(1))
      cat(
        sprintf(
          "[B=1000 serial] %-24s mean_drop=%.12f\n",
          feature_name,
          mean(drops)
        )
      )
      data.table(
        feature = feature_name,
        iteration = seq_len(permutation_count),
        drop = drops,
        permutation_seed = 42L
      )
  })
  permutation_data <- rbindlist(permutation_rows)
  fwrite(
    permutation_data,
    file.path(cfg$result_dir, "factor2_permutation_drops_B1000.tsv"),
    sep = "\t"
  )
  permutation_summary <- permutation_data[, .(
    mean = mean(drop),
    sd = sd(drop),
    q_lo = quantile(drop, 0.025),
    q_hi = quantile(drop, 0.975),
    frac_positive = mean(drop > 0)
  ), by = feature]
  setorder(permutation_summary, -mean)
  fwrite(
    permutation_summary,
    file.path(
      cfg$result_dir,
      sprintf("factor2_permutation_importance_B%d.tsv", permutation_count)
    ),
    sep = "\t"
  )
}

historical_comparison <- fread(cfg$validation$model_comparison)
comparison_diff <- max(
  abs(
    as.matrix(historical_comparison[, -1]) -
      as.matrix(comparison[, -1])
  )
)
historical_performance <- fread(cfg$validation$performance)
performance_diff <- max(
  abs(
    as.matrix(historical_performance[, c("AUC", "precision", "recall", "F1")]) -
      round(as.matrix(performance[, c("AUC", "precision", "recall", "F1")]), 3)
  )
)
permutation_maximum_difference <- NA_real_
permutation_status <- if (permutation_count == 1000L) {
  historical_permutation <- fread(cfg$validation$permutation_b1000)
  setorder(historical_permutation, feature)
  setorder(permutation_summary, feature)
  permutation_maximum_difference <- max(
    abs(
      as.matrix(
        historical_permutation[, c("mean", "sd", "q_lo", "q_hi", "frac_positive")]
      ) -
        as.matrix(
          permutation_summary[, c("mean", "sd", "q_lo", "q_hi", "frac_positive")]
      )
    )
  )
  ifelse(
    permutation_maximum_difference <= 1e-12,
    "numeric_exact",
    "numeric_mismatch"
  )
} else if (permutation_count > 0L) {
  sprintf("deterministic_smoke_B%d", permutation_count)
} else {
  "not_run"
}

audit <- data.table(
  check_id = c(
    "top7_feature_selection",
    "top_n_model_comparison",
    "top7_mlr_loocv",
    "permutation_importance",
    "validation_design"
  ),
  status = c(
    ifelse(
      identical(
        selected_features,
        c(
          "Palmitic_Acid",
          "Linoleic_Acid",
          "Trans-Vaccenic_Acid",
          "Stearic_Acid",
          "Glycerol",
          "17-Octadecynoic_Acid",
          "D-Fructose"
        )
      ),
      "exact",
      "mismatch"
    ),
    ifelse(comparison_diff <= 1e-12, "numeric_exact", "numeric_mismatch"),
    ifelse(performance_diff <= 1e-12, "numeric_exact", "numeric_mismatch"),
    permutation_status,
    "caveat"
  ),
  observed = c(
    paste(selected_features, collapse = ";"),
    sprintf("max_abs_diff=%.3g", comparison_diff),
    sprintf(
      "balanced_accuracy=%.12f; rounded metric max_abs_diff=%.3g",
      balanced_accuracy(predicted),
      performance_diff
    ),
    sprintf(
      "B=%d; max_abs_diff=%.12g",
      permutation_count,
      permutation_maximum_difference
    ),
    paste(
      "MOFA loading selection used all 54 samples before LOOCV;",
      "RF uses OOB while six other models use LOOCV."
    )
  ),
  note = c(
    "Top seven are ranked by absolute Factor2 GCMS loading.",
    "Seeds and package behavior reproduce the historical table.",
    "The MLR itself is deterministic; no random seed is consumed by multinom.",
    "Historical exact comparison is only defined for seed 42 and B=1000.",
    "This is exploratory same-cohort, non-nested performance, not an unbiased external estimate."
  )
)
fwrite(
  audit,
  file.path(cfg$audit_dir, "mofa_factor2_classifier_audit.tsv"),
  sep = "\t"
)

classifier_seed_audit <- data.table(
  pipeline = "60_serum_mofa",
  step = c(
    "classifier_global_R",
    "classifier_seed_derivation_helper",
    "classifier_LDA_LOOCV",
    "classifier_RF_OOB",
    "classifier_SVM_LOOCV",
    "classifier_MLR_LOOCV",
    "classifier_glmnet_LOOCV_CV",
    "classifier_permutation_importance_B1000_historical_serial_stream"
  ),
  seed = c(
    cfg$master_seed,
    cfg$master_seed,
    cfg$master_seed,
    1L,
    cfg$master_seed,
    cfg$master_seed,
    1L,
    42L
  ),
  script_path = c(
    "code/pipelines/60_serum_mofa/40_factor2_classifier.R",
    "code/lib/R/deterministic_rng.R",
    rep(
      "code/pipelines/60_serum_mofa/40_factor2_classifier.R",
      6L
    )
  ),
  language = "R",
  algorithm = c(
    "R RNG Mersenne-Twister/Inversion/Rejection",
    "deterministic modulo-2^31 seed derivation",
    "MASS::lda leave-one-out refit",
    "randomForest::randomForest OOB",
    "e1071::svm radial leave-one-out refit",
    "nnet::multinom leave-one-out refit",
    "glmnet::cv.glmnet multinomial 5-fold inner CV",
    "nnet::multinom historical serial permutation stream"
  ),
  seed_source = c(
    "LIUF_MASTER_SEED default 20260807",
    "LIUF_MASTER_SEED plus audited offset modulo 2^31",
    "LIUF_MASTER_SEED plus held-out index modulo 2^31",
    "historical run_Fig1_model_comparison.R",
    "LIUF_MASTER_SEED plus held-out index modulo 2^31",
    paste(
      "LIUF_MASTER_SEED plus held-out index modulo 2^31;",
      "historical solver deterministic"
    ),
    "historical run_Fig1_model_comparison.R",
    "historical run_mlr7_perm_importance_B1000.R"
  ),
  rng_kind = c(
    "Mersenne-Twister/Inversion/Rejection",
    "none",
    rep("Mersenne-Twister/Inversion/Rejection", 6L)
  ),
  status = c("fixed", "deterministic", rep("fixed", 6L)),
  notes = c(
    "default seed for classifier operations without a historical seed",
    paste(
      "derive_r_seed preserves every accepted default seed while wrapping",
      "upper-bound offsets into the valid R seed domain"
    ),
    paste(
      "fold_seed=derive_r_seed(master_seed, heldout_index);",
      "LDA refit is otherwise deterministic"
    ),
    "set.seed(1); default randomForest ntree/mtry; preserves historical OOB forest",
    "fold_seed=derive_r_seed(master_seed, heldout_index); kernel=radial",
    "fold_seed=derive_r_seed(master_seed, heldout_index); maxit=300/400",
    "set.seed(1); nfolds=5; lambda=lambda.1se; alpha=0/0.5/1",
    paste(
      "RNGkind set explicitly, then set.seed(42) immediately before baseline;",
      "baseline and all feature permutations run in historical serial order;",
      "B=1000"
    )
  )
)
seed_audit_columns <- c(
  "pipeline",
  "step",
  "script_path",
  "language",
  "algorithm",
  "seed",
  "seed_source",
  "rng_kind",
  "status",
  "notes"
)
setcolorder(classifier_seed_audit, seed_audit_columns)
fwrite(
  classifier_seed_audit,
  file.path(cfg$audit_dir, "mofa_classifier_seed_audit.tsv"),
  sep = "\t"
)


classifier_finished_at <- Sys.time()
classifier_runtime <- data.table(
  stage = "40_factor2_classifier.R_historical_serial_B1000",
  started_epoch = as.numeric(classifier_started_at),
  finished_epoch = as.numeric(classifier_finished_at),
  elapsed_seconds = as.numeric(
    difftime(classifier_finished_at, classifier_started_at, units = "secs")
  ),
  permutation_seed = 42L,
  B = permutation_count,
  rng_kind = "Mersenne-Twister/Inversion/Rejection",
  final_model_sha256 = sha256_file(cfg$analysis_model_path),
  raw_input_sha256 = sha256_file(rebuilt_input_path)
)
fwrite(
  classifier_runtime,
  file.path(cfg$result_dir, "factor2_classifier_runtime.tsv"),
  sep = "\t"
)

cat(
  sprintf(
    paste0(
      "Factor2 classifier completed from the selected K15 analysis model; ",
      "model comparison ",
      "diff=%.12g, performance diff=%.12g, permutation status=%s.\n"
    ),
    comparison_diff,
    performance_diff,
    permutation_status
  )
)
