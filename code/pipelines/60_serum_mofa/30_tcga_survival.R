#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
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

clinical_path <- file.path(cfg$tcga_raw_dir, "LIHC_clinicalMatrix.tsv")
rebuilt_axis_store <- file.path(
  cfg$result_dir,
  "incdec_K15_all_axes_rebuilt.rds"
)
require_local_files(
  c(
    rebuilt_axis_store,
    cfg$validation$survival_stats,
    clinical_path
  )
)

factor_number <- function(factor_name) {
  as.integer(sub("^Factor", "", factor_name))
}

prepare_survival_data <- function(factor_name, axis_store) {
  data <- as.data.table(copy(axis_store[[factor_name]]$tcga_df))
  data <- data[
    is.finite(OS_time) & OS_time > 0 &
      !is.na(OS) &
      is.finite(score)
  ]
  if (nrow(data) != 365L) {
    stop(factor_name, " expected 365 survival cases, found ", nrow(data))
  }
  data
}

fit_dichotomized_survival <- function(data, cutpoint) {
  split_data <- copy(data)
  split_data[, group := factor(
    ifelse(score > cutpoint, "High", "Low"),
    levels = c("High", "Low")
  )]
  fit <- survfit(Surv(OS_time, OS) ~ group, data = split_data)
  logrank <- survdiff(Surv(OS_time, OS) ~ group, data = split_data)
  logrank_p <- pchisq(logrank$chisq, df = 1L, lower.tail = FALSE)

  split_data[, group_hr := factor(group, levels = c("Low", "High"))]
  grouped_cox <- summary(coxph(
    Surv(OS_time, OS) ~ group_hr,
    data = split_data
  ))
  grouped_ci <- grouped_cox$conf.int[1L, ]
  group_summary <- split_data[, .(n = .N, events = sum(OS)), by = group]
  setkey(group_summary, group)
  medians <- as.data.table(
    summary(fit)$table,
    keep.rownames = "strata"
  )
  medians <- medians[is.finite(median), .(
    group = sub("^group=", "", strata),
    median
  )]

  list(
    data = split_data,
    fit = fit,
    medians = medians,
    logrank_p = logrank_p,
    grouped_hr = unname(grouped_ci["exp(coef)"]),
    grouped_low_95 = unname(grouped_ci["lower .95"]),
    grouped_high_95 = unname(grouped_ci["upper .95"]),
    n_high = group_summary["High"]$n,
    events_high = group_summary["High"]$events,
    n_low = group_summary["Low"]$n,
    events_low = group_summary["Low"]$events
  )
}

fit_factor <- function(factor_name, axis_store) {
  data <- prepare_survival_data(factor_name, axis_store)
  cutpoint_object <- maxstat::maxstat.test(
    Surv(OS_time, OS) ~ score,
    data = as.data.frame(data),
    smethod = "LogRank",
    pmethod = "none",
    minprop = 0.2,
    maxprop = 0.8
  )
  cutpoint <- as.numeric(cutpoint_object$estimate)
  maxstat_result <- maxstat::maxstat.test(
    Surv(OS_time, OS) ~ score,
    data = as.data.frame(data),
    smethod = "LogRank",
    pmethod = "Lau94",
    minprop = 0.2,
    maxprop = 0.8
  )
  if (!isTRUE(all.equal(
    cutpoint,
    as.numeric(maxstat_result$estimate),
    tolerance = 1e-10
  ))) {
    stop("Cutpoint disagreement for ", factor_name)
  }

  selected_split <- fit_dichotomized_survival(data, cutpoint)
  median_cutpoint <- median(data$score)
  median_split <- fit_dichotomized_survival(data, median_cutpoint)
  continuous_cox <- summary(coxph(Surv(OS_time, OS) ~ score, data = data))
  continuous_ci <- continuous_cox$conf.int["score", ]

  statistics <- data.table(
    factor = factor_name,
    n_survival = nrow(data),
    events = sum(data$OS),
    cutpoint = cutpoint,
    n_high = selected_split$n_high,
    events_high = selected_split$events_high,
    n_low = selected_split$n_low,
    events_low = selected_split$events_low,
    logrank_p_post_selection = selected_split$logrank_p,
    maxstat_lau94_p = min(1, as.numeric(maxstat_result$p.value)),
    grouped_hr_high_vs_low = selected_split$grouped_hr,
    grouped_hr_low_95 = selected_split$grouped_low_95,
    grouped_hr_high_95 = selected_split$grouped_high_95,
    median_cutpoint = median_cutpoint,
    n_high_median = median_split$n_high,
    events_high_median = median_split$events_high,
    n_low_median = median_split$n_low,
    events_low_median = median_split$events_low,
    logrank_p_median_split = median_split$logrank_p,
    grouped_hr_high_vs_low_median = median_split$grouped_hr,
    grouped_hr_median_low_95 = median_split$grouped_low_95,
    grouped_hr_median_high_95 = median_split$grouped_high_95,
    cox_univariable_hr_per_sd = unname(continuous_ci["exp(coef)"]),
    cox_univariable_low_95 = unname(continuous_ci["lower .95"]),
    cox_univariable_high_95 = unname(continuous_ci["upper .95"]),
    cox_univariable_p = unname(
      continuous_cox$coefficients["score", "Pr(>|z|)"]
    ),
    median_high_days = selected_split$medians[group == "High"]$median,
    median_low_days = selected_split$medians[group == "Low"]$median,
    median_split_high_days = median_split$medians[group == "High"]$median,
    median_split_low_days = median_split$medians[group == "Low"]$median
  )
  list(
    statistics = statistics,
    selected_split = selected_split,
    median_split = median_split
  )
}

axis_store <- readRDS(rebuilt_axis_store)
factor_names <- names(axis_store)
factor_names <- factor_names[order(factor_number(factor_names))]
if (!all(vapply(
  axis_store[factor_names],
  function(axis) identical(axis$score_orientation, "native"),
  logical(1)
))) {
  stop("TCGA survival requires native-orientation factor projections.")
}
results <- setNames(
  lapply(factor_names, fit_factor, axis_store = axis_store),
  factor_names
)
statistics <- rbindlist(lapply(results, function(result) result$statistics))
statistics[, logrank_q_bh := p.adjust(logrank_p_post_selection, method = "BH")]
statistics[, maxstat_lau94_q_bh := p.adjust(maxstat_lau94_p, method = "BH")]
statistics[, logrank_median_q_bh := p.adjust(
  logrank_p_median_split,
  method = "BH"
)]
statistics[, cox_univariable_q_bh := p.adjust(
  cox_univariable_p,
  method = "BH"
)]

output_stats <- file.path(cfg$result_dir, "allfactors_tcga_survival_stats.tsv")
fwrite(statistics, output_stats, sep = "\t")

clinical <- fread(clinical_path)
clinical[, OS_time := ifelse(
  vital_status == "DECEASED",
  as.numeric(days_to_death),
  as.numeric(days_to_last_followup)
)]
clinical[, OS := as.integer(vital_status == "DECEASED")]
clinical_survival <- clinical[
  sample_type == "Primary Tumor" &
    is.finite(OS_time) &
    OS_time > 0,
  .(sampleID, OS_time, OS)
]
axis_survival_by_factor <- setNames(lapply(factor_names, function(factor_name) {
  factor_survival <- prepare_survival_data(factor_name, axis_store)[
    ,
    .(sampleID, OS_time, OS)
  ]
  setorder(factor_survival, sampleID)
  factor_survival
}), factor_names)
axis_survival <- axis_survival_by_factor[["Factor1"]]
factor_cohort_exact <- all(vapply(
  axis_survival_by_factor,
  function(factor_survival) identical(factor_survival, axis_survival),
  logical(1)
))
setorder(clinical_survival, sampleID)
matched_survival <- merge(
  axis_survival,
  clinical_survival,
  by = "sampleID",
  suffixes = c("_axis", "_clinical"),
  all.x = TRUE
)
clinical_subset_exact <- nrow(matched_survival) == nrow(axis_survival) &&
  all(is.finite(matched_survival$OS_time_clinical)) &&
  identical(matched_survival$OS_time_axis, matched_survival$OS_time_clinical) &&
  identical(matched_survival$OS_axis, matched_survival$OS_clinical)
clinical_exclusions <- clinical_survival[
  !axis_survival,
  on = "sampleID"
]
fwrite(
  clinical_exclusions,
  file.path(cfg$result_dir, "tcga_survival_projection_exclusions.tsv"),
  sep = "\t"
)

# Validation target is loaded only after statistics are recomputed from the
# production axis store.
historical <- fread(cfg$validation$survival_stats)
setorder(historical, factor)
setorder(statistics, factor)
same_columns <- identical(names(historical), names(statistics))
factor_exact <- identical(historical$factor, statistics$factor)
cohort_exact <- identical(
  historical[, .(n_survival, events)],
  statistics[, .(n_survival, events)]
)
p_value_columns <- c(
  "logrank_p_post_selection",
  "maxstat_lau94_p",
  "cox_univariable_p",
  "logrank_q_bh",
  "maxstat_lau94_q_bh",
  "cox_univariable_q_bh"
)
p_value_max_diff <- max(vapply(
  p_value_columns,
  function(column) {
    max(abs(historical[[column]] - statistics[[column]]), na.rm = TRUE)
  },
  numeric(1)
))
hr_columns <- c(
  "grouped_hr_high_vs_low",
  "cox_univariable_hr_per_sd"
)
hr_magnitude_max_diff <- max(vapply(
  hr_columns,
  function(column) {
    max(
      abs(
        abs(log(historical[[column]])) -
          abs(log(statistics[[column]]))
      ),
      na.rm = TRUE
    )
  },
  numeric(1)
))
ci_pairs <- list(
  c("grouped_hr_low_95", "grouped_hr_high_95"),
  c("cox_univariable_low_95", "cox_univariable_high_95")
)
ci_magnitude_max_diff <- max(vapply(
  ci_pairs,
  function(columns) {
    max(vapply(
      seq_len(nrow(statistics)),
      function(row_index) {
        historical_ci <- sort(abs(log(as.numeric(unlist(
          historical[row_index, ..columns],
          use.names = FALSE
        )))))
        native_ci <- sort(abs(log(as.numeric(unlist(
          statistics[row_index, ..columns],
          use.names = FALSE
        )))))
        max(abs(historical_ci - native_ci))
      },
      numeric(1)
    ))
  },
  numeric(1)
))
unordered_pair_exact <- function(first_column, second_column) {
  all(vapply(
    seq_len(nrow(statistics)),
    function(row_index) {
      historical_values <- sort(c(
        historical[[first_column]][row_index],
        historical[[second_column]][row_index]
      ))
      native_values <- sort(c(
        statistics[[first_column]][row_index],
        statistics[[second_column]][row_index]
      ))
      isTRUE(all.equal(
        historical_values,
        native_values,
        tolerance = 0,
        check.attributes = FALSE
      ))
    },
    logical(1)
  ))
}
best_split_exact <- all(c(
  unordered_pair_exact("n_high", "n_low"),
  unordered_pair_exact("events_high", "events_low"),
  unordered_pair_exact("median_high_days", "median_low_days")
))
numeric_max_diff <- max(
  p_value_max_diff,
  hr_magnitude_max_diff,
  ci_magnitude_max_diff
)
statistics_exact <- same_columns &&
  factor_exact &&
  cohort_exact &&
  best_split_exact &&
  is.finite(numeric_max_diff) &&
  numeric_max_diff <= 1e-12

audit <- data.table(
  check_id = c(
    "local_clinical_vs_axis_store",
    "allfactor_survival_statistics",
    "cutpoint_inference",
    "survival_scope"
  ),
  status = c(
    ifelse(
      clinical_subset_exact && factor_cohort_exact,
      "exact_subset",
      "mismatch"
    ),
    ifelse(statistics_exact, "numeric_exact", "numeric_mismatch"),
    "post_selection_with_adjusted_audit",
    ifelse(
      all(statistics$n_survival == 365L) &&
        all(statistics$events == 130L),
      "exact",
      "mismatch"
    )
  ),
  source_model_path = cfg$analysis_model_path,
  source_model_sha256 = sha256_file(cfg$analysis_model_path),
  source_input_sha256 = sha256_file(cfg$raw_input_path),
  observed = c(
    paste(
      nrow(axis_survival),
      "projected cases;",
      nrow(clinical_exclusions),
      "clinical-only exclusions; all_factor_cohorts_exact=",
      factor_cohort_exact
    ),
    sprintf("max_abs_diff=%.3g", numeric_max_diff),
    "maximally selected log-rank cutpoint with Lau94 and BH values retained",
    paste(unique(statistics$n_survival), "cases,", unique(statistics$events), "events")
  ),
  note = c(
    paste(
      "All 15 factors share identical sample IDs, OS times and events; the 365",
      "projected cases match local TCGA clinical data exactly, and six otherwise",
      "survival-eligible cases lack a factor projection."
    ),
    paste(
      "Production statistics were recomputed from the local native-factor score store.",
      "No stage or progression sign is applied. The aligned frozen target is",
      "checked only through orientation-invariant P values, absolute log-HR/CI",
      "magnitudes, and unordered best-split counts and medians."
    ),
    "The displayed best-cutoff log-rank P is optimistic after cutpoint search.",
    "All 15 factors use the same complete-case TCGA-LIHC cohort."
  )
)
fwrite(
  audit,
  file.path(cfg$audit_dir, "mofa_tcga_survival_audit.tsv"),
  sep = "\t"
)


cat(
  sprintf(
    paste0(
      "TCGA all-factor survival recomputed from historical-model native axes; ",
      "validation max_abs_diff=%.12g, statistics_exact=%s, ",
      "same_columns=%s, factor_exact=%s, cohort_exact=%s, ",
      "best_split_exact=%s, factor_cohort_exact=%s, ",
      "clinical_subset_exact=%s.\n"
    ),
    numeric_max_diff,
    statistics_exact,
    same_columns,
    factor_exact,
    cohort_exact,
    best_split_exact,
    factor_cohort_exact,
    clinical_subset_exact
  )
)
