#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

all_arguments <- commandArgs(trailingOnly = FALSE)
this_script <- normalizePath(
  sub("^--file=", "", all_arguments[grepl("^--file=", all_arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(this_script), "config.R"))
config <- ensure_output_dirs(scrna_cm02_config(this_script))
output_dir <- config$figure3_output_root
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

boundary_audit <- file.path(
  config$audit_root,
  "cm02_historical_downstream_boundary.tsv"
)
if (config$downstream_source == "historical" && !file.exists(boundary_audit)) {
  stop("Run 10_materialize_historical_downstream.R first.")
}

cm02_members <- config$cm02_anchor_states
nmf_basis <- NMF::basis(readRDS(config$analysis_downstream_paths$nmf))
normalized_basis <- sweep(nmf_basis, 2L, colSums(nmf_basis), "/")
cm02_weights <- normalized_basis[cm02_members, "CM02"]

tcga_proportions <- fread(
  config$analysis_downstream_paths$tcga_proportions
)
proportion_matrix <- as.matrix(tcga_proportions[, ..cm02_members])
tcga_proportions[, cm02_activity := as.numeric(
  proportion_matrix %*% cm02_weights
)]

clinical <- fread(config$analysis_downstream_paths$tcga_clinical)
clinical[, OS_time := fifelse(
  vital_status == "DECEASED",
  as.numeric(days_to_death),
  as.numeric(days_to_last_followup)
)]
clinical[, OS := as.integer(vital_status == "DECEASED")]
clinical[, stage_num := as.numeric(factor(
  pathologic_stage,
  levels = c(
    "Stage I", "Stage II", "Stage IIIA", "Stage IIIB", "Stage IIIC",
    "Stage IV", "Stage IVA", "Stage IVB"
  )
))]
clinical[, age_num := as.numeric(age_at_initial_pathologic_diagnosis)]
clinical_tumor <- clinical[
  sample_type == "Primary Tumor" & !is.na(OS_time) & OS_time > 0
]
survival_data <- merge(
  clinical_tumor,
  tcga_proportions[, .(sampleID, cm02_activity)],
  by = "sampleID"
)
if (nrow(survival_data) != 365L || sum(survival_data$OS) != 130L) {
  stop("TCGA survival cohort contract failed.")
}

univariable_fit <- coxph(
  Surv(OS_time, OS) ~ cm02_activity,
  data = survival_data
)
adjusted_fit <- coxph(
  Surv(OS_time, OS) ~ cm02_activity + stage_num + age_num + gender,
  data = survival_data
)
univariable_summary <- summary(univariable_fit)
adjusted_summary <- summary(adjusted_fit)
cox_result <- data.table(
  metric = "CM02 activity (cell-freq)",
  HR_uni = univariable_summary$coefficients[
    "cm02_activity",
    "exp(coef)"
  ],
  CI_uni = sprintf(
    "[%.2f, %.2f]",
    exp(confint(univariable_fit))["cm02_activity", 1L],
    exp(confint(univariable_fit))["cm02_activity", 2L]
  ),
  p_uni = univariable_summary$coefficients[
    "cm02_activity",
    "Pr(>|z|)"
  ],
  HR_adj = adjusted_summary$coefficients[
    "cm02_activity",
    "exp(coef)"
  ],
  CI_adj = sprintf(
    "[%.2f, %.2f]",
    exp(confint(adjusted_fit))["cm02_activity", 1L],
    exp(confint(adjusted_fit))["cm02_activity", 2L]
  ),
  p_adj = adjusted_summary$coefficients[
    "cm02_activity",
    "Pr(>|z|)"
  ],
  n = univariable_fit$n,
  events = univariable_fit$nevent,
  n_adj = adjusted_fit$n,
  events_adj = adjusted_fit$nevent
)
fwrite(
  cox_result,
  file.path(output_dir, "Figure3C_TCGA_CM02_activity_Cox.csv")
)

maxstat_cutpoint <- maxstat::maxstat.test(
  Surv(OS_time, OS) ~ cm02_activity,
  data = survival_data,
  smethod = "LogRank",
  pmethod = "none",
  minprop = 0.2,
  maxprop = 0.8
)
cutpoint <- as.numeric(maxstat_cutpoint$estimate)
categorized <- survival_data
categorized$grp <- factor(
  ifelse(categorized$cm02_activity > cutpoint, "high", "low"),
  levels = c("low", "high")
)
km_fit <- survfit(Surv(OS_time, OS) ~ grp, data = categorized)
group_cox <- coxph(Surv(OS_time, OS) ~ grp, data = categorized)
high_vs_low_hr <- unname(exp(coef(group_cox))[[1L]])
high_vs_low_ci <- unname(exp(confint(group_cox))[1L, ])
logrank_p <- 1 - pchisq(
  survdiff(Surv(OS_time, OS) ~ grp, data = categorized)$chisq,
  df = 1L
)
km_statistics <- data.table(
  cutpoint = cutpoint,
  high_vs_low_hr = high_vs_low_hr,
  ci_lower = high_vs_low_ci[[1L]],
  ci_upper = high_vs_low_ci[[2L]],
  logrank_p = logrank_p,
  low_n = sum(categorized$grp == "low"),
  low_events = sum(categorized$grp == "low" & categorized$OS == 1L),
  high_n = sum(categorized$grp == "high"),
  high_events = sum(categorized$grp == "high" & categorized$OS == 1L)
)
fwrite(
  km_statistics,
  file.path(output_dir, "Figure3C_TCGA_CM02_activity_KM_stats.csv")
)


weight_table <- data.table(
  component = names(cm02_weights),
  normalized_weight = as.numeric(cm02_weights)
)
fwrite(
  weight_table,
  file.path(output_dir, "Figure3_CM02_normalized_weights.csv")
)
message(
  sprintf(
    paste0(
      "Figure 3C replay complete: n=%d, events=%d, ",
      "cut=%.8f, HR(high vs low)=%.6f, p=%.12g"
    ),
    nrow(categorized),
    sum(categorized$OS),
    cutpoint,
    high_vs_low_hr,
    logrank_p
  )
)
