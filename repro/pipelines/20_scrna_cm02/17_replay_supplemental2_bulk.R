#!/usr/bin/env Rscript

# Compute within-cohort bulk trajectory validation tables.

suppressPackageStartupMessages({
  library(dplyr)
})

script_path <- normalizePath(
  sub("^--file=", "", grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "config.R"))
config <- scrna_cm02_config(script_path)

output_dir <- file.path(
  config$result_root,
  "supplemental_figure2"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- config$analysis_downstream_paths$bulk_pseudotime
if (!file.exists(input_path)) {
  stop("Frozen bulk embedding is missing: ", input_path)
}

embedding <- read.csv(
  input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_columns <- c(
  "sampleID", "Group", "Cohort", "Pseudotime", "phate1", "phate2"
)
if (!all(required_columns %in% names(embedding))) {
  stop("Frozen bulk embedding is missing required columns.")
}


format_p <- function(value) {
  if (!is.finite(value)) {
    return("NA")
  }
  if (value < 0.001) {
    return(formatC(value, format = "e", digits = 2))
  }
  formatC(value, format = "f", digits = 3)
}


# Supplemental Figure 2D: within-cohort disease-distance validation.
bulk_within <- embedding %>%
  filter(
    (Cohort == "GSE63898" & Group %in% c("Cirrhosis", "HCC")) |
      (Cohort == "TCGA" & Group %in% c("Adjacent", "HCC"))
  ) %>%
  mutate(
    comparison = ifelse(
      Cohort == "GSE63898",
      "GSE63898: Cirrhosis vs HCC",
      "TCGA: Adjacent vs HCC"
    ),
    Group = factor(
      Group,
      levels = c("Adjacent", "Cirrhosis", "HCC")
    )
  )

comparison_specs <- list(
  list(
    cohort = "GSE63898",
    reference = "Cirrhosis",
    comparison = "GSE63898: Cirrhosis vs HCC"
  ),
  list(
    cohort = "TCGA",
    reference = "Adjacent",
    comparison = "TCGA: Adjacent vs HCC"
  )
)

bulk_stats <- bind_rows(lapply(comparison_specs, function(spec) {
  cohort_data <- bulk_within %>%
    filter(Cohort == spec$cohort)
  hcc_values <- cohort_data$Pseudotime[cohort_data$Group == "HCC"]
  reference_values <- cohort_data$Pseudotime[
    cohort_data$Group == spec$reference
  ]
  test_result <- wilcox.test(
    hcc_values,
    reference_values,
    exact = FALSE
  )
  probability_superiority <- unname(test_result$statistic) /
    (length(hcc_values) * length(reference_values))

  tibble(
    comparison = spec$comparison,
    cohort = spec$cohort,
    reference_group = spec$reference,
    n_reference = length(reference_values),
    n_hcc = length(hcc_values),
    median_reference = median(reference_values),
    median_hcc = median(hcc_values),
    p_value = test_result$p.value,
    rank_biserial = 2 * probability_superiority - 1
  )
})) %>%
  mutate(
    bh_q_value = p.adjust(p_value, method = "BH"),
    annotation = sprintf(
      "n=%d vs %d; BH q=%s; rank-biserial r=%.2f",
      n_reference,
      n_hcc,
      vapply(bh_q_value, format_p, character(1)),
      rank_biserial
    )
  )

bulk_summary <- bulk_within %>%
  group_by(comparison, Cohort, Group) %>%
  summarise(
    n_samples = n(),
    median = median(Pseudotime),
    q25 = quantile(Pseudotime, 0.25),
    q75 = quantile(Pseudotime, 0.75),
    .groups = "drop"
  )

write.csv(
  bulk_within,
  file.path(output_dir, "S2D_bulk_within_cohort_values.csv"),
  row.names = FALSE
)
write.csv(
  bulk_summary,
  file.path(output_dir, "S2D_bulk_within_cohort_summary.csv"),
  row.names = FALSE
)
write.csv(
  bulk_stats,
  file.path(output_dir, "S2D_bulk_within_cohort_statistics.csv"),
  row.names = FALSE
)

# Supplemental Figure 2E-F: frozen-coordinate UMAP views.
embedding_umap <- embedding %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Healthy", "Adjacent", "Cirrhosis", "HCC")
    )
  )

padded_limits <- function(values, padding_fraction) {
  value_range <- range(values, na.rm = TRUE)
  padding <- diff(value_range) * padding_fraction
  value_range + c(-padding, padding)
}

x_limits <- padded_limits(embedding_umap$phate1, 0.04)
y_limits <- padded_limits(embedding_umap$phate2, 0.04)
coordinate_contract <- data.frame(
  axis = c("x", "y"),
  minimum = c(x_limits[1], y_limits[1]),
  maximum = c(x_limits[2], y_limits[2]),
  ratio = c(1, 1),
  n_samples = c(nrow(embedding_umap), nrow(embedding_umap)),
  stringsAsFactors = FALSE
)
write.table(
  coordinate_contract,
  file.path(output_dir, "S2EF_coordinate_contract.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  sprintf(
    paste0(
      "Computed bulk validation from %d samples; x=[%.3f, %.3f], ",
      "y=[%.3f, %.3f]."
    ),
    nrow(embedding_umap),
    x_limits[1],
    x_limits[2],
    y_limits[1],
    y_limits[2]
  )
)
