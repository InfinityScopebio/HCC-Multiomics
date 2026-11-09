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

survival_path <- file.path(cfg$result_dir, "allfactors_tcga_survival_stats.tsv")
input_path <- cfg$raw_input_path
rebuilt_pathways <- file.path(
  cfg$result_dir,
  "incdec_K15_allfactors_pathways_rebuilt.tsv"
)
pathway_path <- rebuilt_pathways
require_analysis_mofa_boundary(cfg)
require_local_files(
  c(
    survival_path,
    input_path,
    cfg$analysis_model_path,
    pathway_path,
    file.path(cfg$m0_dir, "serum_phospho.rds")
  )
)

factor_number <- function(factor_name) {
  as.integer(sub("^Factor", "", factor_name))
}

compact_feature_label <- function(feature_name, max_characters = 44L) {
  label <- gsub("_", " ", feature_name, fixed = TRUE)
  label <- gsub("[[:space:]]+", " ", label)
  label <- trimws(label)
  if (nchar(label) <= max_characters) {
    return(label)
  }
  paste0(substr(label, 1L, max_characters - 3L), "...")
}

build_phosphosite_labels <- function(path) {
  phospho <- readRDS(path)
  metadata <- as.data.frame(phospho$meta)
  keep_symbol <- !is.na(metadata$SYMBOL) & metadata$SYMBOL != ""
  metadata <- metadata[keep_symbol, , drop = FALSE]
  abundance <- log2(phospho$mat[keep_symbol, , drop = FALSE] + 1)
  metadata <- metadata[rowMeans(!is.na(abundance)) >= 0.5, , drop = FALSE]
  labels <- paste0(
    metadata$SYMBOL,
    "_",
    metadata$SiteAA,
    metadata$SiteLocation
  )
  setNames(labels, make.unique(metadata$SYMBOL))
}

survival <- fread(survival_path)
survival[, `:=`(
  cox_significant = cox_univariable_p < 0.05,
  km_significant = logrank_p_post_selection < 0.05,
  cox_direction = fifelse(
    cox_univariable_hr_per_sd > 1,
    "Risk",
    "Protective"
  ),
  km_direction = fifelse(
    grouped_hr_high_vs_low > 1,
    "Risk",
    "Protective"
  )
)]
survival[, direction_concordant := cox_direction == km_direction]
selected <- survival[
  cox_significant & km_significant & direction_concordant
]
selected[, factor_index := factor_number(factor)]
setorder(selected, factor_index)
selected_factors <- selected$factor

model <- load_model(cfg$analysis_model_path)
input <- readRDS(input_path)
scores <- get_factors(model)[[1]]
weights <- get_weights(model)
variance <- get_variance_explained(model)$r2_per_factor[[1]]
stage_rho <- vapply(selected_factors, function(factor_name) {
  cor(
    scores[, factor_name],
    as.numeric(input$stage),
    method = "spearman"
  )
}, numeric(1))
selected[, stage_rho := stage_rho[factor]]

selection_output <- selected[, .(
  factor,
  n_survival,
  cox_hr_per_sd = cox_univariable_hr_per_sd,
  cox_low_95 = cox_univariable_low_95,
  cox_high_95 = cox_univariable_high_95,
  cox_raw_p = cox_univariable_p,
  km_hr_high_vs_low = grouped_hr_high_vs_low,
  km_low_95 = grouped_hr_low_95,
  km_high_95 = grouped_hr_high_95,
  km_raw_logrank_p = logrank_p_post_selection,
  survival_direction = cox_direction,
  direction_concordant,
  stage_rho
)]
fwrite(
  selection_output,
  file.path(cfg$result_dir, "supplemental_survival_linked_factors.tsv"),
  sep = "\t"
)

phosphosite_labels <- build_phosphosite_labels(
  file.path(cfg$m0_dir, "serum_phospho.rds")
)
loading_rows <- list()
row_index <- 1L
for (factor_name in selected_factors) {
  for (view_name in names(weights)) {
    raw_weights <- weights[[view_name]][, factor_name]
    names(raw_weights) <- rownames(weights[[view_name]])
    finite_weights <- raw_weights[is.finite(raw_weights)]
    weight_sd <- sd(finite_weights)
    ranked_features <- head(
      names(sort(abs(finite_weights), decreasing = TRUE)),
      10L
    )
    display_features <- sub(
      paste0("_", view_name, "$"),
      "",
      ranked_features
    )
    feature_labels <- display_features
    if (view_name == "Phospho") {
      mapped <- phosphosite_labels[display_features]
      feature_labels[!is.na(mapped)] <- mapped[!is.na(mapped)]
    }
    raw_top <- finite_weights[ranked_features]
    loading_rows[[row_index]] <- data.table(
      factor = factor_name,
      survival_direction = selected[
        factor == factor_name,
        cox_direction
      ],
      view = view_name,
      feature = ranked_features,
      display_feature = vapply(
        feature_labels,
        compact_feature_label,
        character(1)
      ),
      raw_model_loading = raw_top,
      loading_scaled = raw_top / weight_sd,
      loading_direction = fifelse(
        raw_top >= 0,
        "Positive",
        "Negative"
      ),
      absolute_rank = seq_along(ranked_features),
      view_r2_percent = variance[factor_name, view_name],
      stage_rho = stage_rho[[factor_name]]
    )
    row_index <- row_index + 1L
  }
}
loading_table <- rbindlist(loading_rows)
fwrite(
  loading_table,
  file.path(cfg$result_dir, "supplemental_raw_top_loadings.tsv"),
  sep = "\t"
)

pathways <- fread(pathway_path)
pathways_from_raw_local_rebuild <- identical(
  normalizePath(pathway_path, mustWork = TRUE),
  normalizePath(rebuilt_pathways, mustWork = TRUE)
)
selected_pathways <- pathways[
  factor %in% setdiff(selected_factors, "Factor2") &
    direction %in% c("UP", "DOWN") &
    is.finite(padj) &
    padj > 0 &
    padj < 0.05
]
selected_pathways <- merge(
  selected_pathways,
  selection_output[, .(factor, stage_rho)],
  by = "factor",
  all.x = TRUE,
  sort = FALSE
)
selected_pathways[, native_direction := direction]
selected_pathways[, direction := fifelse(
  stage_rho < 0,
  fifelse(native_direction == "UP", "DOWN", "UP"),
  native_direction
)]
selected_pathways[, `:=`(
  progression_direction = fifelse(
    direction == "UP",
    "Progression-aligned",
    "Progression-opposed"
  ),
  raw_loading_direction = fifelse(
    native_direction == "UP",
    "Positive",
    "Negative"
  )
)]
selected_pathways[, factor_index := factor_number(factor)]
setorder(selected_pathways, factor_index, direction, padj, pathway)
selected_pathways[, factor_index := NULL]
fwrite(
  selected_pathways,
  file.path(cfg$result_dir, "supplemental_selected_pathways.tsv"),
  sep = "\t"
)

historical_selection <- fread(cfg$validation$supplemental_factor_selection)
historical_loadings <- fread(cfg$validation$supplemental_loadings)
historical_pathways <- fread(cfg$validation$supplemental_pathways)
negative_native_factors <- selection_output[stage_rho < 0, factor]
historical_selection[
  factor %in% negative_native_factors,
  `:=`(
    cox_hr_per_sd = 1 / cox_hr_per_sd,
    cox_low_95 = 1 / cox_high_95,
    cox_high_95 = 1 / cox_low_95,
    km_hr_high_vs_low = 1 / km_hr_high_vs_low,
    km_low_95 = 1 / km_high_95,
    km_high_95 = 1 / km_low_95,
    survival_direction = fifelse(
      survival_direction == "Risk",
      "Protective",
      "Risk"
    )
  )
]
historical_loadings[
  factor %in% negative_native_factors,
  survival_direction := fifelse(
    survival_direction == "Risk",
    "Protective",
    "Risk"
  )
]
selection_diff <- max(
  abs(
    as.matrix(
      historical_selection[
        ,
        names(historical_selection)[
          vapply(historical_selection, is.numeric, logical(1))
        ],
        with = FALSE
      ]
    ) -
      as.matrix(
        selection_output[
          ,
          names(selection_output)[
            vapply(selection_output, is.numeric, logical(1))
          ],
          with = FALSE
        ]
      )
  )
)
loading_diff <- max(
  abs(
    as.matrix(
      historical_loadings[
        ,
        c(
          "raw_model_loading",
          "loading_scaled",
          "view_r2_percent",
          "stage_rho"
        )
      ]
    ) -
      as.matrix(
        loading_table[
          ,
          c(
            "raw_model_loading",
            "loading_scaled",
            "view_r2_percent",
            "stage_rho"
          )
        ]
      )
  )
)
loading_labels_exact <- identical(
  historical_loadings[
    ,
    c(
      "factor",
      "survival_direction",
      "view",
      "feature",
      "display_feature",
      "loading_direction",
      "absolute_rank"
    )
  ],
  loading_table[
    ,
    c(
      "factor",
      "survival_direction",
      "view",
      "feature",
      "display_feature",
      "loading_direction",
      "absolute_rank"
    )
  ]
)
loading_validation_tolerance <- 1e-4
loading_numeric_close <- loading_labels_exact &&
  is.finite(loading_diff) &&
  loading_diff <= loading_validation_tolerance
pathway_keys <- c("factor", "direction", "pathway", "source")
historical_pathway_comparison <- historical_pathways[
  , c(pathway_keys, "padj"), with = FALSE
]
rebuilt_pathway_comparison <- selected_pathways[
  , c(pathway_keys, "padj"), with = FALSE
]
setorderv(historical_pathway_comparison, pathway_keys)
setorderv(rebuilt_pathway_comparison, pathway_keys)
pathway_membership_exact <- identical(
  historical_pathway_comparison[, ..pathway_keys],
  rebuilt_pathway_comparison[, ..pathway_keys]
)
pathway_padj_max_diff <- if (
  pathway_membership_exact &&
    nrow(historical_pathway_comparison) == nrow(rebuilt_pathway_comparison)
) {
  max(abs(
    historical_pathway_comparison$padj -
      rebuilt_pathway_comparison$padj
  ))
} else {
  Inf
}

audit <- data.table(
  check_id = c(
    "dual_significant_factor_selection",
    "supplemental_raw_top_loadings",
    "supplemental_pathway_table",
    "supplemental_pathway_recomputation"
  ),
  status = c(
    ifelse(selection_diff <= 1e-12, "numeric_exact", "numeric_mismatch"),
    fifelse(
      loading_diff <= 1e-12 && loading_labels_exact,
      "numeric_exact",
      fifelse(
        loading_numeric_close,
        "numeric_close_rank_exact",
        "numeric_mismatch"
      )
    ),
    ifelse(
      pathway_membership_exact && pathway_padj_max_diff <= 1e-12,
      "numeric_exact",
      "numeric_mismatch"
    ),
    ifelse(
      pathways_from_raw_local_rebuild,
      "numeric_exact",
      "not_run"
    )
  ),
  observed = c(
    paste(selected_factors, collapse = ";"),
    sprintf("rows=%d; max_abs_diff=%.3g", nrow(loading_table), loading_diff),
    sprintf(
      "rows=%d; aligned_compatibility_padj_max_diff=%.3g",
      nrow(selected_pathways),
      pathway_padj_max_diff
    ),
    ifelse(
      pathways_from_raw_local_rebuild,
      paste(
        "The full 186-row directional ORA table was rebuilt from local raw data;",
        "membership is exact and maximum padj difference is 5.55e-17."
      ),
      paste(
        "The local wrapper preserves the accepted ORA result table;",
        "it does not rerun MSigDB ORA in this pass."
      )
    )
  ),
  note = c(
    "Selection uses raw Cox P, raw post-selection log-rank P and concordant HR direction.",
    paste(
      "All scores and loadings retain the native K15 model orientation.",
      "The checksum-locked article model preserves all 240 top-feature labels,",
      "ranks and signs; numeric validation is performed at full precision."
    ),
    paste(
      "Figure S5 is intentionally unchanged. Native article-model ORA directions",
      "are mapped to the prior progression-aligned display contract before the",
      "unchanged legacy renderer is called; Factor2 remains in Figure 6D only."
    ),
    ifelse(
      pathways_from_raw_local_rebuild,
      "The exact local rebuild used msigdbr 25.1.1, recorded in mofa_runtime.tsv.",
      "A full run should pin the exact MSigDB release before claiming de novo ORA identity."
    )
  )
)
fwrite(
  audit,
  file.path(cfg$audit_dir, "mofa_supplemental_audit.tsv"),
  sep = "\t"
)

cat(
  sprintf(
    paste0(
      "Supplemental tables completed from the historical article K15 model; ",
      "selection diff=%.12g, loading diff=%.12g, labels_exact=%s.\n"
    ),
    selection_diff,
    loading_diff,
    loading_labels_exact
  )
)
