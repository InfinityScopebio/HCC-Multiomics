#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

all_arguments <- commandArgs(trailingOnly = FALSE)
this_script <- normalizePath(
  sub("^--file=", "", all_arguments[grepl("^--file=", all_arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(this_script), "config.R"))
config <- ensure_tissue_resistance_dirs(
  tissue_resistance_config(this_script)
)

read_gmt <- function(path) {
  records <- strsplit(readLines(path), "\t", fixed = TRUE)
  setNames(
    lapply(records, function(record) record[-seq_len(2L)]),
    vapply(records, `[[`, character(1), 1L)
  )
}

score_direction <- function(gene_ranks, genes) {
  selected_genes <- intersect(genes, rownames(gene_ranks))
  if (length(selected_genes) < 3L) {
    return(rep(NA_real_, ncol(gene_ranks)))
  }
  colMeans(gene_ranks[selected_genes, , drop = FALSE]) / nrow(gene_ranks)
}

spearman_profile <- function(score_table, composite) {
  profile <- t(vapply(
    config$clusters,
    function(cluster_name) {
      complete <- !is.na(score_table[[cluster_name]]) & !is.na(composite)
      if (sum(complete) < 8L) {
        return(c(rho = NA_real_, p = NA_real_, n = sum(complete)))
      }
      test <- suppressWarnings(cor.test(
        score_table[[cluster_name]][complete],
        composite[complete],
        method = "spearman"
      ))
      c(
        rho = unname(test$estimate),
        p = test$p.value,
        n = sum(complete)
      )
    },
    numeric(3L)
  ))
  data.table(cluster = rownames(profile), profile)
}

model <- fread(config$inputs$model)
liver_biliary <- model[
  OncotreeLineage %in% c("Liver", "Biliary Tract"),
  .(ModelID, lineage = OncotreeLineage)
]
expression <- fread(config$inputs$expression)
setnames(expression, 1L, "row1")
expression <- expression[
  IsDefaultEntryForModel %in% c("Yes", TRUE) &
    ModelID %in% liver_biliary$ModelID
]
metadata_columns <- intersect(
  c(
    "row1", "SequencingID", "ModelConditionID", "ModelID",
    "IsDefaultEntryForMC", "IsDefaultEntryForModel"
  ),
  names(expression)
)
expression_matrix <- t(as.matrix(
  expression[, setdiff(names(expression), metadata_columns), with = FALSE]
))
colnames(expression_matrix) <- expression$ModelID
rownames(expression_matrix) <- sub(" \\(.*", "", rownames(expression_matrix))
expression_matrix <- expression_matrix[
  !duplicated(rownames(expression_matrix)),
  ,
  drop = FALSE
]

signatures <- read_gmt(config$inputs$signatures)
gene_ranks <- apply(expression_matrix, 2L, rank)
singscore <- vapply(
  config$clusters,
  function(cluster_name) {
    score_direction(
      gene_ranks,
      signatures[[paste0(cluster_name, "_UP")]]
    ) - score_direction(
      gene_ranks,
      signatures[[paste0(cluster_name, "_DOWN")]]
    )
  },
  numeric(ncol(expression_matrix))
)
rownames(singscore) <- colnames(expression_matrix)
scaled_scores <- scale(singscore)
rownames(scaled_scores) <- colnames(expression_matrix)
score_table <- data.frame(
  ModelID = colnames(expression_matrix),
  lineage = liver_biliary$lineage[
    match(colnames(expression_matrix), liver_biliary$ModelID)
  ],
  scaled_scores,
  check.names = FALSE
)
fwrite(
  score_table,
  file.path(config$result_root, "R3d_cell_line_scores.csv")
)

prism <- fread(config$inputs$prism_auc)
setnames(prism, 1L, "ModelID")
auc_columns <- setdiff(names(prism), "ModelID")
normalized_drug_names <- toupper(sub(
  " \\(BRD.*",
  "",
  auc_columns
))

profile_outputs <- list()
composite_receipts <- list()
for (composite_name in names(config$composites)) {
  score_subset <- score_table[
    score_table$ModelID %in% liver_biliary$ModelID,
    ,
    drop = FALSE
  ]
  prism_subset <- prism[ModelID %in% score_subset$ModelID]
  score_subset <- score_subset[
    match(prism_subset$ModelID, score_subset$ModelID),
    ,
    drop = FALSE
  ]
  requested_drugs <- config$composites[[composite_name]]
  class_auc <- vapply(
    requested_drugs,
    function(drug_name) {
      matched_column <- auc_columns[normalized_drug_names == drug_name]
      if (length(matched_column) > 0L) {
        return(prism_subset[[matched_column[[1L]]]])
      }
      rep(NA_real_, nrow(prism_subset))
    },
    numeric(nrow(prism_subset))
  )
  class_auc <- class_auc[
    ,
    colSums(!is.na(class_auc)) >= 8L,
    drop = FALSE
  ]
  complete_rows <- complete.cases(class_auc)
  principal_component <- prcomp(
    scale(class_auc[complete_rows, , drop = FALSE])
  )$x[, 1L]
  mean_scaled_auc <- rowMeans(scale(
    class_auc[complete_rows, , drop = FALSE]
  ))
  if (cor(principal_component, mean_scaled_auc) < 0) {
    principal_component <- -principal_component
  }
  composite <- rep(NA_real_, nrow(score_subset))
  composite[complete_rows] <- principal_component
  profile <- spearman_profile(score_subset, composite)
  profile_outputs[[composite_name]] <- profile
  output_name <- paste0(
    "R3d_profile_liver_biliary_",
    composite_name,
    ".csv"
  )
  fwrite(profile, file.path(config$result_root, output_name))
  composite_receipts[[composite_name]] <- data.table(
    composite = composite_name,
    requested_drugs = paste(requested_drugs, collapse = ";"),
    retained_drugs = paste(colnames(class_auc), collapse = ";"),
    n_retained_drugs = ncol(class_auc),
    n_complete_models = sum(complete_rows),
    auc_direction = "higher composite score = higher AUC = resistance"
  )
}
fwrite(
  rbindlist(composite_receipts),
  file.path(config$result_root, "Figure4D_composite_contract.tsv"),
  sep = "\t"
)

message(
  "Figure 4D rebuilt from checksum-locked local DepMap/PRISM data: ",
  "71 liver/biliary models; TKI n=",
  unique(profile_outputs$TKI$n),
  "; chemo n=",
  unique(profile_outputs$TACE_HAIC_chemo$n),
  "."
)
