#!/usr/bin/env Rscript

# Compute composition and ratio-of-enrichment tables from the frozen curated atlas.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Seurat)
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

atlas_path <- config$curated_atlas_paths$final
if (!file.exists(atlas_path)) {
  stop("Frozen curated atlas is missing: ", atlas_path)
}
output_dir <- file.path(config$result_root, "figure2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

atlas <- readRDS(atlas_path)
metadata <- atlas@meta.data
metadata$celltype <- as.character(metadata$celltype)
metadata$sub_celltype <- as.character(metadata$sub_celltype)
metadata$group <- factor(
  metadata$group,
  levels = c("Healthy", "Adjacent", "Cirrhosis", "HCC")
)
if (any(as.character(metadata$sample) %in% config$prohibited_samples)) {
  stop("CC1-CC3 entered the frozen Figure 2 atlas.")
}

pretty_state <- function(values) {
  sub("M_$", "Mφ", values)
}
metadata$sub_celltype <- pretty_state(metadata$sub_celltype)

composition_table <- function(subset_metadata, column, levels, panel) {
  subset_metadata <- subset_metadata[
    subset_metadata[[column]] %in% levels,
    ,
    drop = FALSE
  ]
  subset_metadata[[column]] <- factor(
    subset_metadata[[column]],
    levels = levels
  )
  result <- subset_metadata %>%
    count(.data$group, .data[[column]]) %>%
    group_by(.data$group) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  result$panel <- panel
  result
}

major_levels <- intersect(c("T", "NK", "B", "Plasma", "Macrophage", "cDC1", "cDC2", "pDC", "Neutrophil", "Mast", "Hepatocyte", "Cholangiocyte", "Endothelial", "LSEC", "Fibroblast"), unique(metadata$celltype))
t_levels <- c(
  "CD4+Tn", "MAIT", "CD8+Tem", "CD8+Temra", "CD8+Trm", "gdT", "Treg"
)
b_levels <- c("naiveB", "memoryB")
macrophage_levels <- c(
  "Monocytes", "IL1B+Mφ", "CXCL10+Mφ", "Kupffer", "TREM2+Mφ", "SPP1+Mφ"
)
fibroblast_levels <- c("Pericyte", "SMC", "HSC")

major_data <- composition_table(
  metadata,
  "celltype",
  major_levels,
  "Major cell type"
)
t_data <- composition_table(metadata, "sub_celltype", t_levels, "T cell")
b_data <- composition_table(metadata, "sub_celltype", b_levels, "B cell")
macrophage_data <- composition_table(
  metadata,
  "sub_celltype",
  macrophage_levels,
  "Macrophage"
)
fibroblast_data <- composition_table(
  metadata,
  "sub_celltype",
  fibroblast_levels,
  "Fibroblast"
)

composition_audit <- bind_rows(
  major_data %>% rename(state = celltype),
  t_data %>% rename(state = sub_celltype),
  b_data %>% rename(state = sub_celltype),
  macrophage_data %>% rename(state = sub_celltype),
  fibroblast_data %>% rename(state = sub_celltype)
) %>%
  select(panel, group, state, n, prop)
fwrite(
  as.data.table(composition_audit),
  file.path(output_dir, "Figure2_composition_values.tsv"),
  sep = "\t"
)

# Figure 2Q: accepted stromal/immune RoE heatmap.
stromal_celltypes <- c("Endothelial", "LSEC", "Fibroblast")
immune_celltypes <- c(
  "T", "NK", "B", "Plasma", "Macrophage", "cDC1", "cDC2", "pDC",
  "Mast", "Neutrophil"
)
metadata$lineage <- ifelse(
  metadata$celltype %in% stromal_celltypes,
  "Stromal",
  ifelse(metadata$celltype %in% immune_celltypes, "Immune", NA)
)

roe_matrix <- function(lineage, row_order) {
  subset_metadata <- metadata[
    !is.na(metadata$lineage) & metadata$lineage == lineage,
    ,
    drop = FALSE
  ]
  cells_per_sample <- table(subset_metadata$sample)
  subset_metadata <- subset_metadata[
    subset_metadata$sample %in%
      names(cells_per_sample)[cells_per_sample >= 20],
    ,
    drop = FALSE
  ]
  chi_square <- suppressWarnings(chisq.test(table(
    subset_metadata$sub_celltype,
    subset_metadata$group
  )))
  roe <- as.matrix(chi_square$observed / chi_square$expected)
  roe <- roe[, c("Healthy", "Adjacent", "Cirrhosis", "HCC")]
  rownames(roe) <- pretty_state(rownames(roe))
  roe[intersect(row_order, rownames(roe)), , drop = FALSE]
}

stromal_rows <- c("Endothelial", "LSEC", "Pericyte", "SMC", "HSC")
immune_rows <- c(
  "CD4+Tn", "MAIT", "CD8+Tem", "CD8+Temra", "CD8+Trm", "gdT", "Treg",
  "NK", "naiveB", "memoryB", "Plasma", "Kupffer", "IL1B+Mφ",
  "TREM2+Mφ", "Monocytes", "SPP1+Mφ", "CXCL10+Mφ", "cDC1", "cDC2",
  "pDC", "Mast", "Neutrophil"
)
stromal_roe <- roe_matrix("Stromal", stromal_rows)
immune_roe <- roe_matrix("Immune", immune_rows)
roe_all <- cbind(t(stromal_roe), t(immune_roe))

roe_audit <- as.data.table(as.table(roe_all))
setnames(roe_audit, c("group", "state", "roe"))
roe_audit[, lineage := ifelse(state %in% rownames(stromal_roe), "Stromal", "Immune")]
fwrite(
  roe_audit,
  file.path(output_dir, "Figure2Q_stromal_immune_roe.tsv"),
  sep = "\t"
)

message(
  "Computed composition and RoE from the 125,121-cell frozen atlas."
)
