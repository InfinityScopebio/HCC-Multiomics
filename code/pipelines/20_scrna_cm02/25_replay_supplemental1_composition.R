#!/usr/bin/env Rscript

# Compute per-sample major-cell-type composition from the frozen curated atlas.

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

output_dir <- file.path(config$result_root, "supplemental_figure1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(config$curated_atlas_paths$final)
metadata <- object[[]]
rm(object)
invisible(gc(verbose = FALSE))

disease_levels <- c("Healthy", "Adjacent", "Cirrhosis", "HCC")
source_levels <- c("GSE115469", "GSE136103", "GSE149614", "LiuF")
major_levels <- c(
  "T", "NK", "B", "Plasma", "Macrophage", "cDC1", "cDC2", "pDC",
  "Neutrophil", "Mast", "Hepatocyte", "Cholangiocyte", "Endothelial",
  "LSEC", "Fibroblast"
)

required_columns <- c("sample", "group", "source", "celltype")
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0L) {
  stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
}
if (any(as.character(metadata$sample) %in% config$prohibited_samples)) {
  stop("CC1-CC3 entered the frozen Supplemental Figure 1 atlas.")
}

for (field_spec in list(
  list(field = "group", allowed = disease_levels),
  list(field = "source", allowed = source_levels),
  list(field = "celltype", allowed = major_levels)
)) {
  values <- as.character(metadata[[field_spec$field]])
  unknown <- setdiff(unique(values), field_spec$allowed)
  unknown <- unknown[!is.na(unknown)]
  if (anyNA(values) || length(unknown) > 0L) {
    stop(
      "Missing or unmapped ",
      field_spec$field,
      " values: ",
      paste(unknown, collapse = ", ")
    )
  }
}

composition_raw <- metadata %>%
  count(sample, group, source, celltype, name = "n_cells") %>%
  group_by(sample) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

t_fraction <- composition_raw %>%
  filter(celltype == "T") %>%
  select(sample, t_fraction = fraction)

sample_key <- metadata %>%
  distinct(sample, group, source) %>%
  left_join(t_fraction, by = "sample") %>%
  mutate(
    t_fraction = coalesce(t_fraction, 0),
    group_order = match(group, disease_levels)
  ) %>%
  arrange(group_order, desc(t_fraction), sample)

if (anyDuplicated(sample_key$sample)) {
  stop("Sample IDs map to multiple group/source combinations.")
}
sample_order <- sample_key$sample
composition <- composition_raw %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    group = factor(group, levels = disease_levels),
    source = factor(source, levels = source_levels),
    celltype = factor(celltype, levels = major_levels)
  )

fraction_check <- composition %>%
  group_by(sample) %>%
  summarise(total = sum(fraction), .groups = "drop")
if (any(abs(fraction_check$total - 1) > 1e-10)) {
  stop("Per-sample composition fractions do not sum to one.")
}

fwrite(
  as.data.table(composition),
  file.path(output_dir, "S1K_major_celltype_composition_by_sample.tsv"),
  sep = "\t"
)
fwrite(
  as.data.table(sample_key),
  file.path(output_dir, "S1K_sample_order.tsv"),
  sep = "\t"
)

message(
  "Computed per-sample composition from ",
  format(nrow(metadata), big.mark = ","),
  " frozen cells across ",
  length(sample_order),
  " samples."
)
