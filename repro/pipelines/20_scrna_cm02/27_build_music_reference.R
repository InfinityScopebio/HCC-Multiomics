#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biobase)
  library(data.table)
  library(Seurat)
})

arguments <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(
  sub("^--file=", "", arguments[grepl("^--file=", arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "config.R"))
config <- scrna_cm02_config(script_path)
dir.create(config$raw_music_root, recursive = TRUE, showWarnings = FALSE)

atlas_path <- config$curated_atlas_paths$final
if (!file.exists(atlas_path)) {
  stop("Frozen curated atlas is missing: ", atlas_path)
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42L)
atlas <- readRDS(atlas_path)
metadata <- atlas@meta.data
metadata$sub_celltype <- as.character(metadata$sub_celltype)
metadata$sample <- as.character(metadata$sample)
keep <- !is.na(metadata$sub_celltype) &
  metadata$sub_celltype != "Cholangiocyte"
atlas <- subset(atlas, cells = colnames(atlas)[keep])

metadata <- as.data.table(atlas@meta.data, keep.rownames = "cell")
metadata[, within_sample_state_index := seq_len(.N), by = .(sample, sub_celltype)]
keep_cells <- metadata[within_sample_state_index <= 30L, cell]
reference <- subset(atlas, cells = keep_cells)
reference <- FindVariableFeatures(reference, nfeatures = 2000L, verbose = FALSE)
hvg <- VariableFeatures(reference)
hvg <- hvg[!grepl("^(MT-|RPL|RPS|HB[ABM])", hvg)]

counts <- as.matrix(GetAssayData(
  reference,
  assay = "RNA",
  layer = "counts"
)[hvg, , drop = FALSE])
reference_metadata <- reference@meta.data
phenotype <- AnnotatedDataFrame(data.frame(
  sampleID = as.character(reference_metadata$sample),
  cellType = as.character(reference_metadata$sub_celltype),
  row.names = colnames(counts),
  stringsAsFactors = FALSE
))
expression_set <- ExpressionSet(
  assayData = counts,
  phenoData = phenotype
)
output_path <- file.path(config$raw_music_root, "sc_eset.rds")
saveRDS(expression_set, output_path)

provenance <- data.table(
  input_path = normalizePath(atlas_path, mustWork = TRUE),
  input_sha256 = unname(tools::sha256sum(atlas_path)),
  output_path = normalizePath(output_path, mustWork = TRUE),
  output_sha256 = unname(tools::sha256sum(output_path)),
  seed = 42L,
  n_genes = nrow(counts),
  n_cells = ncol(counts),
  n_cell_states = uniqueN(reference_metadata$sub_celltype),
  n_samples = uniqueN(reference_metadata$sample),
  selection = "first 30 cells per sample and sub_celltype in frozen atlas order"
)
fwrite(
  provenance,
  file.path(config$raw_music_root, "sc_reference_provenance.tsv"),
  sep = "\t"
)
message(
  "MuSiC reference rebuilt: ", nrow(counts), " genes x ",
  ncol(counts), " cells."
)
