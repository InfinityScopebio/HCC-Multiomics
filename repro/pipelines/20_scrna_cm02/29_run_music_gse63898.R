#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biobase)
  library(data.table)
  library(MuSiC)
  library(SingleCellExperiment)
})

arguments <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(
  sub("^--file=", "", arguments[grepl("^--file=", arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "config.R"))
config <- scrna_cm02_config(script_path)
raw_root <- file.path(config$project_root, "repo", "data", "public", "raw", "GSE63898")
series_path <- file.path(raw_root, "GSE63898_series_matrix.txt.gz")
platform_path <- file.path(raw_root, "GPL13667.soft.gz")
reference_path <- file.path(config$raw_music_root, "sc_eset.rds")
missing <- c(series_path, platform_path, reference_path)
if (any(!file.exists(missing))) {
  stop("Missing local GSE63898 inputs: ", paste(missing[!file.exists(missing)], collapse = ", "))
}

read_geo_matrix <- function(path) {
  connection <- gzfile(path, "rt")
  lines <- readLines(connection)
  close(connection)
  table_start <- grep("^!series_matrix_table_begin", lines)
  table_end <- grep("^!series_matrix_table_end", lines)
  matrix <- read.table(
    text = lines[(table_start + 1L):(table_end - 1L)],
    sep = "\t",
    header = TRUE,
    row.names = 1L,
    check.names = FALSE,
    quote = "\""
  )
  source_line <- grep("^!Sample_source_name_ch1", lines, value = TRUE)[1L]
  source <- gsub("\"", "", strsplit(source_line, "\t")[[1L]][-1L])
  list(matrix = as.matrix(matrix), source = source)
}

read_platform <- function(path) {
  connection <- gzfile(path, "rt")
  lines <- readLines(connection)
  close(connection)
  table_start <- grep("^!platform_table_begin", lines)
  table_end <- grep("^!platform_table_end", lines)
  read.table(
    text = lines[(table_start + 1L):(table_end - 1L)],
    sep = "\t",
    header = TRUE,
    quote = "",
    check.names = FALSE,
    comment.char = ""
  )
}

geo <- read_geo_matrix(series_path)
annotation <- read_platform(platform_path)
symbol_candidates <- c(
  "Gene Symbol", "Gene symbol", "Symbol", "GENE_SYMBOL",
  "GeneSymbol", "GeneName", "gene_assignment"
)
symbol_column <- symbol_candidates[symbol_candidates %in% names(annotation)][1L]
if (is.na(symbol_column)) {
  symbol_column <- grep("symbol|gene", names(annotation), value = TRUE, ignore.case = TRUE)[1L]
}
if (!"ID" %in% names(annotation) || is.na(symbol_column)) {
  stop("Could not identify GPL13667 probe or gene-symbol columns.")
}

annotation <- annotation[annotation$ID %in% rownames(geo$matrix), , drop = FALSE]
annotation <- annotation[
  !is.na(annotation[[symbol_column]]) & annotation[[symbol_column]] != "",
  ,
  drop = FALSE
]
annotation$gene <- sub(" ///.*", "", annotation[[symbol_column]])
probe_matrix <- geo$matrix[as.character(annotation$ID), , drop = FALSE]
ordering <- order(rowMeans(probe_matrix, na.rm = TRUE), decreasing = TRUE)
annotation <- annotation[ordering, , drop = FALSE]
probe_matrix <- probe_matrix[ordering, , drop = FALSE]
keep <- !duplicated(annotation$gene)
gene_matrix <- probe_matrix[keep, , drop = FALSE]
rownames(gene_matrix) <- annotation$gene[keep]

annotated_path <- file.path(config$raw_music_root, "GSE63898_sym.rds")
saveRDS(
  list(mat = gene_matrix, group = geo$source, sample = colnames(gene_matrix)),
  annotated_path
)

linear_matrix <- 2^gene_matrix - 1
linear_matrix[linear_matrix < 0] <- 0
reference <- readRDS(reference_path)
single_cell <- SingleCellExperiment(
  assays = list(counts = exprs(reference)),
  colData = pData(reference)
)
common_genes <- intersect(rownames(linear_matrix), rownames(reference))
if (length(common_genes) < 500L) {
  stop("GSE63898 has too few genes shared with the MuSiC reference.")
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42L)
fit <- music_prop(
  bulk.mtx = linear_matrix,
  sc.sce = single_cell,
  clusters = "cellType",
  samples = "sampleID",
  select.ct = unique(single_cell$cellType),
  verbose = FALSE
)
proportions <- as.data.table(fit$Est.prop.weighted, keep.rownames = "sampleID")
proportions[, cohort := "GSE63898"]
output_path <- file.path(config$raw_music_root, "prop_GSE63898.csv")
fwrite(proportions, output_path)
fwrite(
  data.table(
    cohort = "GSE63898",
    n_samples = nrow(proportions),
    n_states = ncol(proportions) - 2L,
    n_symbol_rows = nrow(gene_matrix),
    n_shared_reference_genes = length(common_genes),
    seed = 42L,
    output_sha256 = unname(tools::sha256sum(output_path))
  ),
  file.path(config$raw_music_root, "prop_GSE63898_provenance.tsv"),
  sep = "\t"
)
message("GSE63898 MuSiC complete: ", nrow(proportions), " samples.")
