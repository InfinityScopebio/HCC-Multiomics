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
cohort_arguments <- commandArgs(trailingOnly = TRUE)
if (length(cohort_arguments) != 1L ||
    !cohort_arguments[[1L]] %in% c("local", "TCGA", "GTEx")) {
  stop("Usage: 28_run_music_rnaseq.R <local|TCGA|GTEx>")
}
cohort <- cohort_arguments[[1L]]

reference_path <- file.path(config$raw_music_root, "sc_eset.rds")
if (!file.exists(reference_path)) {
  stop("Run 27_build_music_reference.R first.")
}
reference <- readRDS(reference_path)
single_cell <- SingleCellExperiment(
  assays = list(counts = exprs(reference)),
  colData = pData(reference)
)

prepare_matrix <- function(matrix) {
  if (anyDuplicated(rownames(matrix))) {
    ordering <- order(rowMeans(matrix, na.rm = TRUE), decreasing = TRUE)
    matrix <- matrix[ordering, , drop = FALSE]
    matrix <- matrix[!duplicated(rownames(matrix)), , drop = FALSE]
  }
  storage.mode(matrix) <- "numeric"
  matrix
}

if (cohort == "local") {
  tissue_path <- file.path(
    config$project_root,
    "repo",
    "results",
    "00_preprocess",
    "multiomics",
    "dat",
    "tissue_rna.rds"
  )
  tissue <- readRDS(tissue_path)
  bulk_matrix <- tissue$mat
  rownames(bulk_matrix) <- tissue$feature_id
} else if (cohort == "TCGA") {
  expression_path <- file.path(
    config$project_root,
    "repo",
    "data",
    "public",
    "raw",
    "TCGA-LIHC",
    "LIHC_HiSeqV2.gz"
  )
  expression <- fread(
    cmd = paste("gunzip -c", shQuote(expression_path)),
    header = TRUE
  )
  genes <- expression[[1L]]
  bulk_matrix <- as.matrix(expression[, -1L])
  rownames(bulk_matrix) <- genes
  bulk_matrix <- 2^bulk_matrix - 1
  bulk_matrix[bulk_matrix < 0] <- 0
} else {
  gtex_path <- file.path(
    config$project_root,
    "repo",
    "data",
    "public",
    "raw",
    "GTEx-Liver",
    "GTEx_LIVER_recount3.rds"
  )
  gtex <- readRDS(gtex_path)
  bulk_matrix <- gtex$counts
  row_data <- as.data.frame(gtex$rowData)
  if ("gene_name" %in% names(row_data)) {
    rownames(bulk_matrix) <- row_data$gene_name
    bulk_matrix <- bulk_matrix[
      !is.na(rownames(bulk_matrix)) & rownames(bulk_matrix) != "",
      ,
      drop = FALSE
    ]
  }
}

bulk_matrix <- prepare_matrix(bulk_matrix)
common_genes <- intersect(rownames(reference), rownames(bulk_matrix))
if (length(common_genes) < 500L) {
  stop(cohort, " has too few genes shared with the MuSiC reference: ", length(common_genes))
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42L)
fit <- music_prop(
  bulk.mtx = bulk_matrix,
  sc.sce = single_cell,
  clusters = "cellType",
  samples = "sampleID",
  select.ct = unique(single_cell$cellType),
  verbose = FALSE
)
proportions <- as.data.table(fit$Est.prop.weighted, keep.rownames = "sampleID")
proportions[, cohort := cohort]
output_path <- file.path(config$raw_music_root, paste0("prop_", cohort, ".csv"))
fwrite(proportions, output_path)

fwrite(
  data.table(
    cohort = cohort,
    n_samples = nrow(proportions),
    n_states = ncol(proportions) - 2L,
    n_bulk_genes = nrow(bulk_matrix),
    n_shared_reference_genes = length(common_genes),
    seed = 42L,
    output_sha256 = unname(tools::sha256sum(output_path))
  ),
  file.path(config$raw_music_root, paste0("prop_", cohort, "_provenance.tsv")),
  sep = "\t"
)
message(cohort, " MuSiC complete: ", nrow(proportions), " samples.")
