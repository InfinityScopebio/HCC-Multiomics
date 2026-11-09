#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
  library(Seurat)
})

arguments <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(
  sub("^--file=", "", arguments[grepl("^--file=", arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "config.R"))
config <- scrna_cm02_config(script_path)
output_dir <- file.path(config$raw_to_figure3_root, "hepc_activity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hepatocyte_path <- config$curated_atlas_paths$hepatocyte
progeny_path <- file.path(
  config$project_root,
  "repo",
  "data",
  "reference",
  "PROGENy",
  "progeny_net.rds"
)
collectri_path <- file.path(
  config$project_root,
  "repo",
  "data",
  "reference",
  "CollecTRI",
  "collectri_net.csv"
)
required <- c(hepatocyte_path, progeny_path, collectri_path)
if (any(!file.exists(required))) {
  stop("Missing Hep_C activity inputs: ", paste(required[!file.exists(required)], collapse = ", "))
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42L)
hepatocytes <- readRDS(hepatocyte_path)
DefaultAssay(hepatocytes) <- "RNA"
hepatocytes <- NormalizeData(hepatocytes, verbose = FALSE)
expression <- GetAssayData(hepatocytes, assay = "RNA", layer = "data")

progeny_network <- as.data.frame(readRDS(progeny_path))
if (!"mor" %in% names(progeny_network)) {
  progeny_network$mor <- progeny_network$weight
}
progeny_genes <- intersect(rownames(expression), unique(progeny_network$target))
progeny_activity <- run_mlm(
  mat = as.matrix(expression[progeny_genes, , drop = FALSE]),
  network = progeny_network,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5L
)
progeny_activity <- as.data.table(progeny_activity)
progeny_activity[, sub := hepatocytes@meta.data[condition, "sub_celltype"]]
progeny_activity <- progeny_activity[!is.na(sub) & grepl("^Hep_C", sub)]
progeny_summary <- progeny_activity[, .(
  mean_score = mean(score, na.rm = TRUE)
), by = .(source, sub)]
progeny_matrix <- as.matrix(dcast(
  progeny_summary,
  source ~ sub,
  value.var = "mean_score"
), rownames = "source")
fwrite(
  as.data.table(progeny_matrix, keep.rownames = "pathway"),
  file.path(output_dir, "progeny_HepC.csv")
)

collectri_network <- fread(collectri_path)
# Preserve the article computation exactly: the historical script reused the
# PROGENy-restricted expression matrix for the subsequent CollecTRI ULM.
collectri_genes <- intersect(progeny_genes, unique(collectri_network$target))
collectri_activity <- run_ulm(
  mat = as.matrix(expression[collectri_genes, , drop = FALSE]),
  network = as.data.frame(collectri_network),
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5L
)
collectri_activity <- as.data.table(collectri_activity)
collectri_activity[, sub := hepatocytes@meta.data[condition, "sub_celltype"]]
collectri_activity <- collectri_activity[!is.na(sub) & grepl("^Hep_C", sub)]
collectri_summary <- collectri_activity[, .(
  mean_score = mean(score, na.rm = TRUE)
), by = .(source, sub)]
collectri_matrix <- as.matrix(dcast(
  collectri_summary,
  source ~ sub,
  value.var = "mean_score"
), rownames = "source")
fwrite(
  as.data.table(collectri_matrix, keep.rownames = "TF"),
  file.path(output_dir, "dorothea_HepC.csv")
)

fwrite(
  data.table(
    input = c("hepatocyte_atlas", "PROGENy_network", "CollecTRI_network"),
    path = normalizePath(required, mustWork = TRUE),
    sha256 = unname(tools::sha256sum(required)),
    seed = 42L
  ),
  file.path(output_dir, "activity_input_provenance.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    method = c("PROGENy", "CollecTRI"),
    expression_gene_scope = c(
      "PROGENy network targets present in Hep_C expression",
      "CollecTRI targets intersected with the historical PROGENy-restricted gene set"
    ),
    n_input_genes = c(length(progeny_genes), length(collectri_genes))
  ),
  file.path(output_dir, "activity_method_contract.tsv"),
  sep = "\t"
)
message(
  "Hep_C activities rebuilt: ", nrow(progeny_matrix), " pathways and ",
  nrow(collectri_matrix), " TFs."
)
