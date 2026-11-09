#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
})

script_path <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

protein_path <- data_path("tissue_protein.rds")
network_path <- reference_path("PROGENy", "progeny_net.rds")
assert_files(c(protein_path, network_path))

protein <- readRDS(protein_path)
protein_matrix <- build_gene_matrix(
  protein$mat, protein$meta$Gene, take = "max"
)
network <- readRDS(network_path)
stopifnot(
  all(c("source", "target", "weight") %in% colnames(network)),
  nrow(network) == 7000L
)

activity_long <- decoupleR::run_mlm(
  mat = protein_matrix,
  net = network,
  .source = "source",
  .target = "target",
  .mor = "weight",
  minsize = 5
)
activity_wide <- dcast(
  as.data.table(activity_long),
  source ~ condition,
  value.var = "score"
)
pathways <- activity_wide$source
activity_matrix <- as.matrix(activity_wide[, -"source"])
rownames(activity_matrix) <- pathways

activity_output <- as.data.table(activity_matrix, keep.rownames = "pathway")
fwrite(
  activity_output,
  result_path("progeny", "PROGENy_Protein_per_sample.csv")
)

metadata <- paired_tissue_metadata()
paired <- paired_activity_test(activity_matrix, metadata$tissue)
setnames(paired, "item", "pathway")
fwrite(
  paired,
  result_path("progeny", "PROGENy_paired_HCC_vs_Adjacent.csv")
)

plot_input <- paired[, .(
  pathway,
  tissue_Prot_delta = delta,
  tissue_Prot_padj = p_adj
)]
fwrite(
  plot_input,
  result_path("progeny", "PROGENy_tissue_plot_input.csv")
)

message(
  "PROGENy complete: ", nrow(activity_matrix), " pathways x ",
  ncol(activity_matrix), " samples."
)
