#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(uwot)
})

arguments <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(
  sub("^--file=", "", arguments[grepl("^--file=", arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(script_path), "config.R"))
config <- scrna_cm02_config(script_path)
dir.create(config$raw_trajectory_root, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("local", "TCGA", "GTEx", "GSE63898")
proportion_paths <- file.path(
  config$raw_music_root,
  paste0("prop_", cohorts, ".csv")
)
if (any(!file.exists(proportion_paths))) {
  stop(
    "Missing MuSiC proportions: ",
    paste(proportion_paths[!file.exists(proportion_paths)], collapse = ", ")
  )
}
combined <- rbindlist(lapply(proportion_paths, fread), fill = TRUE)
combined[is.na(combined)] <- 0
combined[, group := NA_character_]

sample_sheet_path <- file.path(
  config$project_root,
  "repo",
  "results",
  "00_preprocess",
  "multiomics",
  "dat",
  "sample_sheet.csv"
)
sample_sheet <- fread(sample_sheet_path)
local_group <- setNames(
  sample_sheet[source == "tissue", group_4],
  sample_sheet[source == "tissue", sample_id]
)
clinical <- fread(config$raw_downstream_paths$tcga_clinical)
tcga_type <- setNames(clinical$sample_type, clinical$sampleID)
gse63898 <- readRDS(file.path(config$raw_music_root, "GSE63898_sym.rds"))
gse63898_group <- setNames(gse63898$group, gse63898$sample)

combined[cohort == "GTEx", group := "Healthy"]
combined[cohort == "GSE63898", group := fcase(
  gse63898_group[sampleID] == "cirrhosis", "Cirrhosis",
  gse63898_group[sampleID] == "hepatocellular carcinoma", "HCC",
  default = NA_character_
)]
combined[cohort == "local", group := fcase(
  local_group[sampleID] == "HCC_tumor", "HCC",
  local_group[sampleID] == "Adjacent", "Adjacent",
  default = NA_character_
)]
combined[cohort == "TCGA", group := fcase(
  tcga_type[sampleID] == "Primary Tumor", "HCC",
  tcga_type[sampleID] == "Solid Tissue Normal", "Adjacent",
  tcga_type[sampleID] == "Recurrent Tumor", "HCC",
  default = NA_character_
)]

fwrite(
  combined,
  file.path(config$raw_trajectory_root, "combined_freq.csv")
)

trajectory_data <- combined[!is.na(group)]
if (nrow(trajectory_data) != 1087L) {
  stop("Expected 1,087 grouped bulk samples, found ", nrow(trajectory_data))
}
cm02_members <- config$cm02_anchor_states
missing_members <- setdiff(cm02_members, names(trajectory_data))
if (length(missing_members) > 0L) {
  stop("MuSiC output lacks CM02 members: ", paste(missing_members, collapse = ", "))
}

frequency_matrix <- as.matrix(trajectory_data[, ..cm02_members])
rownames(frequency_matrix) <- trajectory_data$sampleID
scaled_matrix <- frequency_matrix
for (cohort in unique(trajectory_data$cohort)) {
  cohort_rows <- which(trajectory_data$cohort == cohort)
  if (length(cohort_rows) >= 5L) {
    scaled_matrix[cohort_rows, ] <- scale(
      frequency_matrix[cohort_rows, , drop = FALSE]
    )
  }
}
scaled_matrix[is.na(scaled_matrix)] <- 0
scaled_matrix <- pmin(pmax(scaled_matrix, -3), 3)

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42L)
umap <- uwot::umap(
  scaled_matrix,
  n_neighbors = 50L,
  n_components = 2L,
  metric = "euclidean",
  min_dist = 0.05,
  spread = 1.5,
  n_epochs = 500L,
  ret_nn = TRUE
)
embedding_matrix <- umap$embedding
rownames(embedding_matrix) <- rownames(scaled_matrix)
neighbors <- umap$nn$euclidean$idx
edges <- do.call(rbind, lapply(seq_len(nrow(neighbors)), function(index) {
  cbind(index, neighbors[index, -1L])
}))
graph <- simplify(graph_from_edgelist(edges, directed = FALSE))
cluster <- as.integer(membership(cluster_louvain(graph)))

embedding <- data.table(
  sampleID = rownames(embedding_matrix),
  phate1 = embedding_matrix[, 1L],
  phate2 = embedding_matrix[, 2L],
  Cluster = as.character(cluster)
)
metadata <- trajectory_data[, .(
  sampleID,
  Group = group,
  Cohort = cohort
)]
embedding <- merge(embedding, metadata, by = "sampleID")

healthy_rows <- which(trajectory_data$group == "Healthy")
healthy_centroid <- colMeans(scaled_matrix[healthy_rows, , drop = FALSE])
pseudotime <- sqrt(rowSums(sweep(scaled_matrix, 2L, healthy_centroid, "-")^2))
pseudotime <- (pseudotime - min(pseudotime)) /
  (max(pseudotime) - min(pseudotime))
embedding[, Pseudotime := pseudotime[sampleID]]

fwrite(
  embedding[, .(sampleID, phate1, phate2, Cluster, Group, Cohort)],
  file.path(config$raw_trajectory_root, "umap_CM02_embedding.csv")
)
fwrite(
  embedding,
  file.path(config$raw_trajectory_root, "umap_CM02_pseudotime.csv")
)
saveRDS(
  list(
    scaled_cm02_frequencies = scaled_matrix,
    healthy_centroid = healthy_centroid,
    umap_neighbors = neighbors,
    seed = 42L
  ),
  file.path(config$raw_trajectory_root, "trajectory_intermediates.rds")
)


fwrite(
  trajectory_data[, .N, by = .(cohort, group)][order(cohort, group)],
  file.path(config$raw_trajectory_root, "trajectory_cohort_group_counts.tsv"),
  sep = "\t"
)
message("Bulk CM02 trajectory rebuilt for ", nrow(embedding), " samples.")
