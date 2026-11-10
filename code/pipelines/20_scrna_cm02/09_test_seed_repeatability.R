#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clue)
  library(data.table)
  library(digest)
  library(foreach)
  library(Matrix)
  library(NMF)
  library(Seurat)
})

arguments <- commandArgs(trailingOnly = FALSE)
this_script <- normalizePath(
  sub("^--file=", "", arguments[grepl("^--file=", arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(this_script), "config.R"))
config <- ensure_output_dirs(scrna_cm02_config(this_script))
setwd(config$project_root)

object_specifications <- list(
  atlas = list(
    path = config$atlas_paths$uncurated,
    stage = "atlas",
    reduction = "integrated.rpca",
    dimensions = 1:30
  ),
  epithelial = list(
    path = file.path(
      config$result_root,
      "epithelial_uncurated_raw_derived.rds"
    ),
    stage = "epithelial",
    reduction = "integrated.rpca",
    dimensions = 1:15
  ),
  t_cell = list(
    path = file.path(
      config$result_root,
      "t_cell_uncurated_raw_derived.rds"
    ),
    stage = "t_cell",
    reduction = "integrated.rpca",
    dimensions = 1:20
  ),
  macrophage = list(
    path = file.path(
      config$result_root,
      "macrophage_uncurated_raw_derived.rds"
    ),
    stage = "macrophage",
    reduction = "harmony",
    dimensions = 1:20
  ),
  fibroblast = list(
    path = file.path(
      config$result_root,
      "fibroblast_uncurated_raw_derived.rds"
    ),
    stage = "fibroblast",
    reduction = "harmony",
    dimensions = 1:15
  ),
  b_cell = list(
    path = file.path(
      config$result_root,
      "b_cell_uncurated_raw_derived.rds"
    ),
    stage = "b_cell",
    reduction = "harmony",
    dimensions = 1:15
  ),
  hepatocyte_initial = list(
    path = file.path(
      config$result_root,
      "hepatocyte_initial_uncurated_raw_derived.rds"
    ),
    stage = "hepatocyte_initial",
    reduction = "integrated.rpca",
    dimensions = 1:50
  ),
  hepatocyte_clean = list(
    path = file.path(
      config$result_root,
      "hepatocyte_clean_uncurated_raw_derived.rds"
    ),
    stage = "hepatocyte_clean",
    reduction = "integrated.rpca",
    dimensions = 1:50
  )
)

sparse_checksum <- function(value) {
  digest::digest(
    list(
      dimensions = dim(value),
      dimnames = dimnames(value),
      i = value@i,
      p = value@p,
      x = value@x
    ),
    algo = "sha256",
    serialize = TRUE
  )
}

numeric_checksum <- function(value) {
  digest::digest(
    list(
      dimensions = dim(value),
      dimnames = dimnames(value),
      values = as.numeric(value)
    ),
    algo = "sha256",
    serialize = TRUE
  )
}

edge_jaccard <- function(left, right) {
  left_triplet <- summary(left)
  right_triplet <- summary(right)
  left_key <- paste(left_triplet$i, left_triplet$j, sep = ":")
  right_key <- paste(right_triplet$i, right_triplet$j, sep = ":")
  length(intersect(left_key, right_key)) / length(union(left_key, right_key))
}

adjusted_rand_index <- function(left, right) {
  contingency <- table(left, right)
  choose2 <- function(values) values * (values - 1) / 2
  cell_pairs <- sum(choose2(contingency))
  row_pairs <- sum(choose2(rowSums(contingency)))
  column_pairs <- sum(choose2(colSums(contingency)))
  total_pairs <- choose2(sum(contingency))
  expected <- row_pairs * column_pairs / total_pairs
  denominator <- (row_pairs + column_pairs) / 2 - expected
  if (denominator == 0) {
    return(ifelse(identical(left, right), 1, 0))
  }
  (cell_pairs - expected) / denominator
}

procrustes_correlation <- function(left, right) {
  left_scaled <- scale(left)
  right_scaled <- scale(right)
  singular <- svd(crossprod(left_scaled, right_scaled))
  aligned <- left_scaled %*% singular$u %*% t(singular$v)
  cor(as.numeric(aligned), as.numeric(right_scaled))
}

run_neighbor_cluster_umap <- function(object, specification) {
  stage <- specification$stage
  set_recorded_seed(config, stage, "FindNeighbors")
  object <- FindNeighbors(
    object,
    reduction = specification$reduction,
    dims = specification$dimensions,
    nn.method = "annoy",
    n.trees = 50L,
    annoy.metric = "euclidean",
    graph.name = c("repeat_nn", "repeat_snn"),
    verbose = FALSE
  )
  set_recorded_seed(config, stage, "FindClusters")
  object <- FindClusters(
    object,
    graph.name = "repeat_snn",
    cluster.name = "repeat_cluster",
    resolution = if (stage %in% c("atlas", "t_cell", "hepatocyte_initial",
      "hepatocyte_clean")) {
      0.5
    } else if (stage == "macrophage") {
      0.4
    } else {
      0.3
    },
    algorithm = 1L,
    n.start = 10L,
    n.iter = 10L,
    random.seed = 0L,
    verbose = FALSE
  )
  set_recorded_seed(config, stage, "RunUMAP")
  object <- RunUMAP(
    object,
    reduction = specification$reduction,
    dims = specification$dimensions,
    umap.method = "uwot",
    n.neighbors = if (stage == "atlas") 30L else 15L,
    min.dist = if (stage == "atlas") 0.3 else 0.05,
    spread = if (stage == "atlas") 1 else 1.5,
    seed.use = 42L,
    reduction.name = "repeat_umap",
    reduction.key = "repeatUMAP_",
    verbose = FALSE
  )
  list(
    graph = as(object[["repeat_snn"]], "dgCMatrix"),
    cluster = as.character(object$repeat_cluster),
    umap = Embeddings(object, "repeat_umap")
  )
}

rows <- list()
for (object_name in names(object_specifications)) {
  specification <- object_specifications[[object_name]]
  if (!file.exists(specification$path)) {
    stop("Missing repeatability input: ", specification$path)
  }
  message("Repeatability test: ", object_name)
  first_object <- readRDS(specification$path)
  assert_no_cc(
    first_object$sample,
    paste0(object_name, " replicate 1 input")
  )
  first_input_checksum <- numeric_checksum(
    Embeddings(first_object, specification$reduction)
  )
  first <- run_neighbor_cluster_umap(first_object, specification)
  rm(first_object)
  invisible(gc())

  second_object <- readRDS(specification$path)
  assert_no_cc(
    second_object$sample,
    paste0(object_name, " replicate 2 input")
  )
  second_input_checksum <- numeric_checksum(
    Embeddings(second_object, specification$reduction)
  )
  if (!identical(first_input_checksum, second_input_checksum)) {
    stop(object_name, ": saved reduction changed between replicates.")
  }
  second <- run_neighbor_cluster_umap(second_object, specification)
  rm(second_object)
  invisible(gc())

  graph_difference <- drop0(
    as(first$graph - second$graph, "dgCMatrix")
  )
  graph_max_difference <- if (length(graph_difference@x) == 0L) {
    0
  } else {
    max(abs(graph_difference@x))
  }
  umap_max_difference <- max(abs(first$umap - second$umap))
  cluster_agreement <- mean(first$cluster == second$cluster)
  rows[[object_name]] <- rbindlist(
    list(
      data.table(
        stage = specification$stage,
        algorithm = "FindNeighbors/Annoy SNN",
        seed = config$master_seed,
        replicate1_checksum = sparse_checksum(first$graph),
        replicate2_checksum = sparse_checksum(second$graph),
        max_abs_diff = graph_max_difference,
        agreement = NA_real_,
        exact_match = graph_max_difference == 0 &&
          sparse_checksum(first$graph) == sparse_checksum(second$graph),
        semantic_metric = "nonzero_edge_jaccard",
        semantic_value = edge_jaccard(first$graph, second$graph),
        status = if (
          graph_max_difference == 0 &&
            sparse_checksum(first$graph) == sparse_checksum(second$graph)
        ) {
          "fixed"
        } else {
          "blocked"
        },
        notes = paste0(
          "Independent disk reloads; input reduction checksum=",
          first_input_checksum
        )
      ),
      data.table(
        stage = specification$stage,
        algorithm = "FindClusters/Louvain",
        seed = 0L,
        replicate1_checksum = digest::digest(
          first$cluster,
          algo = "sha256"
        ),
        replicate2_checksum = digest::digest(
          second$cluster,
          algo = "sha256"
        ),
        max_abs_diff = NA_real_,
        agreement = cluster_agreement,
        exact_match = identical(first$cluster, second$cluster),
        semantic_metric = "adjusted_rand_index",
        semantic_value = adjusted_rand_index(
          first$cluster,
          second$cluster
        ),
        status = if (cluster_agreement == 1) "fixed" else "blocked",
        notes = paste0(
          "Independent disk reloads and independently rebuilt SNN; ",
          "input reduction checksum=",
          first_input_checksum
        )
      ),
      data.table(
        stage = specification$stage,
        algorithm = "RunUMAP/uwot",
        seed = 42L,
        replicate1_checksum = numeric_checksum(first$umap),
        replicate2_checksum = numeric_checksum(second$umap),
        max_abs_diff = umap_max_difference,
        agreement = NA_real_,
        exact_match = umap_max_difference == 0 &&
          numeric_checksum(first$umap) == numeric_checksum(second$umap),
        semantic_metric = "procrustes_correlation",
        semantic_value = procrustes_correlation(
          first$umap,
          second$umap
        ),
        status = if (
          umap_max_difference == 0 &&
            numeric_checksum(first$umap) == numeric_checksum(second$umap)
        ) {
          "fixed"
        } else {
          "blocked"
        },
        notes = paste0(
          "Independent disk reloads; input reduction checksum=",
          first_input_checksum
        )
      )
    ),
    fill = TRUE
  )
  rm(first, second, graph_difference)
  invisible(gc())
}

frequency_path <- config$covarnet_paths$frequency
foreach::registerDoSEQ()
run_nmf <- function() {
  frequency <- readRDS(frequency_path)$mat_norm_filtered
  input_checksum <- numeric_checksum(frequency)
  set_recorded_seed(config, "covarnet", "nsNMF_K5")
  result <- NMF::nmf(
    frequency,
    rank = 5L,
    method = "nsNMF",
    seed = rep(77L, 6L),
    nrun = 30L,
    .options = "v",
    .pbackend = "seq"
  )
  list(result = result, input_checksum = input_checksum)
}
first_nmf_run <- run_nmf()
second_nmf_run <- run_nmf()
if (!identical(
  first_nmf_run$input_checksum,
  second_nmf_run$input_checksum
)) {
  stop("Saved NMF input changed between independent replicates.")
}
first_nmf <- first_nmf_run$result
second_nmf <- second_nmf_run$result
first_basis <- NMF::basis(first_nmf)
second_basis <- NMF::basis(second_nmf)
first_coefficient <- NMF::coef(first_nmf)
second_coefficient <- NMF::coef(second_nmf)
nmf_max_difference <- max(
  abs(first_basis - second_basis),
  abs(first_coefficient - second_coefficient)
)
first_basis_norm <- sweep(
  first_basis,
  2L,
  sqrt(colSums(first_basis^2)),
  "/"
)
second_basis_norm <- sweep(
  second_basis,
  2L,
  sqrt(colSums(second_basis^2)),
  "/"
)
nmf_similarity <- crossprod(first_basis_norm, second_basis_norm)
nmf_assignment <- clue::solve_LSAP(nmf_similarity, maximum = TRUE)
nmf_matched_cosine <- mean(
  nmf_similarity[cbind(seq_len(nrow(nmf_similarity)), nmf_assignment)]
)
rows[["covarnet"]] <- data.table(
  stage = "covarnet",
  algorithm = "NMF/nsNMF K=5",
  seed = 77L,
  replicate1_checksum = numeric_checksum(
    matrix(c(first_basis, first_coefficient), ncol = 1L)
  ),
  replicate2_checksum = numeric_checksum(
    matrix(c(second_basis, second_coefficient), ncol = 1L)
  ),
  max_abs_diff = nmf_max_difference,
  agreement = NA_real_,
  exact_match = nmf_max_difference == 0,
  semantic_metric = "mean_matched_basis_cosine",
  semantic_value = nmf_matched_cosine,
  status = if (nmf_max_difference == 0) "fixed" else "blocked",
  notes = paste0(
    "Independent disk reloads; sequential foreach backend; input checksum=",
    first_nmf_run$input_checksum
  )
)

repeatability <- rbindlist(rows, fill = TRUE)
fwrite(
  repeatability,
  file.path(config$audit_root, "scrna_seed_repeatability.tsv"),
  sep = "\t"
)

publish_strict_seed_audit(config, repeatability)

if (any(repeatability$status == "blocked")) {
  warning(
    "One or more random algorithms drifted across identical duplicate runs; ",
    "see scrna_seed_repeatability.tsv."
  )
} else {
  message("PASS: all tested Annoy/Louvain/UMAP/NMF duplicate runs are exact.")
}
