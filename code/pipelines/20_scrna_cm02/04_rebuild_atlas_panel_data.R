#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clue)
  library(data.table)
  library(future)
  library(harmony)
  library(Matrix)
  library(RANN)
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
source(file.path(config$project_root, "code", "lib", "R", "runtime_flags.R"))
options(future.globals.maxSize = 60 * 1024^3)
future::plan(future::sequential)
assert_production_input(config, config$atlas_paths$integrated)

panel_dir <- file.path(config$result_root, "atlas_panel_data")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
atlas <- readRDS(config$atlas_paths$integrated)
assert_no_cc(atlas$sample, "major atlas before subclustering")

run_subcluster <- function(
    object,
    stage,
    method,
    dimensions,
    resolution,
    minimum_cells_per_sample,
    k_weight = 30L,
    pca_components = max(dimensions) + 5L) {
  checkpoint_path <- file.path(
    config$result_root,
    paste0(stage, "_uncurated_raw_derived.rds")
  )
  if (
    resume_enabled("LIUF_RESUME_SUBCLUSTERS") &&
      file.exists(checkpoint_path)
  ) {
    message("Resuming ", stage, " from pre-curation checkpoint.")
    checkpoint <- readRDS(checkpoint_path)
    assert_no_cc(
      checkpoint$sample,
      paste0(stage, " pre-curation checkpoint")
    )
    return(checkpoint)
  }
  sample_counts <- table(object$sample)
  retained_samples <- names(sample_counts)[sample_counts >= minimum_cells_per_sample]
  object <- subset(object, subset = sample %in% retained_samples)
  if (ncol(object) == 0L) {
    stop(stage, ": no cells after sample-size filter.")
  }
  assert_no_cc(object$sample, paste0(stage, " subcluster input"))
  object[["RNA"]] <- split(object[["RNA"]], f = object$sample)
  object <- FindVariableFeatures(object, verbose = FALSE)
  object <- ScaleData(object, verbose = FALSE)

  set_recorded_seed(config, stage, "RunPCA")
  object <- RunPCA(
    object,
    npcs = pca_components,
    seed.use = 42L,
    verbose = FALSE
  )
  if (method == "rpca") {
    set_recorded_seed(config, stage, "RPCAIntegration")
    object <- IntegrateLayers(
      object = object,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.rpca",
      dims = dimensions,
      k.weight = k_weight,
      verbose = FALSE
    )
    reduction <- "integrated.rpca"
  } else if (method == "harmony") {
    set_recorded_seed(config, stage, "HarmonyIntegration")
    object <- IntegrateLayers(
      object = object,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      group.by.vars = "sample",
      theta = 3,
      max_iter = 30,
      verbose = FALSE
    )
    reduction <- "harmony"
  } else {
    stop("Unknown integration method: ", method)
  }
  object[["RNA"]] <- JoinLayers(object[["RNA"]])

  set_recorded_seed(config, stage, "FindNeighbors")
  object <- FindNeighbors(
    object,
    reduction = reduction,
    dims = dimensions,
    nn.method = "annoy",
    n.trees = 50L,
    annoy.metric = "euclidean",
    verbose = FALSE
  )
  set_recorded_seed(config, stage, "FindClusters")
  object <- FindClusters(
    object,
    resolution = resolution,
    algorithm = 1L,
    n.start = 10L,
    n.iter = 10L,
    random.seed = 0L,
    verbose = FALSE
  )
  set_recorded_seed(config, stage, "RunUMAP")
  object <- RunUMAP(
    object,
    reduction = reduction,
    dims = dimensions,
    umap.method = "uwot",
    n.neighbors = 15L,
    min.dist = 0.05,
    spread = 1.5,
    seed.use = 42L,
    verbose = FALSE
  )
  assert_no_cc(object$sample, paste0(stage, " new pre-curation checkpoint"))
  saveRDS(object, checkpoint_path, compress = FALSE)
  object
}

write_all_cluster_marker_evidence <- function(object, stage) {
  observed <- sort(unique(as.character(object$seurat_clusters)))
  normalized <- LayerData(object, assay = "RNA", layer = "data")
  candidate_genes <- intersect(
    VariableFeatures(object),
    rownames(normalized)
  )
  candidate_matrix <- normalized[candidate_genes, , drop = FALSE]
  cluster_ids <- as.character(object$seurat_clusters)
  all_evidence <- lapply(observed, function(cluster_id) {
    inside <- cluster_ids == cluster_id
    outside <- !inside
    average_inside <- Matrix::rowMeans(
      candidate_matrix[, inside, drop = FALSE]
    )
    average_outside <- Matrix::rowMeans(
      candidate_matrix[, outside, drop = FALSE]
    )
    marker_table <- data.table(
      stage = stage,
      cluster = cluster_id,
      gene = candidate_genes,
      mean_log_normalized_inside = average_inside,
      mean_log_normalized_outside = average_outside,
      mean_difference = average_inside - average_outside,
      detected_fraction_inside = Matrix::rowMeans(
        candidate_matrix[, inside, drop = FALSE] > 0
      ),
      detected_fraction_outside = Matrix::rowMeans(
        candidate_matrix[, outside, drop = FALSE] > 0
      )
    )
    head(
      marker_table[
        order(-mean_difference, -detected_fraction_inside)
      ],
      50L
    )
  })
  all_evidence <- rbindlist(all_evidence)
  fwrite(
    all_evidence,
    file.path(
      config$audit_root,
      paste0("scrna_", stage, "_all_cluster_top_markers.tsv")
    ),
    sep = "\t"
  )
  all_evidence
}

annotate_clusters <- function(object, annotation_file, stage) {
  annotation <- read_annotation(config, annotation_file)
  overlay_path <- file.path(
    config$pipeline_dir,
    paste0("annotation_", stage, "_raw_derived_overlay.csv")
  )
  if (file.exists(overlay_path)) {
    assert_production_input(config, overlay_path)
    overlay <- fread(
      overlay_path,
      header = FALSE,
      col.names = c("cluster", "label", "evidence")
    )
    overlay[, cluster := as.character(cluster)]
    if (anyDuplicated(overlay$cluster)) {
      stop(stage, ": curation overlay contains duplicate clusters.")
    }
    annotation <- rbind(
      annotation[!cluster %in% overlay$cluster],
      overlay,
      fill = TRUE
    )
    fwrite(
      overlay,
      file.path(
        config$audit_root,
        paste0("scrna_", stage, "_curation_overlay.tsv")
      ),
      sep = "\t"
    )
  }
  observed <- sort(unique(as.character(object$seurat_clusters)))
  all_evidence <- write_all_cluster_marker_evidence(object, stage)
  missing <- setdiff(observed, annotation$cluster)
  if (length(missing) > 0L) {
    fwrite(
      all_evidence[cluster %in% missing],
      file.path(
        config$audit_root,
        paste0("scrna_", stage, "_unmapped_cluster_markers.tsv")
      ),
      sep = "\t"
    )
    stop(
      stage,
      " raw-derived clusters are not covered by ",
      annotation_file,
      ": ",
      paste(missing, collapse = ", ")
    )
  }
  mapping <- setNames(annotation$label, annotation$cluster)
  object$sub_celltype <- unname(
    mapping[as.character(object$seurat_clusters)]
  )
  object <- subset(object, subset = sub_celltype != "rm")
  object
}

filter_by_marker_rules <- function(object, rules, threshold = 0) {
  retained <- rep(TRUE, ncol(object))
  names(retained) <- colnames(object)
  labels <- as.character(object$sub_celltype)
  for (label in names(rules)) {
    label_cells <- colnames(object)[labels == label]
    genes <- intersect(rules[[label]], rownames(object))
    if (length(label_cells) == 0L || length(genes) == 0L) {
      next
    }
    expression <- FetchData(
      object,
      vars = genes,
      cells = label_cells,
      layer = "data"
    )
    pass <- rowSums(expression > threshold) > 0L
    retained[label_cells] <- pass
  }
  subset(object, cells = names(retained)[retained])
}

filter_umap_isolates <- function(object, quantile_cutoff = 0.99) {
  coordinates <- Embeddings(object, "umap")
  neighbor_count <- min(31L, nrow(coordinates))
  if (neighbor_count < 3L) {
    return(object)
  }
  neighbor_fit <- RANN::nn2(coordinates, k = neighbor_count)
  median_distance <- apply(
    neighbor_fit$nn.dists[, 2:neighbor_count, drop = FALSE],
    1L,
    median
  )
  threshold <- quantile(
    median_distance,
    quantile_cutoff,
    na.rm = TRUE,
    names = FALSE
  )
  subset(
    object,
    cells = colnames(object)[median_distance <= threshold]
  )
}

write_subcluster_audit <- function(object, stage) {
  audit <- data.table(
    stage = stage,
    cluster = as.character(object$seurat_clusters),
    sub_celltype = as.character(object$sub_celltype),
    sample = as.character(object$sample),
    group = as.character(object$group)
  )[
    ,
    .(n_cells = .N, n_samples = uniqueN(sample)),
    by = .(stage, cluster, sub_celltype, group)
  ]
  fwrite(
    audit,
    file.path(config$audit_root, paste0("scrna_", stage, "_annotation.tsv")),
    sep = "\t"
  )
}

update_parent_state <- function(atlas_object, parent, subobject) {
  parent_cells <- colnames(atlas_object)[
    as.character(atlas_object$celltype) == parent
  ]
  keep_cells <- intersect(parent_cells, colnames(subobject))
  remove_cells <- setdiff(parent_cells, keep_cells)
  if (length(remove_cells) > 0L) {
    atlas_object <- subset(
      atlas_object,
      cells = setdiff(colnames(atlas_object), remove_cells)
    )
  }
  labels <- setNames(
    as.character(subobject$sub_celltype),
    colnames(subobject)
  )
  atlas_object$sub_celltype[keep_cells] <- labels[keep_cells]
  atlas_object
}

message("Reconstructing epithelial major-state curation.")
epithelial <- subset(atlas, subset = celltype == "Cholangiocyte")
epithelial <- run_subcluster(
  epithelial,
  stage = "epithelial",
  method = "rpca",
  dimensions = 1:15,
  resolution = 0.3,
  minimum_cells_per_sample = 30L,
  k_weight = 30L,
  pca_components = 30L
)
invisible(write_all_cluster_marker_evidence(epithelial, "epithelial"))
epithelial_overlay_path <- file.path(
  config$pipeline_dir,
  "annotation_epithelial_raw_derived_overlay.csv"
)
assert_production_input(config, epithelial_overlay_path)
epithelial_annotation <- fread(
  epithelial_overlay_path,
  header = FALSE,
  col.names = c("cluster", "label", "evidence")
)
epithelial_annotation[, cluster := as.character(cluster)]
if (anyDuplicated(epithelial_annotation$cluster)) {
  stop("Epithelial curation overlay contains duplicate clusters.")
}
epithelial_map <- setNames(
  epithelial_annotation$label,
  epithelial_annotation$cluster
)
observed_epithelial <- unique(as.character(epithelial$seurat_clusters))
if (!setequal(observed_epithelial, names(epithelial_map))) {
  stop(
    "Epithelial raw-derived clusters do not exactly match the expert overlay; ",
    "missing=",
    paste(setdiff(observed_epithelial, names(epithelial_map)), collapse = ","),
    "; extra=",
    paste(setdiff(names(epithelial_map), observed_epithelial), collapse = ",")
  )
}
epithelial$celltype_new <- unname(
  epithelial_map[as.character(epithelial$seurat_clusters)]
)
epithelial_audit <- data.table(
  cluster = epithelial_annotation$cluster,
  expert_label = epithelial_annotation$label,
  evidence = epithelial_annotation$evidence,
  source = "raw-derived marker-based expert curation",
  historical_cluster_ids_used = FALSE,
  status = "observed"
)
fwrite(
  epithelial_audit,
  file.path(config$audit_root, "scrna_epithelial_curation.tsv"),
  sep = "\t"
)
saveRDS(epithelial, config$atlas_paths$epithelial, compress = FALSE)
epithelial_keep <- colnames(epithelial)[epithelial$celltype_new != "rm"]
epithelial_remove <- setdiff(
  colnames(atlas)[as.character(atlas$celltype) == "Cholangiocyte"],
  epithelial_keep
)
atlas <- subset(
  atlas,
  cells = setdiff(colnames(atlas), epithelial_remove)
)
atlas$celltype[epithelial_keep] <- epithelial$celltype_new[
  match(epithelial_keep, colnames(epithelial))
]
atlas$sub_celltype[epithelial_keep] <- atlas$celltype[epithelial_keep]
rm(epithelial)
invisible(gc())

message("Reconstructing T-cell states.")
t_cell <- subset(atlas, subset = celltype == "T")
t_cell <- run_subcluster(
  t_cell,
  stage = "t_cell",
  method = "rpca",
  dimensions = 1:20,
  resolution = 0.5,
  minimum_cells_per_sample = 30L,
  k_weight = 30L,
  pca_components = 25L
)
t_cell <- annotate_clusters(t_cell, "annotation_t.csv", "t_cell")
t_cell <- filter_by_marker_rules(
  t_cell,
  list(
    `CD4+Tn` = c("IL7R", "LEF1"),
    MAIT = c("KLRB1", "SLC4A10"),
    `CD8+Tem` = c("GZMK", "CD8A", "CD8B"),
    `CD8+Temra` = c("GZMB", "GNLY", "FGFBP2"),
    gdT = c("TRDC", "TRGC1", "TRGC2"),
    `CD8+Trm` = c("ZNF683", "XCL1", "XCL2"),
    Treg = c("FOXP3", "CTLA4")
  )
)
t_cell <- filter_umap_isolates(t_cell)
write_subcluster_audit(t_cell, "t_cell")
saveRDS(t_cell, config$atlas_paths$t_cell, compress = FALSE)
atlas <- update_parent_state(atlas, "T", t_cell)
rm(t_cell)
invisible(gc())

message("Reconstructing macrophage states.")
macrophage <- subset(atlas, subset = celltype == "Macrophage")
macrophage <- run_subcluster(
  macrophage,
  stage = "macrophage",
  method = "harmony",
  dimensions = 1:20,
  resolution = 0.4,
  minimum_cells_per_sample = 30L,
  pca_components = 25L
)
macrophage <- annotate_clusters(
  macrophage,
  "annotation_mp.csv",
  "macrophage"
)
macrophage <- filter_by_marker_rules(
  macrophage,
  list(
    `SPP1+M_` = "SPP1",
    `TREM2+M_` = "TREM2",
    `IL1B+M_` = "IL1B",
    `CXCL10+M_` = c("CXCL10", "CXCL9"),
    Kupffer = c("CD5L", "MARCO"),
    Monocytes = c("LYZ", "FCGR1A")
  ),
  threshold = 0.5
)
macrophage <- filter_umap_isolates(macrophage)
write_subcluster_audit(macrophage, "macrophage")
saveRDS(macrophage, config$atlas_paths$macrophage, compress = FALSE)
atlas <- update_parent_state(atlas, "Macrophage", macrophage)
rm(macrophage)
invisible(gc())

message("Reconstructing fibroblast states.")
fibroblast <- subset(atlas, subset = celltype == "Fibroblast")
fibroblast <- run_subcluster(
  fibroblast,
  stage = "fibroblast",
  method = "harmony",
  dimensions = 1:15,
  resolution = 0.3,
  minimum_cells_per_sample = 30L,
  pca_components = 20L
)
fibroblast <- annotate_clusters(
  fibroblast,
  "annotation_fib.csv",
  "fibroblast"
)
fibroblast <- filter_by_marker_rules(
  fibroblast,
  list(
    Pericyte = c("RGS5", "ABCC9", "KCNJ8"),
    SMC = c("MYH11", "CNN1", "MYL9"),
    HSC = c("COL1A1", "CTHRC1", "LUM")
  )
)
fibroblast <- filter_umap_isolates(fibroblast)
write_subcluster_audit(fibroblast, "fibroblast")
saveRDS(fibroblast, config$atlas_paths$fibroblast, compress = FALSE)
atlas <- update_parent_state(atlas, "Fibroblast", fibroblast)
rm(fibroblast)
invisible(gc())

message("Reconstructing B-cell states.")
b_cell <- subset(atlas, subset = celltype == "B")
b_cell <- run_subcluster(
  b_cell,
  stage = "b_cell",
  method = "harmony",
  dimensions = 1:15,
  resolution = 0.3,
  minimum_cells_per_sample = 30L,
  pca_components = 20L
)
b_cell <- annotate_clusters(b_cell, "annotation_b.csv", "b_cell")
b_cell <- filter_by_marker_rules(
  b_cell,
  list(
    naiveB = c("CD79A", "MS4A1", "CD79B"),
    memoryB = c("CD79A", "MS4A1", "CD79B")
  )
)
b_cell <- filter_umap_isolates(b_cell)
write_subcluster_audit(b_cell, "b_cell")
saveRDS(b_cell, config$atlas_paths$b_cell, compress = FALSE)
atlas <- update_parent_state(atlas, "B", b_cell)
rm(b_cell)
invisible(gc())

message("Reconstructing initial hepatocyte states and v5 cholangiocyte curation.")
hepatocyte_initial <- subset(atlas, subset = celltype == "Hepatocyte")
hepatocyte_initial <- run_subcluster(
  hepatocyte_initial,
  stage = "hepatocyte_initial",
  method = "rpca",
  dimensions = 1:50,
  resolution = 0.5,
  minimum_cells_per_sample = 51L,
  k_weight = 30L,
  pca_components = 55L
)
hepatocyte_initial <- annotate_clusters(
  hepatocyte_initial,
  "annotation_hep.csv",
  "hepatocyte_initial"
)
hepatocyte_initial <- filter_umap_isolates(hepatocyte_initial)
write_subcluster_audit(hepatocyte_initial, "hepatocyte_initial")
saveRDS(
  hepatocyte_initial,
  config$atlas_paths$hepatocyte_initial,
  compress = FALSE
)

cholangiocyte_like_states <- "Cholangiocyte"
cholangiocyte_like_cells <- colnames(hepatocyte_initial)[
  as.character(hepatocyte_initial$sub_celltype) %in%
    cholangiocyte_like_states
]
clean_hepatocyte_cells <- setdiff(
  colnames(hepatocyte_initial),
  cholangiocyte_like_cells
)
fwrite(
  data.table(
    cell = colnames(hepatocyte_initial),
    initial_state = as.character(hepatocyte_initial$sub_celltype),
    v5_action = ifelse(
      colnames(hepatocyte_initial) %in% cholangiocyte_like_cells,
      "reclassify_as_Cholangiocyte",
      "retain_for_clean_hepatocyte_RPCA"
    ),
    curation_source = paste(
      "raw-derived marker-based initial hepatocyte curation;",
      "only clusters explicitly labelled Cholangiocyte are reclassified"
    )
  ),
  file.path(config$audit_root, "scrna_hepatocyte_v5_curation.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

hepatocyte_clean <- subset(
  hepatocyte_initial,
  cells = clean_hepatocyte_cells
)
hepatocyte_clean$sub_celltype <- NULL
hepatocyte_clean <- run_subcluster(
  hepatocyte_clean,
  stage = "hepatocyte_clean",
  method = "rpca",
  dimensions = 1:50,
  resolution = 0.5,
  minimum_cells_per_sample = 51L,
  k_weight = 30L,
  pca_components = 55L
)
invisible(
  write_all_cluster_marker_evidence(hepatocyte_clean, "hepatocyte_clean")
)

clean_overlay_path <- file.path(
  config$pipeline_dir,
  "annotation_hepatocyte_clean_raw_derived_overlay.csv"
)
assert_production_input(config, clean_overlay_path)
clean_annotation <- fread(
  clean_overlay_path,
  header = FALSE,
  col.names = c("cluster", "action", "evidence")
)
clean_annotation[, cluster := as.character(cluster)]
observed_clean_clusters <- sort(
  unique(as.character(hepatocyte_clean$seurat_clusters))
)
if (
  anyDuplicated(clean_annotation$cluster) ||
    !setequal(observed_clean_clusters, clean_annotation$cluster)
) {
  stop(
    "Clean hepatocyte clusters do not exactly match the raw-derived overlay."
  )
}
fwrite(
  clean_annotation[
    ,
    `:=`(
      curation_source = "raw-derived marker-based expert curation",
      historical_cluster_ids_used = FALSE
    )
  ],
  file.path(
    config$audit_root,
    "scrna_hepatocyte_clean_curation_overlay.tsv"
  ),
  sep = "\t"
)
clean_action <- setNames(
  clean_annotation$action,
  clean_annotation$cluster
)
hepatocyte_clean$curation_action <- unname(
  clean_action[as.character(hepatocyte_clean$seurat_clusters)]
)
hepatocyte_clean <- subset(
  hepatocyte_clean,
  subset = curation_action == "retain"
)

hep_signatures <- list(
  Hep_C1 = c("PON3", "UGT2B7", "SLC27A5", "F10", "PON1", "PLG"),
  Hep_C2 = c("FTL", "TXN", "UQCRB", "ATP5J2", "FTH1", "APOC1"),
  Hep_C3 = c("CYP2C9", "CYP3A4", "FGA", "FGB", "FGG", "TTR", "ALB"),
  Hep_C4 = c("NDUFA4L2", "ADM", "NDRG1", "PLOD2", "EGLN3", "BNIP3L"),
  Hep_C5 = c("TP53INP1", "CYP8B1", "G6PC", "PDK4", "GPAM", "IGF2"),
  Hep_C6 = c("BIRC5", "PTTG1", "CDKN3", "CENPF", "MKI67", "TOP2A"),
  Hep_C7 = c("GINS2", "CLSPN", "DTL", "MCM10", "CDC6", "TYMS", "MCM4")
)
observed_clean_clusters <- sort(
  unique(as.character(hepatocyte_clean$seurat_clusters))
)
if (length(observed_clean_clusters) < length(hep_signatures)) {
  stop(
    "Clean hepatocyte RPCA produced ",
    length(observed_clean_clusters),
    " retained clusters; seven-signature expert assignment requires at least 7. ",
    "This is an explicit curation boundary."
  )
}
signature_genes <- unique(unlist(hep_signatures))
signature_genes <- intersect(signature_genes, rownames(hepatocyte_clean))
normalized_clean <- LayerData(
  hepatocyte_clean,
  assay = "RNA",
  layer = "data"
)
clean_cluster_ids <- as.character(hepatocyte_clean$seurat_clusters)
average <- vapply(
  observed_clean_clusters,
  function(cluster_id) {
    Matrix::rowMeans(
      normalized_clean[
        signature_genes,
        clean_cluster_ids == cluster_id,
        drop = FALSE
      ]
    )
  },
  numeric(length(signature_genes))
)
rownames(average) <- signature_genes
colnames(average) <- observed_clean_clusters
rm(normalized_clean)
invisible(gc())
scaled <- t(scale(t(as.matrix(average))))
scaled[!is.finite(scaled)] <- 0
score_matrix <- sapply(
  hep_signatures,
  function(genes) {
    available <- intersect(genes, rownames(scaled))
    if (length(available) == 0L) {
      return(rep(-Inf, ncol(scaled)))
    }
    colMeans(scaled[available, , drop = FALSE])
  }
)
rownames(score_matrix) <- colnames(scaled)
assignment_scores <- t(score_matrix)
assignment_scores <- assignment_scores - min(assignment_scores)
core_assignment <- clue::solve_LSAP(
  assignment_scores,
  maximum = TRUE
)
core_clusters <- rownames(score_matrix)[as.integer(core_assignment)]
clean_mapping <- setNames(colnames(score_matrix), core_clusters)
extra_clusters <- setdiff(rownames(score_matrix), core_clusters)
if (length(extra_clusters) > 0L) {
  extra_mapping <- vapply(
    extra_clusters,
    function(cluster_id) {
      colnames(score_matrix)[which.max(score_matrix[cluster_id, ])]
    },
    character(1)
  )
  clean_mapping <- c(clean_mapping, extra_mapping)
}
hepatocyte_clean$sub_celltype <- unname(
  clean_mapping[as.character(hepatocyte_clean$seurat_clusters)]
)
hepatocyte_clean <- filter_umap_isolates(hepatocyte_clean)
hepatocyte_clean$celltype <- "Hepatocyte"
write_subcluster_audit(hepatocyte_clean, "hepatocyte_clean")
score_audit <- as.data.table(
  score_matrix,
  keep.rownames = "cluster"
)
score_audit <- melt(
  score_audit,
  id.vars = "cluster",
  variable.name = "state",
  value.name = "signature_score"
)
score_audit[
  ,
  `:=`(
    assigned = clean_mapping[as.character(cluster)] == as.character(state),
    assignment_mode = ifelse(
      as.character(cluster) %in% core_clusters,
      "one_to_one_core",
      "argmax_merge"
    )
  )
]
fwrite(
  score_audit,
  file.path(config$audit_root, "scrna_hepatocyte_signature_assignment.tsv"),
  sep = "\t"
)
saveRDS(
  hepatocyte_clean,
  config$atlas_paths$hepatocyte,
  compress = FALSE
)

all_hepatocyte_parent <- colnames(atlas)[
  as.character(atlas$celltype) == "Hepatocyte"
]
atlas$celltype[cholangiocyte_like_cells] <- "Cholangiocyte"
atlas$sub_celltype[cholangiocyte_like_cells] <- "Cholangiocyte"
hepatocyte_keep <- colnames(hepatocyte_clean)
hepatocyte_remove <- setdiff(
  all_hepatocyte_parent,
  c(cholangiocyte_like_cells, hepatocyte_keep)
)
if (length(hepatocyte_remove) > 0L) {
  atlas <- subset(
    atlas,
    cells = setdiff(colnames(atlas), hepatocyte_remove)
  )
}
hep_labels <- setNames(
  as.character(hepatocyte_clean$sub_celltype),
  hepatocyte_keep
)
atlas$celltype[hepatocyte_keep] <- "Hepatocyte"
atlas$sub_celltype[hepatocyte_keep] <- hep_labels[hepatocyte_keep]
rm(hepatocyte_initial, hepatocyte_clean)
invisible(gc())

atlas$group <- factor(
  as.character(atlas$group),
  levels = config$group_levels
)
assert_no_cc(atlas$sample, "final raw-derived atlas")
saveRDS(atlas, config$atlas_paths$final, compress = FALSE)

metadata <- as.data.table(atlas[[]], keep.rownames = "cell")
metadata[
  ,
  `:=`(
    group = as.character(group),
    celltype = as.character(celltype),
    sub_celltype = as.character(sub_celltype),
    sample = as.character(sample)
  )
]
sample_summary <- metadata[
  ,
  .(n_cells = .N),
  by = .(sample, group, source, celltype, sub_celltype)
]
sample_summary[
  ,
  proportion_within_sample := n_cells / sum(n_cells),
  by = sample
]
fwrite(
  sample_summary,
  file.path(panel_dir, "figure2_cell_state_frequency_by_sample.tsv"),
  sep = "\t"
)

group_summary <- metadata[
  ,
  .(n_cells = .N),
  by = .(group, celltype, sub_celltype)
]
group_summary[
  ,
  proportion_within_group := n_cells / sum(n_cells),
  by = group
]
fwrite(
  group_summary,
  file.path(panel_dir, "figure2_composition_values.tsv"),
  sep = "\t"
)

stromal_states <- c("Endothelial", "LSEC", "Pericyte", "SMC", "HSC")
immune_states <- c(
  "CD4+Tn", "MAIT", "CD8+Tem", "CD8+Temra", "gdT", "CD8+Trm",
  "Treg", "NK", "naiveB", "memoryB", "Plasma", "Kupffer",
  "IL1B+M_", "TREM2+M_", "Monocytes", "SPP1+M_", "CXCL10+M_",
  "cDC1", "cDC2", "pDC", "Mast", "Neutrophil"
)
roe_metadata <- metadata[
  sub_celltype %in% c(stromal_states, immune_states)
]
contingency <- table(roe_metadata$sub_celltype, roe_metadata$group)
expected <- outer(rowSums(contingency), colSums(contingency)) /
  sum(contingency)
roe <- contingency / expected
roe_table <- as.data.table(as.table(roe))
setnames(roe_table, c("sub_celltype", "group", "roe"))
roe_table[
  ,
  lineage := ifelse(sub_celltype %in% stromal_states, "Stromal", "Immune")
]
fwrite(
  roe_table,
  file.path(panel_dir, "figure2_stromal_immune_roe.tsv"),
  sep = "\t"
)

embedding <- as.data.table(
  Embeddings(atlas, "umap"),
  keep.rownames = "cell"
)
setnames(embedding, names(embedding)[2:3], c("UMAP_1", "UMAP_2"))
embedding <- merge(
  embedding,
  metadata[
    ,
    .(cell, sample, group, source, celltype, sub_celltype)
  ],
  by = "cell",
  all.x = TRUE,
  sort = FALSE
)
fwrite(
  embedding,
  file.path(panel_dir, "figure2_atlas_umap_metadata.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
message(
  "PASS: raw-derived subclusters, final atlas, Figure 2 composition, RoE, ",
  "and UMAP metadata completed."
)
