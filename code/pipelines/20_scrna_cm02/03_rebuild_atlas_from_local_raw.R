#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(data.table)
  library(future)
  library(Matrix)
  library(org.Hs.eg.db)
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
set_recorded_seed(config, "raw_ingest", "matrix_construction")
uncurated_path <- config$atlas_paths$uncurated
resume_uncurated <- resume_enabled("LIUF_RESUME_ATLAS")

for (path in unlist(config$raw_paths, use.names = FALSE)) {
  assert_production_input(config, path)
}

observed_internal <- sort(
  basename(list.dirs(
    config$raw_paths$internal_root,
    recursive = FALSE,
    full.names = TRUE
  ))
)
if (!identical(observed_internal, sort(config$internal_samples))) {
  stop(
    "Internal raw directories must be exactly HCC1-HCC3/NC1-NC3; observed: ",
    paste(observed_internal, collapse = ", ")
  )
}
assert_no_cc(observed_internal, "internal raw directories")

collapse_ensembl_to_symbol <- function(counts) {
  feature_names <- rownames(counts)
  ensembl <- sub("\\..*$", "", feature_names)
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  symbols <- unname(symbols)
  keep <- !is.na(symbols) & nzchar(symbols) & !duplicated(symbols)
  output <- counts[keep, , drop = FALSE]
  rownames(output) <- symbols[keep]
  output
}

read_internal_sample <- function(sample_id) {
  matrix_dir <- file.path(
    config$raw_paths$internal_root,
    sample_id,
    "filtered_feature_bc_matrix"
  )
  counts <- Read10X(
    data.dir = matrix_dir,
    gene.column = 1L,
    unique.features = TRUE
  )
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop("No Gene Expression matrix for ", sample_id)
    }
    counts <- counts[["Gene Expression"]]
  }
  counts <- collapse_ensembl_to_symbol(counts)
  suffix <- unname(config$internal_suffix[[sample_id]])
  colnames(counts) <- paste0(colnames(counts), "_", suffix)
  object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 0L,
    min.features = 0L
  )
  object$sample <- sample_id
  object$group <- ifelse(startsWith(sample_id, "HCC"), "HCC", "Adjacent")
  object$source <- "LiuF"
  object$library <- "internal"
  object$gsm <- NA_character_
  object
}

read_gse136103 <- function(directory) {
  matrix_files <- list.files(
    directory,
    pattern = "_matrix\\.mtx\\.gz$",
    full.names = TRUE
  )
  library_table <- data.frame(
    matrix = matrix_files,
    base = sub("_matrix\\.mtx\\.gz$", "", basename(matrix_files)),
    stringsAsFactors = FALSE
  )
  library_table$gsm <- sub("_.*", "", library_table$base)
  library_table$sample <- sub(
    "^GSM[0-9]+_([^_]+).*",
    "\\1",
    library_table$base
  )
  library_table$library <- sub(
    "^GSM[0-9]+_[^_]+_?",
    "",
    library_table$base
  )
  library_table <- library_table[
    grepl("^(healthy|cirrhotic)", library_table$sample),
    ,
    drop = FALSE
  ]
  if (nrow(library_table) == 0L) {
    stop("No healthy/cirrhotic GSE136103 libraries found.")
  }
  objects <- lapply(seq_len(nrow(library_table)), function(index) {
    record <- library_table[index, ]
    prefix <- file.path(directory, record$base)
    counts <- ReadMtx(
      mtx = paste0(prefix, "_matrix.mtx.gz"),
      cells = paste0(prefix, "_barcodes.tsv.gz"),
      features = paste0(prefix, "_genes.tsv.gz"),
      feature.column = 2L
    )
    object <- CreateSeuratObject(
      counts = counts,
      project = record$sample,
      min.cells = 0L,
      min.features = 0L
    )
    object$sample <- record$sample
    object$library <- record$library
    object$gsm <- record$gsm
    object$group <- ifelse(
      startsWith(record$sample, "healthy"),
      "Healthy",
      "Cirrhosis"
    )
    object$source <- "GSE136103"
    object
  })
  names(objects) <- library_table$base
  output <- merge(
    objects[[1L]],
    y = objects[-1L],
    add.cell.ids = paste0(library_table$sample, "_", library_table$library)
  )
  JoinLayers(output)
}

read_gse115469 <- function(path) {
  input <- fread(path, data.table = FALSE)
  gene_names <- input[[1L]]
  counts <- Matrix(as.matrix(input[, -1L, drop = FALSE]), sparse = TRUE)
  rownames(counts) <- gene_names
  object <- CreateSeuratObject(
    counts = counts,
    project = "GSE115469",
    min.cells = 0L,
    min.features = 0L
  )
  object$sample <- sub("_.*", "", colnames(object))
  object$library <- "TLH"
  object$gsm <- NA_character_
  object$group <- "Healthy"
  object$source <- "GSE115469"
  object
}

read_gse149614 <- function(count_path, metadata_path) {
  metadata <- fread(metadata_path, data.table = FALSE)
  metadata <- metadata[metadata$site %in% c("Tumor", "Normal"), , drop = FALSE]
  input <- fread(count_path, data.table = FALSE)
  gene_names <- input[[1L]]
  input <- input[, -1L, drop = FALSE]
  selected <- colnames(input) %in% metadata$Cell
  input <- input[, selected, drop = FALSE]
  counts <- Matrix(as.matrix(input), sparse = TRUE)
  rownames(counts) <- gene_names
  object <- CreateSeuratObject(
    counts = counts,
    project = "GSE149614",
    min.cells = 0L,
    min.features = 0L
  )
  metadata <- metadata[match(colnames(object), metadata$Cell), , drop = FALSE]
  if (anyNA(metadata$Cell)) {
    stop("GSE149614 metadata alignment failed.")
  }
  object$sample <- metadata$sample
  object$patient <- metadata$patient
  object$stage <- metadata$stage
  object$virus <- metadata$virus
  object$library <- metadata$site
  object$group <- ifelse(metadata$site == "Tumor", "HCC", "Adjacent")
  object$source <- "GSE149614"
  object
}

if (resume_uncurated && file.exists(uncurated_path)) {
  message("Resuming from raw-derived pre-curation atlas checkpoint.")
  combined <- readRDS(uncurated_path)
  assert_no_cc(combined$sample, "pre-curation atlas checkpoint")
} else {
  message("Reading six internal Cell Ranger matrices.")
  internal_objects <- lapply(config$internal_samples, read_internal_sample)
  names(internal_objects) <- config$internal_samples
  internal <- merge(internal_objects[[1L]], y = internal_objects[-1L])
  internal <- JoinLayers(internal)
  assert_no_cc(internal$sample, "raw-derived internal object")
  rm(internal_objects)
  invisible(gc())

  message("Reading three public scRNA datasets from local files.")
  gse136103 <- read_gse136103(config$raw_paths$gse136103)
  gse115469 <- read_gse115469(config$raw_paths$gse115469)
  gse149614 <- read_gse149614(
    config$raw_paths$gse149614_counts,
    config$raw_paths$gse149614_metadata
  )

  combined <- merge(internal, y = list(gse136103, gse115469, gse149614))
  combined <- JoinLayers(combined)
  combined$group <- factor(
    as.character(combined$group),
    levels = config$group_levels
  )
  assert_no_cc(combined$sample, "merged raw atlas")
  rm(internal, gse136103, gse115469, gse149614)
  invisible(gc())

  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
  combined[["percent.HB"]] <- PercentageFeatureSet(
    combined,
    pattern = "^HB[^P]"
  )
  pre_qc <- data.table(
    sample = as.character(combined$sample),
    group = as.character(combined$group),
    source = as.character(combined$source)
  )[, .(pre_qc_cells = .N), by = .(sample, group, source)]
  combined <- subset(
    combined,
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 5000 &
      percent.mt < 10
  )
  assert_no_cc(combined$sample, "post-QC raw atlas")
  post_qc <- data.table(
    sample = as.character(combined$sample),
    group = as.character(combined$group),
    source = as.character(combined$source)
  )[, .(post_qc_cells = .N), by = .(sample, group, source)]
  qc_audit <- merge(
    pre_qc,
    post_qc,
    by = c("sample", "group", "source"),
    all = TRUE
  )
  qc_audit[is.na(qc_audit)] <- 0L
  fwrite(
    qc_audit,
    file.path(config$audit_root, "scrna_raw_qc_counts.tsv"),
    sep = "\t"
  )
  saveRDS(combined, config$atlas_paths$qc, compress = FALSE)

  message("Running full sample-layer RPCA integration.")
  combined[["RNA"]] <- split(combined[["RNA"]], f = combined$sample)
  combined <- NormalizeData(combined, verbose = FALSE)
  combined <- FindVariableFeatures(combined, verbose = FALSE)
  combined <- ScaleData(combined, verbose = FALSE)

  set_recorded_seed(config, "atlas", "RunPCA")
  combined <- RunPCA(
    combined,
    npcs = 50L,
    seed.use = 42L,
    verbose = FALSE
  )
  set_recorded_seed(config, "atlas", "RPCAIntegration")
  combined <- IntegrateLayers(
    object = combined,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    dims = 1:30,
    verbose = FALSE
  )
  combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

  set_recorded_seed(config, "atlas", "FindNeighbors")
  combined <- FindNeighbors(
    combined,
    reduction = "integrated.rpca",
    dims = 1:30,
    nn.method = "annoy",
    n.trees = 50L,
    annoy.metric = "euclidean",
    verbose = FALSE
  )
  set_recorded_seed(config, "atlas", "FindClusters")
  combined <- FindClusters(
    combined,
    resolution = 0.5,
    algorithm = 1L,
    n.start = 10L,
    n.iter = 10L,
    random.seed = 0L,
    verbose = FALSE
  )
  set_recorded_seed(config, "atlas", "RunUMAP")
  combined <- RunUMAP(
    combined,
    dims = 1:30,
    reduction = "integrated.rpca",
    umap.method = "uwot",
    n.neighbors = 30L,
    min.dist = 0.3,
    seed.use = 42L,
    verbose = FALSE
  )
  assert_no_cc(combined$sample, "new pre-curation atlas checkpoint")
  saveRDS(combined, uncurated_path, compress = FALSE)
}

historical_annotation <- read_annotation(config, "annotation.csv")
overlay_path <- file.path(
  config$pipeline_dir,
  "annotation_global_raw_derived_overlay.csv"
)
if (!file.exists(overlay_path)) {
  stop(
    "Current raw-derived clusters require an explicit global curation overlay: ",
    overlay_path
  )
}
assert_production_input(config, overlay_path)
annotation <- fread(
  overlay_path,
  header = FALSE,
  col.names = c("cluster", "label", "evidence")
)
annotation[, cluster := as.character(cluster)]
if (anyDuplicated(annotation$cluster)) {
  stop("Global raw-derived curation overlay contains duplicate clusters.")
}
observed_clusters <- sort(unique(as.character(combined$seurat_clusters)))
missing_clusters <- setdiff(observed_clusters, annotation$cluster)
extra_clusters <- setdiff(annotation$cluster, observed_clusters)
if (length(missing_clusters) > 0L || length(extra_clusters) > 0L) {
  canonical_markers <- c(
    "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1",
    "S100A8", "S100A9", "FCGR3B", "CSF3R",
    "LYZ", "CD68", "C1QC", "MARCO", "CD163",
    "ALB", "APOA1", "HNF4A", "CYP3A4", "AFP",
    "PECAM1", "VWF", "KDR", "EMCN",
    "EPCAM", "KRT7", "KRT19", "KRT8",
    "COL1A1", "COL3A1", "DCN", "RGS5", "ACTA2",
    "CD79A", "MS4A1", "CD37", "FCER1A", "CD1C",
    "CLEC9A", "BATF3", "MZB1", "JCHAIN", "GZMB",
    "TPSAB1", "KIT", "MKI67"
  )
  canonical_markers <- intersect(canonical_markers, rownames(combined))
  canonical_average <- AverageExpression(
    combined,
    assays = "RNA",
    features = canonical_markers,
    group.by = "seurat_clusters",
    layer = "data",
    verbose = FALSE
  )$RNA
  canonical_long <- as.data.table(
    as.table(as.matrix(canonical_average)),
    keep.rownames = FALSE
  )
  setnames(canonical_long, c("gene", "cluster", "average_expression"))
  fwrite(
    canonical_long,
    file.path(config$audit_root, "scrna_global_cluster_marker_evidence.tsv"),
    sep = "\t"
  )

  normalized <- LayerData(combined, assay = "RNA", layer = "data")
  candidate_genes <- intersect(VariableFeatures(combined), rownames(normalized))
  candidate_matrix <- normalized[candidate_genes, , drop = FALSE]
  cluster_ids <- as.character(combined$seurat_clusters)
  top_marker_rows <- lapply(missing_clusters, function(cluster_id) {
    inside <- cluster_ids == cluster_id
    outside <- !inside
    average_inside <- Matrix::rowMeans(candidate_matrix[, inside, drop = FALSE])
    average_outside <- Matrix::rowMeans(candidate_matrix[, outside, drop = FALSE])
    fraction_inside <- Matrix::rowMeans(
      candidate_matrix[, inside, drop = FALSE] > 0
    )
    fraction_outside <- Matrix::rowMeans(
      candidate_matrix[, outside, drop = FALSE] > 0
    )
    evidence <- data.table(
      cluster = cluster_id,
      gene = candidate_genes,
      mean_log_normalized_inside = average_inside,
      mean_log_normalized_outside = average_outside,
      mean_difference = average_inside - average_outside,
      detected_fraction_inside = fraction_inside,
      detected_fraction_outside = fraction_outside
    )
    head(
      evidence[order(-mean_difference, -detected_fraction_inside)],
      50L
    )
  })
  fwrite(
    rbindlist(top_marker_rows, fill = TRUE),
    file.path(config$audit_root, "scrna_unmapped_cluster_top_markers.tsv"),
    sep = "\t"
  )
  fwrite(
    data.table(
      boundary = "global_cluster_annotation",
      observed_cluster = observed_clusters,
      mapped = observed_clusters %in% annotation$cluster,
      curation_source = "raw-derived expert overlay"
    ),
    file.path(config$audit_root, "scrna_curation_boundary.tsv"),
    sep = "\t"
  )
  stop(
    "Raw-derived global overlay mismatch; missing=",
    paste(missing_clusters, collapse = ","),
    "; extra=",
    paste(extra_clusters, collapse = ",")
  )
}
allowed_labels <- union(
  unique(historical_annotation$label),
  "mixed_curation"
)
unknown_labels <- setdiff(unique(annotation$label), allowed_labels)
if (length(unknown_labels) > 0L) {
  stop(
    "Raw-derived global overlay uses labels outside the historical vocabulary: ",
    paste(unknown_labels, collapse = ", ")
  )
}
fwrite(
  annotation[
    ,
    `:=`(
      curation_source = "raw-derived marker-based expert curation",
      historical_cluster_ids_used = FALSE
    )
  ],
  file.path(config$audit_root, "scrna_global_curation_overlay.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    boundary = "global_cluster_annotation",
    observed_cluster = observed_clusters,
    mapped = TRUE,
    curation_source = "raw-derived expert overlay",
    historical_cluster_ids_used = FALSE
  ),
  file.path(config$audit_root, "scrna_curation_boundary.tsv"),
  sep = "\t"
)
cluster_map <- setNames(annotation$label, annotation$cluster)
combined$celltype <- unname(
  cluster_map[as.character(combined$seurat_clusters)]
)

normalized <- LayerData(combined, assay = "RNA", layer = "data")
apply_mixed_cluster_gate <- function(
    object,
    normalized_matrix,
    cluster_id,
    target_label,
    markers,
    minimum_detected = 2L) {
  cluster_cells <- colnames(object)[
    as.character(object$seurat_clusters) == cluster_id
  ]
  available_markers <- intersect(markers, rownames(normalized_matrix))
  if (length(cluster_cells) == 0L) {
    stop("Mixed-cluster gate cluster is absent: ", cluster_id)
  }
  if (length(available_markers) < minimum_detected) {
    stop(
      "Mixed-cluster gate has insufficient available markers for cluster ",
      cluster_id
    )
  }
  detected <- Matrix::colSums(
    normalized_matrix[
      available_markers,
      cluster_cells,
      drop = FALSE
    ] > 0
  )
  assigned <- ifelse(detected >= minimum_detected, target_label, "rm")
  object$celltype[cluster_cells] <- assigned
  audit <- data.table(
    cell = cluster_cells,
    cluster = cluster_id,
    target_label = target_label,
    detected_marker_count = as.integer(detected),
    minimum_detected = minimum_detected,
    assigned_label = assigned,
    available_markers = paste(available_markers, collapse = ";"),
    rule = paste0(
      "retain as ",
      target_label,
      " when >= ",
      minimum_detected,
      " lineage markers have normalized expression > 0"
    )
  )
  list(object = object, audit = audit)
}

mixed_gate_specs <- list(
  list(
    cluster = "14",
    label = "Hepatocyte",
    markers = c("ALB", "APOA2", "TTR", "FGA", "FGB", "FGG", "CYP3A4", "HNF4A")
  ),
  list(
    cluster = "20",
    label = "Mast",
    markers = c("TPSAB1", "TPSB2", "KIT", "MS4A2", "CPA3", "HDC")
  ),
  list(
    cluster = "27",
    label = "T",
    markers = c("CD3D", "CD3E", "TRAC", "CD247", "TRBC1", "TRBC2")
  )
)
mixed_gate_audit <- vector("list", length(mixed_gate_specs))
for (index in seq_along(mixed_gate_specs)) {
  specification <- mixed_gate_specs[[index]]
  gated <- apply_mixed_cluster_gate(
    combined,
    normalized,
    cluster_id = specification$cluster,
    target_label = specification$label,
    markers = specification$markers,
    minimum_detected = 2L
  )
  combined <- gated$object
  mixed_gate_audit[[index]] <- gated$audit
}
rm(normalized)
invisible(gc())
fwrite(
  rbindlist(mixed_gate_audit),
  file.path(
    config$audit_root,
    "scrna_global_mixed_cluster_cell_gates.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)
if (any(combined$celltype == "mixed_curation")) {
  stop("Mixed-cluster curation left unresolved cells.")
}
combined <- subset(combined, subset = celltype != "rm")
combined$major <- combined$celltype
combined$sub_celltype <- combined$celltype
combined$group <- factor(
  as.character(combined$group),
  levels = config$group_levels
)
assert_no_cc(combined$sample, "annotated raw-derived major atlas")

cluster_audit <- data.table(
  cluster = as.character(combined$seurat_clusters),
  label = as.character(combined$celltype),
  sample = as.character(combined$sample)
)[
  ,
  .(n_cells = .N, n_samples = uniqueN(sample)),
  by = .(cluster, label)
]
cluster_audit[, cluster_numeric := as.integer(cluster)]
setorder(cluster_audit, cluster_numeric)
cluster_audit[, cluster_numeric := NULL]
fwrite(
  cluster_audit,
  file.path(config$audit_root, "scrna_global_cluster_annotation.tsv"),
  sep = "\t"
)
saveRDS(combined, config$atlas_paths$integrated, compress = FALSE)
message(
  "PASS: raw counts -> QC -> full RPCA -> major annotation completed; ",
  ncol(combined),
  " cells retained."
)
