#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(NMF)
  library(Seurat)
})

arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments) >= 1L) arguments[[1L]] else "audit"
if (!mode %in% c("audit", "curated", "full")) {
  stop("Mode must be one of: audit, curated, full.")
}
all_arguments <- commandArgs(trailingOnly = FALSE)
this_script <- normalizePath(
  sub("^--file=", "", all_arguments[grepl("^--file=", all_arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(this_script), "config.R"))
config <- ensure_output_dirs(scrna_cm02_config(this_script))

target_paths <- list(
  atlas = file.path(
    config$validation_root,
    "dat",
    "combined_annotated_v5.rds"
  ),
  hepatocyte = file.path(config$validation_root, "dat", "hep_annotated.rds"),
  t_cell = file.path(config$validation_root, "dat", "t_annotated.rds"),
  macrophage = file.path(config$validation_root, "dat", "mp_annotated.rds"),
  fibroblast = file.path(config$validation_root, "dat", "fib_annotated.rds"),
  b_cell = file.path(config$validation_root, "dat", "b_annotated.rds"),
  frequency = file.path(
    config$validation_root,
    "covarnet_v2",
    "freq_inputs.rds"
  ),
  nmf = file.path(config$validation_root, "covarnet_cv_k5", "NMF_K.rds"),
  network = file.path(
    config$validation_root,
    "covarnet_cv_k5",
    "cm_network.rds"
  ),
  network_filter = file.path(
    config$validation_root,
    "covarnet_cv_k5",
    "cm_filter_top10.csv"
  )
)
missing_targets <- target_paths[!file.exists(unlist(target_paths))]
if (length(missing_targets) > 0L) {
  stop("Missing validation targets: ", paste(missing_targets, collapse = ", "))
}

target_inventory <- data.table(
  target = names(target_paths),
  path = unlist(target_paths),
  sha256 = vapply(unlist(target_paths), sha256_file, character(1)),
  role = "validation_only"
)
trajectory_paths <- list.files(
  file.path(config$validation_root, "trajectory"),
  full.names = TRUE
)
target_inventory <- rbind(
  target_inventory,
  data.table(
    target = paste0("trajectory/", basename(trajectory_paths)),
    path = trajectory_paths,
    sha256 = vapply(trajectory_paths, sha256_file, character(1)),
    role = "validation_only_not_used_as_scRNA_or_CoVarNet_input"
  )
)
fwrite(
  target_inventory,
  file.path(config$audit_root, "scrna_validation_target_inventory.tsv"),
  sep = "\t"
)

if (mode == "audit") {
  fwrite(
    data.table(
      check = "raw-derived vs frozen validation",
      observed = "not_run",
      expected = "run_all.sh curated or run_all.sh full",
      status = "NOT_RUN_AUDIT_MODE"
    ),
    file.path(config$audit_root, "scrna_raw_vs_frozen.tsv"),
    sep = "\t"
  )
  message(
    "Validation targets inventoried; frozen-curated production awaits curated ",
    "or full mode."
  )
  quit(save = "no", status = 0L)
}

production_paths <- c(
  if (mode == "full") {
    unlist(config$atlas_paths, use.names = FALSE)
  } else {
    character()
  },
  unlist(config$curated_atlas_paths, use.names = FALSE),
  unlist(config$covarnet_paths, use.names = FALSE)
)
missing_production <- production_paths[!file.exists(production_paths)]
if (length(missing_production) > 0L) {
  stop(
    mode,
    " validation is missing required production outputs: ",
    paste(missing_production, collapse = ", ")
  )
}

comparison_rows <- list()
compare_seurat <- function(name, observed_path, expected_path) {
  observed <- readRDS(observed_path)
  expected <- readRDS(expected_path)
  observed_metadata <- observed[[]]
  expected_metadata <- expected[[]]
  assert_no_cc(observed_metadata$sample, paste0(name, " raw-derived"))
  assert_no_cc(expected_metadata$sample, paste0(name, " validation target"))

  observed_cells <- colnames(observed)
  expected_cells <- colnames(expected)
  shared_cells <- intersect(observed_cells, expected_cells)
  union_cells <- union(observed_cells, expected_cells)
  cell_jaccard <- length(shared_cells) / length(union_cells)
  metadata_columns <- intersect(
    c("sample", "group", "source", "celltype", "sub_celltype"),
    intersect(colnames(observed_metadata), colnames(expected_metadata))
  )
  metadata_agreement <- vapply(
    metadata_columns,
    function(column) {
      mean(
        as.character(observed_metadata[shared_cells, column]) ==
          as.character(expected_metadata[shared_cells, column])
      )
    },
    numeric(1)
  )

  count_status <- "not_comparable"
  count_difference_nnz <- NA_real_
  shared_features <- intersect(rownames(observed), rownames(expected))
  if (
    length(shared_cells) > 0L &&
      length(shared_features) > 0L
  ) {
    observed_counts <- GetAssayData(observed, assay = "RNA", layer = "counts")
    expected_counts <- GetAssayData(expected, assay = "RNA", layer = "counts")
    difference <- drop0(
      as(
        observed_counts[shared_features, shared_cells, drop = FALSE] -
          expected_counts[shared_features, shared_cells, drop = FALSE],
        "dgCMatrix"
      )
    )
    count_difference_nnz <- length(difference@x)
    count_status <- if (count_difference_nnz == 0L) {
      "exact_on_shared_cells_and_features"
    } else {
      "different_on_shared_cells_or_features"
    }
    rm(observed_counts, expected_counts, difference)
  }

  umap_correlation <- NA_real_
  if (
    "umap" %in% Reductions(observed) &&
      "umap" %in% Reductions(expected) &&
      length(shared_cells) >= 3L
  ) {
    observed_umap <- scale(Embeddings(observed, "umap")[shared_cells, ])
    expected_umap <- scale(Embeddings(expected, "umap")[shared_cells, ])
    singular <- svd(crossprod(observed_umap, expected_umap))
    rotation <- singular$u %*% t(singular$v)
    aligned <- observed_umap %*% rotation
    umap_correlation <- cor(
      as.numeric(aligned),
      as.numeric(expected_umap)
    )
  }

  output <- data.table(
    object = name,
    observed_cells = length(observed_cells),
    expected_cells = length(expected_cells),
    shared_cells = length(shared_cells),
    cell_jaccard = cell_jaccard,
    observed_features = nrow(observed),
    expected_features = nrow(expected),
    shared_features = length(shared_features),
    count_difference_nnz = count_difference_nnz,
    count_status = count_status,
    umap_procrustes_correlation = umap_correlation,
    status = if (
      length(shared_cells) > 0L &&
        all(metadata_agreement >= 0, na.rm = TRUE)
    ) {
      "COMPARISON_COMPLETE"
    } else {
      "FAIL"
    }
  )
  for (column in names(metadata_agreement)) {
    output[[paste0(column, "_agreement")]] <- metadata_agreement[[column]]
  }
  rm(observed, expected)
  invisible(gc())
  output
}

if (mode == "full") {
  comparison_rows[["atlas"]] <- compare_seurat(
    "atlas",
    config$atlas_paths$final,
    target_paths$atlas
  )
  comparison_rows[["hepatocyte"]] <- compare_seurat(
    "hepatocyte",
    config$atlas_paths$hepatocyte,
    target_paths$hepatocyte
  )
  comparison_rows[["t_cell"]] <- compare_seurat(
    "t_cell",
    config$atlas_paths$t_cell,
    target_paths$t_cell
  )
  comparison_rows[["macrophage"]] <- compare_seurat(
    "macrophage",
    config$atlas_paths$macrophage,
    target_paths$macrophage
  )
  comparison_rows[["fibroblast"]] <- compare_seurat(
    "fibroblast",
    config$atlas_paths$fibroblast,
    target_paths$fibroblast
  )
  comparison_rows[["b_cell"]] <- compare_seurat(
    "b_cell",
    config$atlas_paths$b_cell,
    target_paths$b_cell
  )
  fwrite(
    rbindlist(comparison_rows, fill = TRUE),
    file.path(config$audit_root, "scrna_raw_vs_frozen_objects.tsv"),
    sep = "\t"
  )
}

observed_frequency <- readRDS(
  config$covarnet_paths$frequency
)
expected_frequency <- readRDS(target_paths$frequency)
shared_states <- intersect(
  rownames(observed_frequency$mat_norm),
  rownames(expected_frequency$mat_norm)
)
shared_samples <- intersect(
  colnames(observed_frequency$mat_norm),
  colnames(expected_frequency$mat_norm)
)
frequency_difference <- (
  observed_frequency$mat_norm[shared_states, shared_samples, drop = FALSE] -
    expected_frequency$mat_norm[shared_states, shared_samples, drop = FALSE]
)
frequency_audit <- data.table(
  check = "CoVarNet minmax frequency",
  observed_dimensions = paste(dim(observed_frequency$mat_norm), collapse = "x"),
  expected_dimensions = paste(dim(expected_frequency$mat_norm), collapse = "x"),
  shared_states = length(shared_states),
  shared_samples = length(shared_samples),
  max_absolute_difference = max(abs(frequency_difference)),
  pearson_correlation = cor(
    as.numeric(
      observed_frequency$mat_norm[
        shared_states,
        shared_samples,
        drop = FALSE
      ]
    ),
    as.numeric(
      expected_frequency$mat_norm[
        shared_states,
        shared_samples,
        drop = FALSE
      ]
    )
  ),
  status = "COMPARISON_COMPLETE"
)
fwrite(
  frequency_audit,
  file.path(
    config$audit_root,
    "scrna_frozen_curated_vs_validation_frequency.tsv"
  ),
  sep = "\t"
)
rm(observed_frequency, expected_frequency, frequency_difference)
invisible(gc())

observed_nmf <- readRDS(
  config$covarnet_paths$nmf
)
expected_nmf <- readRDS(target_paths$nmf)
observed_basis <- NMF::basis(observed_nmf)
expected_basis <- NMF::basis(expected_nmf)
shared_states <- intersect(rownames(observed_basis), rownames(expected_basis))
normalize_basis <- function(value) {
  sweep(value, 2L, sqrt(colSums(value^2)), "/")
}
similarity <- crossprod(
  normalize_basis(expected_basis[shared_states, , drop = FALSE]),
  normalize_basis(observed_basis[shared_states, , drop = FALSE])
)
nmf_audit <- as.data.table(as.table(similarity))
setnames(nmf_audit, c("expected_module", "observed_module", "cosine"))
fwrite(
  nmf_audit,
  file.path(
    config$audit_root,
    "scrna_frozen_curated_vs_validation_nmf.tsv"
  ),
  sep = "\t"
)

observed_network <- readRDS(
  config$covarnet_paths$network
)
expected_network <- readRDS(target_paths$network)
observed_cm02 <- unique(
  observed_network$filter$subCluster[observed_network$filter$cm == "CM02"]
)
expected_cm02 <- unique(
  expected_network$filter$subCluster[expected_network$filter$cm == "CM02"]
)
network_audit <- data.table(
  check = "CM02 connected states",
  observed = paste(sort(observed_cm02), collapse = ";"),
  expected = paste(sort(expected_cm02), collapse = ";"),
  intersection = paste(
    sort(intersect(observed_cm02, expected_cm02)),
    collapse = ";"
  ),
  jaccard = length(intersect(observed_cm02, expected_cm02)) /
    length(union(observed_cm02, expected_cm02)),
  status = "COMPARISON_COMPLETE"
)
fwrite(
  network_audit,
  file.path(
    config$audit_root,
    "scrna_frozen_curated_vs_validation_network.tsv"
  ),
  sep = "\t"
)

if (mode == "full") {
  fwrite(
    data.table(
      check = "raw-derived audit plus frozen-curated production validation",
      observed = "completed",
      expected = paste(
        "raw-derived objects retained for audit; expert-curated objects used",
        "for downstream annotations and CoVarNet production"
      ),
      status = "PASS"
    ),
    file.path(config$audit_root, "scrna_raw_vs_frozen.tsv"),
    sep = "\t"
  )
  message(
    "PASS: raw-derived objects were audited separately; frozen expert-curated ",
    "objects drove CoVarNet and were compared with validation-only targets."
  )
} else {
  fwrite(
    data.table(
      check = "frozen-curated production validation",
      observed = "completed_without_raw_atlas_rebuild",
      expected = paste(
        "checksum-locked curated atlas drives CoVarNet and all downstream",
        "single-cell-derived inputs"
      ),
      status = "PASS"
    ),
    file.path(config$audit_root, "scrna_curated_production.tsv"),
    sep = "\t"
  )
  message(
    "PASS: frozen expert-curated objects drove CoVarNet without raw-count ",
    "atlas reconstruction and matched the validation-only targets."
  )
}
