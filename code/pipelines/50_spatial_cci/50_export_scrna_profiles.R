#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

source(file.path(dirname(normalizePath(
  sub(
    "^--file=",
    "",
    commandArgs(trailingOnly = FALSE)[
      grepl("^--file=", commandArgs(trailingOnly = FALSE))
    ]
  ),
  mustWork = TRUE
)), "common.R"))

seed <- configure_reproducibility()
root <- project_root()
result_root <- spatial_result_root()
arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[grepl("^--output-root=", arguments)]
if (length(output_argument) > 1L) {
  stop("At most one --output-root argument is allowed.")
}
atlas_path <- assert_production_path(
  file.path(
    root,
    "data",
    "curation",
    "scRNA",
    "combined_annotated_v5.rds"
  )
)
hep_path <- assert_production_path(
  file.path(root, "data", "curation", "scRNA", "hep_annotated.rds")
)
curation_audit_path <- assert_production_path(
  file.path(root, "audit", "scrna_frozen_curation_objects.tsv")
)
ligrec_path <- assert_production_path(
  file.path(result_root, "commot", "cellchat_ligrec_full.csv")
)
msigdb_zip_path <- assert_production_path(
  file.path(
    root,
    "data",
    "reference",
    "MSigDB",
    "msigdb_v2023.2.Hs_files_to_download_locally.zip"
  )
)
expected_msigdb_sha256 <-
  "0aa92b11ad56204905ced89581a65affe7ce422c11c9660851c2ec6829a1c33a"
output_root <- if (length(output_argument) == 1L) {
  normalizePath(
    sub("^--output-root=", "", output_argument),
    mustWork = FALSE
  )
} else {
  file.path(result_root, "cci")
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
if (!all(file.exists(c(
  atlas_path,
  hep_path,
  curation_audit_path,
  ligrec_path,
  msigdb_zip_path
)))) {
  stop("Frozen-curated scRNA objects, CellChat DB, and MSigDB are required.")
}
if (!identical(sha256_file(msigdb_zip_path), expected_msigdb_sha256)) {
  stop("MSigDB v2023.2 archive checksum mismatch.")
}
curation_audit <- fread(curation_audit_path)
validate_scrna_curation_audit(curation_audit, atlas_path)
hep_audit <- curation_audit[object == "hepatocyte"]
if (nrow(hep_audit) != 1L ||
    !identical(sha256_file(hep_path), as.character(hep_audit$sha256[[1L]]))) {
  stop("Frozen-curated hepatocyte object changed after its audit was written.")
}

atlas <- readRDS(atlas_path)
hep <- readRDS(hep_path)
if (!inherits(atlas, "Seurat") || !inherits(hep, "Seurat")) {
  stop("Frozen-curated scRNA objects are not Seurat objects.")
}
DefaultAssay(atlas) <- "RNA"
atlas <- tryCatch(JoinLayers(atlas), error = function(error) atlas)
DefaultAssay(hep) <- "RNA"
hep <- tryCatch(JoinLayers(hep), error = function(error) hep)
metadata <- atlas[[]]
if (!all(c("sub_celltype", "group", "sample", "source") %in%
         colnames(metadata))) {
  stop("Frozen-curated scRNA atlas lacks required annotations.")
}
if (any(as.character(metadata$sample) %in% c("CC1", "CC2", "CC3"))) {
  stop("Frozen-curated scRNA atlas contains prohibited CC1-CC3 samples.")
}
labels <- as.character(metadata$sub_celltype)
hep_labels <- setNames(as.character(hep$sub_celltype), colnames(hep))
hep_match <- match(colnames(atlas), names(hep_labels))
labels[!is.na(hep_match)] <- hep_labels[hep_match[!is.na(hep_match)]]
labels <- clean_rctd_label(labels)
counts_all <- GetAssayData(atlas, assay = "RNA", layer = "counts")
valid <- !is.na(labels) &
  labels != "" &
  labels != "rm"
label_counts <- table(labels[valid])
valid_types <- names(label_counts)[label_counts >= 25L]
keep <- valid & labels %in% valid_types
labels <- droplevels(factor(labels[keep]))
counts <- counts_all[, keep, drop = FALSE]

ligrec <- fread(ligrec_path)
complexes <- c(as.character(ligrec$ligand), as.character(ligrec$receptor))
lr_genes <- sort(unique(unlist(strsplit(complexes, "_", fixed = TRUE))))
genes <- intersect(lr_genes, rownames(counts))
library_size <- Matrix::colSums(counts)
library_size[library_size <= 0] <- 1
lr_counts <- as(counts[genes, , drop = FALSE], "dgCMatrix")
celltypes <- levels(labels)
mean_cpm <- matrix(
  0,
  nrow = length(genes),
  ncol = length(celltypes),
  dimnames = list(genes, celltypes)
)
detection <- mean_cpm
for (celltype in celltypes) {
  indices <- which(labels == celltype)
  cell_counts <- lr_counts[, indices, drop = FALSE]
  scaling <- Matrix::Diagonal(x = 1e4 / library_size[indices])
  mean_cpm[, celltype] <- Matrix::rowMeans(cell_counts %*% scaling)
  detection[, celltype] <- Matrix::rowMeans(cell_counts > 0)
}
write.csv(
  data.frame(gene = rownames(mean_cpm), mean_cpm, check.names = FALSE),
  file.path(output_root, "rctd_reference_lr_mean_cpm.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene = rownames(detection), detection, check.names = FALSE),
  file.path(output_root, "rctd_reference_lr_pct.csv"),
  row.names = FALSE
)
fwrite(
  data.table(
    celltype = celltypes,
    n_cells = as.integer(table(labels)[celltypes])
  ),
  file.path(output_root, "rctd_reference_celltype_counts.tsv"),
  sep = "\t"
)

hcc_hepc4 <- colnames(hep)[
  hep$sub_celltype == "Hep_C4" & hep$group == "HCC"
]
if (length(hcc_hepc4) < 100L) {
  stop("Expected at least 100 HCC Hep_C4 cells; found ", length(hcc_hepc4))
}
hcc_expression <- GetAssayData(hep, assay = "RNA", layer = "data")[
  , hcc_hepc4, drop = FALSE
]
background <- rownames(hcc_expression)[
  Matrix::rowSums(hcc_expression > 0) >= 0.10 * length(hcc_hepc4)
]

pkn_genes <- c(
  "FGFR2", "TGFBR1", "TGFBR2", "FRS2", "SHC1", "GRB2", "SOS1",
  "HRAS", "NRAS", "KRAS", "RAF1", "BRAF", "ARAF", "MAP2K1", "MAP2K2",
  "MAPK1", "MAPK3", "MYC", "ELK1", "FOS", "JUN", "ETS1", "CREB1"
)
missing_pkn <- setdiff(pkn_genes, rownames(hcc_expression))
if (length(missing_pkn) > 0L) {
  stop(
    "Missing PKN genes in frozen-curated atlas: ",
    paste(missing_pkn, collapse = ", ")
  )
}
pkn_expression_matrix <- hcc_expression[pkn_genes, , drop = FALSE]
pkn_expression <- data.table(
  gene = pkn_genes,
  Hep_C4_HCC_avg_expr = as.numeric(Matrix::rowMeans(pkn_expression_matrix)),
  Hep_C4_HCC_pct_positive = as.numeric(
    Matrix::rowMeans(pkn_expression_matrix > 0) * 100
  ),
  n_Hep_C4_HCC_cells = length(hcc_hepc4)
)
fwrite(
  pkn_expression,
  file.path(output_root, "step5_HepC4_HCC_expression.csv")
)

msigdb_member <- paste0(
  "msigdb_v2023.2.Hs_files_to_download_locally/",
  "msigdb_v2023.2.Hs_GMTs/msigdb.v2023.2.Hs.symbols.gmt"
)
msigdb_connection <- unz(msigdb_zip_path, msigdb_member, open = "r")
msigdb_lines <- readLines(msigdb_connection, warn = FALSE)
close(msigdb_connection)
msigdb_parts <- strsplit(msigdb_lines, "\t", fixed = TRUE)
program_terms <- c(
  MAPK = "KEGG_MAPK_SIGNALING_PATHWAY",
  mTOR = "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  NFkB = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  Hypoxia = "HALLMARK_HYPOXIA"
)
target_rows <- list()
official_rows <- list()
for (program in names(program_terms)) {
  term <- program_terms[[program]]
  term_index <- which(vapply(
    msigdb_parts,
    function(values) length(values) >= 3L && identical(values[[1L]], term),
    logical(1L)
  ))
  if (length(term_index) != 1L) {
    stop("MSigDB v2023.2 term count is not one for ", term)
  }
  official <- sort(unique(toupper(msigdb_parts[[term_index]][-(1:2)])))
  if (length(official) == 0L) {
    stop("MSigDB v2023.2 lacks ", term)
  }
  expressed <- intersect(official, background)
  official_rows[[program]] <- data.table(
    program = program,
    gene = official,
    gene_set_name = term,
    gene_set_source = paste0(
      "data/reference/MSigDB/",
      basename(msigdb_zip_path),
      "@sha256:",
      expected_msigdb_sha256
    ),
    target_filter = "official_pathway_unfiltered",
    official_gene_set_n = length(official),
    hepc4_expressed_background_n = length(background),
    expressed_target_n = length(expressed)
  )
  target_rows[[program]] <- data.table(
    program = program,
    gene = expressed,
    gene_set_name = term,
    gene_set_source = paste0(
      "data/reference/MSigDB/",
      basename(msigdb_zip_path),
      "@sha256:",
      expected_msigdb_sha256
    ),
    target_filter = "official_pathway_intersect_hepc4_expressed_background_no_DE",
    official_gene_set_n = length(official),
    hepc4_expressed_background_n = length(background),
    expressed_target_n = length(expressed)
  )
}
fwrite(
  rbindlist(target_rows),
  file.path(output_root, "target_sets_4pathway_expressed_long.csv")
)
fwrite(
  rbindlist(official_rows),
  file.path(output_root, "target_sets_4pathway_official_long.csv")
)
fwrite(
  data.table(
    dependency = c(
      "frozen-curated atlas",
      "CellChat LR profiles",
      "HCC Hep_C4 expression",
      "four expressed target sets"
    ),
    status = "completed",
    n = c(
      ncol(atlas),
      length(genes),
      length(hcc_hepc4),
      nrow(rbindlist(target_rows))
    ),
    seed = seed
  ),
  file.path(output_root, "scrna_profile_status.tsv"),
  sep = "\t"
)
cat(
  "scRNA profiles exported from frozen-curated counts: ",
  length(genes),
  " LR genes; ",
  length(hcc_hepc4),
  " HCC Hep_C4 cells.\n",
  sep = ""
)
