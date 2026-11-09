#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(foreach)
  library(igraph)
  library(NMF)
  library(psych)
  library(reshape2)
  library(Seurat)
})

arguments <- commandArgs(trailingOnly = FALSE)
this_script <- normalizePath(
  sub("^--file=", "", arguments[grepl("^--file=", arguments)]),
  mustWork = TRUE
)
source(file.path(dirname(this_script), "config.R"))
source(file.path(dirname(this_script), "covarnet_v030_vendor.R"))
config <- ensure_output_dirs(scrna_cm02_config(this_script))
setwd(config$project_root)
assert_production_input(config, config$curated_atlas_paths$final)
validate_frozen_curation_audit(config)

output_dir <- file.path(config$result_root, "covarnet")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
atlas <- readRDS(config$curated_atlas_paths$final)
metadata <- atlas[[]]
metadata$sub_celltype <- as.character(metadata$sub_celltype)
metadata$sample <- as.character(metadata$sample)
metadata$group <- as.character(metadata$group)
assert_no_cc(metadata$sample, "frozen-curated CoVarNet atlas")

hepatocyte_states <- paste0("Hep_C", 1:7)
stromal_states <- c("HSC", "Pericyte", "SMC", "Endothelial", "LSEC")
myeloid_states <- c(
  "Kupffer", "IL1B+M_", "TREM2+M_", "SPP1+M_", "CXCL10+M_",
  "Monocytes", "cDC1", "cDC2", "pDC", "Mast", "Neutrophil"
)
lymphocyte_states <- c(
  "CD4+Tn", "CD8+Tem", "CD8+Temra", "CD8+Trm", "MAIT", "gdT",
  "Treg", "NK", "naiveB", "memoryB", "Plasma"
)
map_major <- function(states) {
  ifelse(
    states %in% hepatocyte_states,
    "Hepatocyte",
    ifelse(
      states %in% stromal_states,
      "Stromal",
      ifelse(
        states %in% myeloid_states,
        "Myeloid",
        ifelse(states %in% lymphocyte_states, "Lymphocyte", NA_character_)
      )
    )
  )
}

metadata$majorCluster <- map_major(metadata$sub_celltype)
metadata <- metadata[!is.na(metadata$majorCluster), , drop = FALSE]
meta <- data.frame(
  sampleID = metadata$sample,
  majorCluster = metadata$majorCluster,
  subCluster = metadata$sub_celltype,
  site = factor(metadata$group, levels = config$group_levels),
  stringsAsFactors = FALSE
)
assert_no_cc(meta$sampleID, "frozen-curated CoVarNet frequency metadata")

raw_frequency <- vendor_freq_calculate(meta)
normalized_frequency <- vendor_freq_normalize(
  raw_frequency,
  normalize = "minmax"
)
normalized_frequency[is.na(normalized_frequency)] <- 1
keep_rows <- rowSums(normalized_frequency, na.rm = TRUE) > 0
keep_columns <- colSums(normalized_frequency, na.rm = TRUE) > 0
raw_frequency <- raw_frequency[keep_rows, c(TRUE, keep_columns), drop = FALSE]
normalized_frequency <- normalized_frequency[
  keep_rows,
  keep_columns,
  drop = FALSE
]

state_statistics <- data.table(
  subCluster = meta$subCluster
)[
  ,
  .(n_cells = .N),
  by = subCluster
]
state_statistics[
  ,
  pct_total := 100 * n_cells / sum(n_cells)
]
sample_presence <- data.table(
  subCluster = meta$subCluster,
  sampleID = meta$sampleID
)[
  ,
  .(n_samples = uniqueN(sampleID)),
  by = subCluster
]
state_statistics <- merge(
  state_statistics,
  sample_presence,
  by = "subCluster"
)
state_statistics[
  ,
  var_norm := apply(normalized_frequency, 1L, var)[subCluster]
]
raw_frequency_numeric <- as.matrix(
  raw_frequency[
    ,
    setdiff(colnames(raw_frequency), "majorCluster"),
    drop = FALSE
  ]
)
state_statistics[
  ,
  cv_in_major := (
    apply(raw_frequency_numeric, 1L, sd) /
      pmax(rowMeans(raw_frequency_numeric), 1e-9)
  )[subCluster]
]
state_statistics[
  ,
  protected_cm02_anchor := subCluster %in% config$cm02_anchor_states
]
state_statistics[
  ,
  drop_low_abundance := pct_total < 0.5 & n_samples < 20
]
state_statistics[
  ,
  drop_low_cv := cv_in_major < 1
]
state_statistics[
  ,
  dropped := drop_low_abundance | drop_low_cv
]
state_statistics[
  ,
  rule := fifelse(
    drop_low_abundance,
    "pct_total<0.5_and_n_samples<20",
    fifelse(
      drop_low_cv,
      "cv_in_major<1",
      "retained"
    )
  )
]
fwrite(
  state_statistics,
  file.path(config$audit_root, "scrna_covarnet_state_filter.tsv"),
  sep = "\t"
)
dropped_states <- state_statistics[dropped == TRUE, subCluster]
filtered_frequency <- normalized_frequency[
  setdiff(rownames(normalized_frequency), dropped_states),
  ,
  drop = FALSE
]
filtered_raw <- raw_frequency[
  rownames(filtered_frequency),
  ,
  drop = FALSE
]
saveRDS(
  list(
    meta = meta,
    mat_raw = raw_frequency,
    mat_norm = normalized_frequency,
    mat_norm_filtered = filtered_frequency,
    dropped_states = dropped_states
  ),
  config$covarnet_paths$frequency
)
fwrite(
  as.data.table(
    filtered_frequency,
    keep.rownames = "subCluster"
  ),
  file.path(output_dir, "frequency_filtered_frozen_curated.tsv"),
  sep = "\t"
)

set_recorded_seed(config, "covarnet", "foreach_backend")
foreach::registerDoSEQ()
set_recorded_seed(config, "covarnet", "nsNMF_K5")
nmf_result <- NMF::nmf(
  filtered_frequency,
  rank = 5L,
  method = "nsNMF",
  seed = rep(77L, 6L),
  nrun = 30L,
  .options = "v",
  .pbackend = "seq"
)
final_module_names <- sprintf("CM%02d", seq_len(5L))
colnames(NMF::basis(nmf_result)) <- final_module_names
rownames(NMF::coef(nmf_result)) <- final_module_names

normalized_basis <- sweep(
  NMF::basis(nmf_result),
  2L,
  colSums(NMF::basis(nmf_result)),
  "/"
)
jaccard <- function(left, right) {
  length(intersect(left, right)) / length(union(left, right))
}
top_count <- min(10L, nrow(normalized_basis))
top_states <- lapply(seq_len(ncol(normalized_basis)), function(index) {
  names(sort(normalized_basis[, index], decreasing = TRUE))[seq_len(top_count)]
})
anchor_overlap <- vapply(
  top_states,
  function(states) jaccard(states, config$cm02_anchor_states),
  numeric(1)
)
cm02_raw_index <- which.max(anchor_overlap)
if (cm02_raw_index != 2L) {
  stop(
    "Historical component-order contract failed: CM02 anchor best matched ",
    final_module_names[[cm02_raw_index]],
    " instead of CM02."
  )
}
saveRDS(nmf_result, config$covarnet_paths$nmf)

module_label_audit <- data.table(
  raw_module = final_module_names,
  raw_index = seq_len(5L),
  top_states = vapply(top_states, paste, collapse = ";", character(1)),
  anchor_jaccard = anchor_overlap,
  assigned_module = final_module_names,
  semantic_rule = ifelse(
    seq_len(5L) == cm02_raw_index,
    "maximum_top10_Jaccard_to_predeclared_CM02_anchor_states",
    "historical_NMF_component_order_preserved"
  )
)
fwrite(
  module_label_audit,
  file.path(config$audit_root, "scrna_covarnet_module_labeling.tsv"),
  sep = "\t"
)

pairwise_correlation <- vendor_pair_correlation(
  filtered_raw,
  method = "pearson"
)
fwrite(
  pairwise_correlation,
  file.path(output_dir, "pair_correlation_frozen_curated.tsv"),
  sep = "\t"
)
network <- vendor_cm_network(
  nmf_result,
  pairwise_correlation,
  top_n = 10L,
  corr = -2,
  pval_fdr = 0.05
)
saveRDS(
  network,
  config$covarnet_paths$network
)
fwrite(
  network$filter,
  file.path(output_dir, "cm_filter_top10_frozen_curated.tsv"),
  sep = "\t"
)
fwrite(
  network$each$edge,
  file.path(output_dir, "cm_each_edges_frozen_curated.tsv"),
  sep = "\t"
)

cm02_basis <- sweep(
  NMF::basis(nmf_result),
  2L,
  colSums(NMF::basis(nmf_result)),
  "/"
)[, "CM02"]
cm02_network_states <- network$filter$subCluster[
  network$filter$cm == "CM02"
]
cm02_scope <- data.table(
  state = names(cm02_basis),
  normalized_nmf_loading = as.numeric(cm02_basis),
  nmf_loading_rank = rank(-cm02_basis, ties.method = "first"),
  in_nmf_top10 = names(cm02_basis) %in%
    names(sort(cm02_basis, decreasing = TRUE))[seq_len(top_count)],
  in_network = names(cm02_basis) %in% cm02_network_states,
  predeclared_anchor = names(cm02_basis) %in% config$cm02_anchor_states
)
setorder(cm02_scope, nmf_loading_rank)
fwrite(
  cm02_scope,
  file.path(config$audit_root, "scrna_cm02_scope.tsv"),
  sep = "\t"
)

coefficient <- NMF::coef(nmf_result)
membership <- data.table(
  sampleID = colnames(coefficient),
  CMT = sub(
    "^CM",
    "CMT",
    apply(
      coefficient,
      2L,
      function(values) rownames(coefficient)[which.max(values)]
    )
  )
)
fwrite(
  membership,
  file.path(output_dir, "cm_membership_frozen_curated.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    component = "CoVarNet vendor",
    upstream_version = "v0.3.0",
    upstream_commit = "d4c52a4b2ea821dc4b75861c3772e21232d1ee51",
    vendor_file = "code/pipelines/20_scrna_cm02/covarnet_v030_vendor.R",
    historical_fdr_behavior = paste(
      "Historical fdr=1 argument was ignored by upstream; reproduced as",
      "explicit pval_fdr=0.05"
    )
  ),
  file.path(config$audit_root, "scrna_covarnet_vendor.tsv"),
  sep = "\t"
)
message(
  "PASS: frozen expert-curated frequencies -> strict filtering -> pinned ",
  "K=5 nsNMF -> vendored CoVarNet network completed."
)
