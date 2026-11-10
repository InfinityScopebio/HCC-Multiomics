#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CoVarNet)
  library(data.table)
  library(dplyr)
  library(foreach)
  library(igraph)
  library(Matrix)
  library(NMF)
  library(psych)
  library(reshape2)
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
atlas_path <- assert_production_path(
  file.path(
    root,
    "data",
    "curation",
    "scRNA",
    "combined_annotated_v5.rds"
  )
)
curation_audit_path <- assert_production_path(
  file.path(root, "audit", "scrna_frozen_curation_objects.tsv")
)
upstream_model_path <- assert_production_path(
  file.path(
    root,
    "results",
    "20_scrna_cm02",
    "covarnet",
    "cm_network_frozen_curated.rds"
  )
)
vendor_path <- assert_production_path(
  file.path(
    root,
    "code",
    "pipelines",
    "20_scrna_cm02",
    "covarnet_v030_vendor.R"
  )
)
source(vendor_path)
covarnet_root <- file.path(result_root, "covarnet_reference")
dir.create(covarnet_root, recursive = TRUE, showWarnings = FALSE)
model_path <- file.path(
  covarnet_root,
  "cm_network_frozen_curated_reference_historical_filter.rds"
)
output_root <- file.path(result_root, "cm_recovery")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_audit_root(), recursive = TRUE, showWarnings = FALSE)

required_paths <- c(
  atlas_path,
  curation_audit_path,
  upstream_model_path,
  vendor_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing frozen-curated CoVarNet rebuild dependencies: ",
    paste(missing_paths, collapse = ", ")
  )
}
curation_audit <- fread(curation_audit_path)
validate_scrna_curation_audit(curation_audit, atlas_path)

map_major_cluster <- function(states) {
  hepatocyte <- paste0("Hep_C", 1:7)
  stromal <- c("HSC", "Pericyte", "SMC", "Endothelial", "LSEC")
  myeloid <- c(
    "Kupffer", "IL1B+M_", "TREM2+M_", "SPP1+M_", "CXCL10+M_",
    "Monocytes", "cDC1", "cDC2", "pDC", "Mast", "Neutrophil"
  )
  lymphocyte <- c(
    "CD4+Tn", "CD8+Tem", "CD8+Temra", "CD8+Trm", "MAIT", "gdT",
    "Treg", "NK", "naiveB", "memoryB", "Plasma"
  )
  ifelse(
    states %in% hepatocyte,
    "Hepatocyte",
    ifelse(
      states %in% stromal,
      "Stromal",
      ifelse(
        states %in% myeloid,
        "Myeloid",
        ifelse(states %in% lymphocyte, "Lymphocyte", NA_character_)
      )
    )
  )
}

jaccard <- function(left, right) {
  length(intersect(left, right)) / length(union(left, right))
}

run_covarnet_nmf <- function(frequency) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(77L)
  foreach::registerDoSEQ()
  NMF::nmf(
    frequency,
    rank = 5L,
    method = "nsNMF",
    seed = rep(77L, 6L),
    nrun = 30L,
    .options = "v",
    .pbackend = "seq"
  )
}

build_covarnet_reference <- function() {
  atlas <- readRDS(atlas_path)
  if (!inherits(atlas, "Seurat")) {
    stop("Frozen-curated scRNA reference is not a Seurat object.")
  }
  metadata <- atlas[[]]
  required_metadata <- c("sub_celltype", "sample", "group")
  missing_metadata <- setdiff(required_metadata, colnames(metadata))
  if (length(missing_metadata) > 0L) {
    stop(
      "Frozen-curated scRNA reference lacks metadata: ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  if (any(as.character(metadata$sample) %in% c("CC1", "CC2", "CC3"))) {
    stop("Frozen-curated scRNA reference includes prohibited CC1-CC3 samples.")
  }
  metadata$sub_celltype <- as.character(metadata$sub_celltype)
  metadata$sample <- as.character(metadata$sample)
  metadata$group <- as.character(metadata$group)
  metadata$majorCluster <- map_major_cluster(metadata$sub_celltype)
  metadata <- metadata[!is.na(metadata$majorCluster), , drop = FALSE]
  meta <- data.frame(
    sampleID = metadata$sample,
    majorCluster = metadata$majorCluster,
    subCluster = metadata$sub_celltype,
    site = factor(
      metadata$group,
      levels = c("Healthy", "Cirrhosis", "Adjacent", "HCC")
    ),
    stringsAsFactors = FALSE
  )

  raw_frequency <- vendor_freq_calculate(meta)
  normalized_frequency <- vendor_freq_normalize(
    raw_frequency,
    normalize = "minmax"
  )
  normalized_frequency[is.na(normalized_frequency)] <- 1
  keep_rows <- rowSums(normalized_frequency, na.rm = TRUE) > 0
  keep_columns <- colSums(normalized_frequency, na.rm = TRUE) > 0
  raw_frequency <- raw_frequency[
    keep_rows,
    c(TRUE, keep_columns),
    drop = FALSE
  ]
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
    protected_cm02_anchor := subCluster %in% cm02_covarnet_types()
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
    file.path(spatial_audit_root(), "spatial_covarnet_state_filter.tsv"),
    sep = "\t"
  )

  retained_states <- rownames(normalized_frequency)[!rownames(normalized_frequency) %in% state_statistics[dropped == TRUE, subCluster]]
  missing_anchors <- setdiff(cm02_covarnet_types(), retained_states)
  if (length(missing_anchors) > 0L) {
    stop(
      "Anchor-preserving CoVarNet filter lost canonical states: ",
      paste(missing_anchors, collapse = ", ")
    )
  }
  filtered_frequency <- normalized_frequency[
    retained_states,
    ,
    drop = FALSE
  ]
  filtered_raw <- raw_frequency[retained_states, , drop = FALSE]
  saveRDS(
    list(
      meta = meta,
      mat_raw = raw_frequency,
      mat_norm = normalized_frequency,
      mat_norm_filtered = filtered_frequency,
      retained_states = retained_states
    ),
    file.path(
      covarnet_root,
      "frequency_frozen_curated_reference_historical_filter.rds"
    )
  )

  first <- run_covarnet_nmf(filtered_frequency)
  second <- run_covarnet_nmf(filtered_frequency)
  basis_delta <- max(abs(NMF::basis(first) - NMF::basis(second)))
  coefficient_delta <- max(abs(NMF::coef(first) - NMF::coef(second)))
  nmf_exact <- identical(NMF::basis(first), NMF::basis(second)) &&
    identical(NMF::coef(first), NMF::coef(second))
  nmf_repeatability <- data.table(
    component = c("basis", "coefficient"),
    seed = 77L,
    run1 = "in_process_run1",
    run2 = "in_process_run2",
    maximum_numeric_delta = c(basis_delta, coefficient_delta),
    status = if (nmf_exact) "exact" else "non_deterministic",
    rng_kind = "Mersenne-Twister; Inversion; Rejection",
    backend = "foreach::doSEQ; NMF::nsNMF nrun=30"
  )
  fwrite(
    nmf_repeatability,
    file.path(spatial_audit_root(), "spatial_covarnet_nmf_repeatability.tsv"),
    sep = "\t"
  )
  if (!nmf_exact) {
    stop(
      "Two fixed-seed CoVarNet NMF runs were not bitwise exact; basis delta=",
      basis_delta,
      ", coefficient delta=",
      coefficient_delta
    )
  }

  final_module_names <- sprintf("CM%02d", seq_len(5L))
  colnames(NMF::basis(first)) <- final_module_names
  rownames(NMF::coef(first)) <- final_module_names

  normalized_basis <- sweep(
    NMF::basis(first),
    2L,
    colSums(NMF::basis(first)),
    "/"
  )
  top_count <- min(10L, nrow(normalized_basis))
  top_states <- lapply(seq_len(ncol(normalized_basis)), function(index) {
    names(sort(normalized_basis[, index], decreasing = TRUE))[seq_len(top_count)]
  })
  anchor_overlap <- vapply(
    top_states,
    function(states) jaccard(states, cm02_covarnet_types()),
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

  pairwise_correlation <- vendor_pair_correlation(
    filtered_raw,
    method = "pearson"
  )
  network <- vendor_cm_network(
    first,
    pairwise_correlation,
    top_n = 10L,
    corr = -2,
    pval_fdr = 0.05
  )
  upstream_network <- readRDS(upstream_model_path)
  upstream_basis <- NMF::basis(upstream_network$raw)
  upstream_coefficient <- NMF::coef(upstream_network$raw)
  basis_compatible <- identical(dimnames(upstream_basis), dimnames(NMF::basis(first)))
  coefficient_compatible <- identical(
    dimnames(upstream_coefficient),
    dimnames(NMF::coef(first))
  )
  upstream_basis_delta <- if (basis_compatible) {
    max(abs(upstream_basis - NMF::basis(first)))
  } else {
    Inf
  }
  upstream_coefficient_delta <- if (coefficient_compatible) {
    max(abs(upstream_coefficient - NMF::coef(first)))
  } else {
    Inf
  }
  filter_keys <- function(value) {
    sort(paste(value$cm, value$subCluster, sep = "\r"))
  }
  filter_exact <- identical(
    filter_keys(upstream_network$filter),
    filter_keys(network$filter)
  )
  upstream_equivalence <- data.table(
    check = c("NMF_basis", "NMF_coefficient", "network_filter"),
    maximum_numeric_delta = c(
      upstream_basis_delta,
      upstream_coefficient_delta,
      NA_real_
    ),
    exact = c(
      isTRUE(upstream_basis_delta == 0),
      isTRUE(upstream_coefficient_delta == 0),
      filter_exact
    ),
    status = ifelse(
      c(
        isTRUE(upstream_basis_delta == 0),
        isTRUE(upstream_coefficient_delta == 0),
        filter_exact
      ),
      "PASS",
      "FAIL"
    )
  )
  fwrite(
    upstream_equivalence,
    file.path(
      spatial_audit_root(),
      "spatial_covarnet_upstream_equivalence.tsv"
    ),
    sep = "\t"
  )
  if (any(upstream_equivalence$status != "PASS")) {
    stop(
      "Spatial CoVarNet reference differs from the Figure 3 production model."
    )
  }
  observed_cm02 <- network$filter[
    network$filter$cm == "CM02",
    "subCluster"
  ]
  if (!setequal(observed_cm02, cm02_covarnet_types())) {
    stop(
      "Anchor-preserving frozen-curated CM02 is not canonical. Missing: ",
      paste(setdiff(cm02_covarnet_types(), observed_cm02), collapse = ", "),
      "; unexpected: ",
      paste(setdiff(observed_cm02, cm02_covarnet_types()), collapse = ", ")
    )
  }
  saveRDS(network, model_path)
  fwrite(
    network$filter,
    file.path(covarnet_root, "cm_filter_anchor_preserved.tsv"),
    sep = "\t"
  )
  cm02_basis <- sweep(
    NMF::basis(first),
    2L,
    colSums(NMF::basis(first)),
    "/"
  )[, "CM02"]
  scope <- data.table(
    state = names(cm02_basis),
    normalized_nmf_loading = as.numeric(cm02_basis),
    nmf_loading_rank = rank(-cm02_basis, ties.method = "first"),
    in_network = names(cm02_basis) %in% observed_cm02,
    predeclared_anchor = names(cm02_basis) %in% cm02_covarnet_types()
  )
  setorder(scope, nmf_loading_rank)
  fwrite(
    scope,
    file.path(spatial_audit_root(), "spatial_cm02_scope_rebuild.tsv"),
    sep = "\t"
  )
  network
}

reference <- build_covarnet_reference()
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
if (!is.list(reference) || is.null(reference$raw)) {
  stop("CoVarNet model lacks its raw NMF factor object.")
}
factor_states <- rownames(NMF::basis(reference$raw))
if (length(factor_states) == 0L) {
  stop("CoVarNet NMF basis has no state names.")
}
factor_indices <- match(factor_states, rownames(reference$ann))
if (anyNA(factor_indices)) {
  stop(
    "CoVarNet annotation lacks NMF basis states: ",
    paste(factor_states[is.na(factor_indices)], collapse = ", ")
  )
}
reference$ann <- reference$ann[
  factor_indices,
  ,
  drop = FALSE
]
rownames(reference$ann) <- factor_states

recover_section <- function(section) {
  weight_path <- active_rctd_weight_path(section)
  if (!file.exists(weight_path)) {
    stop("Missing production RCTD weights: ", weight_path)
  }
  weights <- readRDS(weight_path)
  colnames(weights) <- restore_covarnet_labels(colnames(weights))
  missing_states <- setdiff(factor_states, colnames(weights))
  if (length(missing_states) > 0L) {
    stop(
      section,
      ": RCTD weights lack CoVarNet factor states: ",
      paste(missing_states, collapse = ", ")
    )
  }
  proportions <- weights[, factor_states, drop = FALSE]
  rownames(proportions) <- paste(section, weights$barcode, sep = "::")
  dummy <- Matrix(
    1,
    nrow = 2L,
    ncol = nrow(proportions),
    sparse = TRUE,
    dimnames = list(c("dummy_1", "dummy_2"), rownames(proportions))
  )
  spatial <- CreateSeuratObject(
    counts = dummy,
    meta.data = proportions
  )
  normalized <- freq_normalize_st(
    spatial,
    reference,
    normalize = "minmax"
  )
  recovered <- sc_cm_recover(reference, normalized)
  observed_cm_names <- rownames(recovered)
  canonical_cm_order <- sprintf(
    "CM%02d",
    seq_len(length(observed_cm_names))
  )
  if (!setequal(observed_cm_names, canonical_cm_order)) {
    stop(
      section,
      ": recovered CM names do not match canonical order ",
      paste(canonical_cm_order, collapse = ", "),
      "; observed ",
      paste(observed_cm_names, collapse = ", ")
    )
  }
  expected_spot_ids <- paste(section, weights$barcode, sep = "::")
  observed_spot_ids <- colnames(recovered)
  if (anyDuplicated(expected_spot_ids) || anyDuplicated(observed_spot_ids)) {
    stop(section, ": duplicated spot IDs prevent deterministic alignment.")
  }
  if (!setequal(observed_spot_ids, expected_spot_ids)) {
    stop(
      section,
      ": recovered spot IDs differ from production RCTD barcodes."
    )
  }
  recovered <- recovered[
    canonical_cm_order,
    expected_spot_ids,
    drop = FALSE
  ]
  recovered_table <- as.data.frame(t(recovered), check.names = FALSE)
  colnames(recovered_table) <- canonical_cm_order
  score_matrix <- as.matrix(
    recovered_table[, canonical_cm_order, drop = FALSE]
  )
  storage.mode(score_matrix) <- "double"
  if (any(!is.finite(score_matrix))) {
    stop(section, ": recovered CM scores contain non-finite values.")
  }
  first_index <- max.col(score_matrix, ties.method = "first")
  last_index <- max.col(score_matrix, ties.method = "last")
  first_cm <- canonical_cm_order[first_index]
  last_cm <- canonical_cm_order[last_index]
  row_maximum <- apply(score_matrix, 1L, max)
  exact_tie_membership <- score_matrix ==
    matrix(
      row_maximum,
      nrow = nrow(score_matrix),
      ncol = ncol(score_matrix)
    )
  exact_tie_count <- rowSums(exact_tie_membership)
  sorted_scores <- t(apply(score_matrix, 1L, sort, decreasing = TRUE))
  top_margin <- sorted_scores[, 1L] - sorted_scores[, 2L]
  spot_id <- expected_spot_ids

  tied_index <- which(exact_tie_count > 1L)
  tied_cms <- character(length(tied_index))
  if (length(tied_index) > 0L) {
    tied_cms <- apply(
      exact_tie_membership[tied_index, , drop = FALSE],
      1L,
      function(values) {
        paste(canonical_cm_order[values], collapse = "|")
      }
    )
  }
  tie_spots <- data.frame(
    section = rep(section, length(tied_index)),
    barcode = weights$barcode[tied_index],
    spot_id = spot_id[tied_index],
    tie_count = exact_tie_count[tied_index],
    tied_cms = tied_cms,
    canonical_first = first_cm[tied_index],
    canonical_last = last_cm[tied_index],
    cm02_status_changes = (
      first_cm[tied_index] == "CM02"
    ) != (
      last_cm[tied_index] == "CM02"
    ),
    top_margin = top_margin[tied_index],
    stringsAsFactors = FALSE
  )

  tolerances <- c(0, 1e-12, 1e-10, 1e-8, 1e-6)
  score_gap <- matrix(
    row_maximum,
    nrow = nrow(score_matrix),
    ncol = ncol(score_matrix)
  ) - score_matrix
  cm02_index <- match("CM02", canonical_cm_order)
  if (is.na(cm02_index)) {
    stop(section, ": canonical CM order does not contain CM02.")
  }
  sensitivity <- bind_rows(lapply(tolerances, function(tolerance) {
    candidate <- score_gap <= tolerance
    ambiguous <- rowSums(candidate) > 1L
    data.frame(
      section = section,
      tolerance = tolerance,
      production_tie_policy = "canonical_first",
      sensitivity_tie_policy = "canonical_last",
      canonical_cm_order = paste(canonical_cm_order, collapse = "|"),
      n_spots = nrow(score_matrix),
      n_ambiguous = sum(ambiguous),
      n_ambiguous_involving_cm02 = sum(
        ambiguous & candidate[, cm02_index]
      ),
      n_exact_ties = sum(exact_tie_count > 1L),
      n_first_last_changed = sum(first_cm != last_cm),
      n_cm02_status_changed = sum(
        (first_cm == "CM02") != (last_cm == "CM02")
      ),
      n_cm02_positive_first = sum(first_cm == "CM02"),
      n_cm02_positive_last = sum(last_cm == "CM02"),
      stringsAsFactors = FALSE
    )
  }))

  recovered_table$section <- section
  recovered_table$patient <- sub("[A-Z]$", "", section)
  recovered_table$region <- substring(section, nchar(section))
  recovered_table$barcode <- weights$barcode
  recovered_table$row <- weights$row
  recovered_table$col <- weights$col
  recovered_table$maxCM <- first_cm
  list(
    recovered = recovered_table,
    tie_spots = tie_spots,
    sensitivity = sensitivity
  )
}

set.seed(seed)
section_results <- lapply(all_sections(), recover_section)
recovered <- bind_rows(lapply(section_results, `[[`, "recovered"))
tie_spots <- bind_rows(lapply(section_results, `[[`, "tie_spots"))
sensitivity <- bind_rows(lapply(section_results, `[[`, "sensitivity"))
overall_sensitivity <- sensitivity %>%
  group_by(tolerance) %>%
  summarise(
    section = "ALL",
    production_tie_policy = "canonical_first",
    sensitivity_tie_policy = "canonical_last",
    canonical_cm_order = paste(sprintf("CM%02d", 1:5), collapse = "|"),
    n_spots = sum(n_spots),
    n_ambiguous = sum(n_ambiguous),
    n_ambiguous_involving_cm02 = sum(n_ambiguous_involving_cm02),
    n_exact_ties = sum(n_exact_ties),
    n_first_last_changed = sum(n_first_last_changed),
    n_cm02_status_changed = sum(n_cm02_status_changed),
    n_cm02_positive_first = sum(n_cm02_positive_first),
    n_cm02_positive_last = sum(n_cm02_positive_last),
    .groups = "drop"
  ) %>%
  select(
    section,
    tolerance,
    everything()
  )
sensitivity <- bind_rows(sensitivity, overall_sensitivity)
cm_columns <- grep("^CM[0-9]+$", colnames(recovered), value = TRUE)
if (!"CM02" %in% cm_columns) {
  stop("Recovered output does not contain CM02.")
}
audit_root <- spatial_audit_root()
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)
fwrite(
  as.data.table(tie_spots),
  file.path(audit_root, "spatial_cm_max_tie_spots.tsv"),
  sep = "\t"
)
fwrite(
  as.data.table(sensitivity),
  file.path(audit_root, "spatial_cm_max_tie_sensitivity.tsv"),
  sep = "\t"
)
saveRDS(recovered, file.path(output_root, "CM_recover.rds"))
saveRDS(recovered, file.path(output_root, "CM_recover_with_barcode.rds"))
fwrite(
  as.data.table(recovered),
  file.path(output_root, "CM_recover.csv.gz")
)
fwrite(
  as.data.table(recovered),
  file.path(output_root, "CM_recover_with_barcode.csv.gz")
)

fractions <- as.data.table(recovered)[
  ,
  .(
    n_spots = .N,
    n_cm02_positive = sum(maxCM == "CM02"),
    cm02_fraction = mean(maxCM == "CM02")
  ),
  by = .(section, patient, region)
]
fractions[, is_tumor_scope := section %in% tumor_sections()]
fwrite(fractions, file.path(output_root, "CM02_fraction_by_section.tsv"), sep = "\t")
fwrite(
  data.table(
    run_variant = spatial_run_variant(),
    rctd_source = rctd_source(),
    n_spots = nrow(recovered),
    n_cm02_positive = sum(recovered$maxCM == "CM02"),
    article_equivalent = identical(rctd_source(), "historical_frozen")
  ),
  file.path(output_root, "rctd_recovery_summary.tsv"),
  sep = "\t"
)
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
cat(
  "CM recovery completed: ",
  nrow(recovered),
  " spots; ",
  sum(recovered$maxCM == "CM02"),
  " maxCM == CM02; ",
  nrow(tie_spots),
  " exact max-score ties resolved by canonical-first policy.\n",
  sep = ""
)
