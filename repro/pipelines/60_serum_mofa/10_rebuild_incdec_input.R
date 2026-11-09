#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(impute)
})

script_dir <- dirname(
  normalizePath(
    sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
      grepl("^--file=", commandArgs(trailingOnly = FALSE))
    ]),
    mustWork = TRUE
  )
)
source(file.path(script_dir, "config.R"))
cfg <- load_pipeline_config()
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(cfg$master_seed)

paths <- c(
  protein = file.path(cfg$m0_dir, "serum_protein.rds"),
  phospho = file.path(cfg$m0_dir, "serum_phospho.rds"),
  lcms = file.path(cfg$m0_dir, "serum_lcms.rds"),
  gcms = file.path(cfg$m0_dir, "serum_gcms.rds"),
  samples = file.path(cfg$m0_dir, "sample_sheet.csv"),
  validation_target = cfg$validation$input
)
require_local_files(paths)

prepare_matrix <- function(mat, ids, minimum_detection = 0.5) {
  transformed <- log2(mat + 1)
  valid_ids <- !is.na(ids) & ids != ""
  transformed <- transformed[valid_ids, , drop = FALSE]
  ids <- ids[valid_ids]
  detected <- rowMeans(!is.na(transformed)) >= minimum_detection
  transformed <- transformed[detected, , drop = FALSE]
  ids <- ids[detected]
  sample_medians <- apply(transformed, 2, median, na.rm = TRUE)
  transformed <- sweep(
    transformed,
    2,
    sample_medians - median(sample_medians, na.rm = TRUE),
    "-"
  )
  transformed <- impute.knn(as.matrix(transformed), k = 10)$data
  list(mat = transformed, ids = ids)
}

clean_feature_names <- function(mat) {
  feature_names <- iconv(
    rownames(mat),
    from = "UTF-8",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  feature_names <- gsub("[^A-Za-z0-9_().:/+-]", "_", feature_names)
  rownames(mat) <- make.unique(gsub("_+", "_", feature_names))
  mat
}

strict_monotonic <- function(mat, group) {
  apply(mat, 1, function(values) {
    group_means <- tapply(values, group, mean)[
      c("Healthy", "HepFibrosis", "Cirrhosis", "HCC")
    ]
    if (any(is.na(group_means)) || var(values, na.rm = TRUE) < 1e-10) {
      return(FALSE)
    }
    all(diff(group_means) > 0) || all(diff(group_means) < 0)
  })
}

serum_protein <- readRDS(paths[["protein"]])
serum_phospho <- readRDS(paths[["phospho"]])
serum_lcms <- readRDS(paths[["lcms"]])
serum_gcms <- readRDS(paths[["gcms"]])
sample_sheet <- fread(paths[["samples"]])
serum_samples <- sample_sheet[source == "serum"]

protein_prepared <- prepare_matrix(
  serum_protein$mat,
  serum_protein$meta$SYMBOL
)
protein <- protein_prepared$mat
protein_symbols <- protein_prepared$ids
protein <- rowsum(protein, protein_symbols) /
  as.numeric(table(protein_symbols)[rownames(rowsum(protein, protein_symbols))])

phospho_prepared <- prepare_matrix(
  serum_phospho$mat,
  serum_phospho$meta$SYMBOL
)
phospho <- phospho_prepared$mat
rownames(phospho) <- make.unique(phospho_prepared$ids)

lcms <- log2(serum_lcms$mat + 1)
lcms_names <- serum_lcms$meta$Metabolites
lcms_keep <- !is.na(lcms_names) & lcms_names != "" &
  rowMeans(!is.na(lcms)) >= 0.5
lcms <- lcms[lcms_keep, , drop = FALSE]
lcms_names <- lcms_names[lcms_keep]
lcms_medians <- apply(lcms, 2, median, na.rm = TRUE)
lcms <- sweep(lcms, 2, lcms_medians - median(lcms_medians, na.rm = TRUE), "-")
lcms <- impute.knn(as.matrix(lcms), k = 10)$data
rownames(lcms) <- make.unique(as.character(lcms_names))

gcms_prepared <- prepare_matrix(
  serum_gcms$mat,
  serum_gcms$meta$Metabolites
)
gcms <- gcms_prepared$mat
rownames(gcms) <- make.unique(as.character(gcms_prepared$ids))

protein <- clean_feature_names(protein)
phospho <- clean_feature_names(phospho)
lcms <- clean_feature_names(lcms)
gcms <- clean_feature_names(gcms)

common <- Reduce(
  intersect,
  list(
    colnames(protein),
    colnames(phospho),
    colnames(lcms),
    colnames(gcms)
  )
)
common <- common[common %in% serum_samples$sample_id]
group <- factor(
  serum_samples$group_4[match(common, serum_samples$sample_id)],
  levels = c("Healthy", "HepFibrosis", "Cirrhosis", "HCC")
)
stage <- c(Healthy = 0, HepFibrosis = 1, Cirrhosis = 2, HCC = 3)[
  as.character(group)
]

protein <- protein[, common, drop = FALSE]
phospho <- phospho[, common, drop = FALSE]
lcms <- lcms[, common, drop = FALSE]
gcms <- gcms[, common, drop = FALSE]

views <- list(
  Protein = protein[strict_monotonic(protein, group), , drop = FALSE],
  Phospho = phospho[strict_monotonic(phospho, group), , drop = FALSE],
  LCMS = lcms[strict_monotonic(lcms, group), , drop = FALSE],
  GCMS = gcms[strict_monotonic(gcms, group), , drop = FALSE]
)
views <- views[vapply(views, nrow, integer(1)) >= 3L]
rebuilt <- list(views = views, common = common, grp = group, stage = stage)

rebuilt_path <- file.path(cfg$result_dir, "incdec_input_rebuilt.rds")
saveRDS(rebuilt, rebuilt_path)
historical <- readRDS(paths[["validation_target"]])

comparison_rows <- lapply(names(historical$views), function(view_name) {
  expected <- historical$views[[view_name]]
  observed <- rebuilt$views[[view_name]]
  dimensions_identical <- identical(dim(expected), dim(observed))
  row_order_identical <- identical(rownames(expected), rownames(observed))
  column_order_identical <- identical(colnames(expected), colnames(observed))
  values_identical <- identical(expected, observed)
  max_abs_diff <- if (
    dimensions_identical && row_order_identical && column_order_identical
  ) {
    max(abs(expected - observed), na.rm = TRUE)
  } else {
    NA_real_
  }
  data.table(
    component = view_name,
    expected_rows = nrow(expected),
    rebuilt_rows = nrow(observed),
    expected_columns = ncol(expected),
    rebuilt_columns = ncol(observed),
    dimensions_identical,
    row_order_identical,
    column_order_identical,
    values_identical,
    max_abs_diff
  )
})
comparison <- rbindlist(comparison_rows)
metadata_identical <- identical(
  historical[c("common", "grp", "stage")],
  rebuilt[c("common", "grp", "stage")]
)
comparison <- rbind(
  comparison,
  data.table(
    component = "sample_metadata",
    expected_rows = length(historical$common),
    rebuilt_rows = length(rebuilt$common),
    expected_columns = NA_integer_,
    rebuilt_columns = NA_integer_,
    dimensions_identical = length(historical$common) == length(rebuilt$common),
    row_order_identical = identical(historical$common, rebuilt$common),
    column_order_identical = NA,
    values_identical = metadata_identical,
    max_abs_diff = if (metadata_identical) 0 else NA_real_
  )
)
comparison[, status := ifelse(values_identical, "exact", "mismatch")]

fwrite(
  comparison,
  file.path(cfg$result_dir, "incdec_input_comparison.tsv"),
  sep = "\t"
)
fwrite(
  comparison,
  file.path(cfg$audit_dir, "mofa_incdec_input_comparison.tsv"),
  sep = "\t"
)

if (any(comparison$status != "exact")) {
  stop("Rebuilt strict-monotonic MOFA input differs from the frozen input.")
}
cat("Rebuilt incdec_input is numerically exact for all four views and metadata.\n")
