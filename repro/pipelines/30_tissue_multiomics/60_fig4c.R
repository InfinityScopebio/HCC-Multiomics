#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_path <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

protein <- readRDS(data_path("tissue_protein.rds"))
metabolite <- readRDS(data_path("tissue_metabolite.rds"))
sample_sheet <- read.csv(data_path("sample_sheet.csv"))
tissue <- sample_sheet[sample_sheet$source == "tissue", ]

prepare_matrix <- function(matrix) {
  transformed <- log2(matrix + 1)
  medians <- apply(transformed, 2, median, na.rm = TRUE)
  transformed <- sweep(
    transformed,
    2,
    medians - median(medians, na.rm = TRUE),
    "-"
  )
  transformed[is.na(transformed)] <- min(
    transformed, na.rm = TRUE
  ) / 2
  transformed
}

protein_matrix <- prepare_matrix(protein$mat)
metabolite_matrix <- prepare_matrix(metabolite$mat)
protein_genes <- protein$meta$Gene
metabolite_names <- metabolite$meta$Metabolites
metabolite_platform <- as.character(metabolite$meta$`Ion mode`)

common_samples <- intersect(
  colnames(protein_matrix), colnames(metabolite_matrix)
)
common_samples <- common_samples[
  tissue$group_4[
    match(common_samples, tissue$sample_id)
  ] %in% c("HCC", "Adjacent")
]
group <- tissue$group_4[match(common_samples, tissue$sample_id)]
pair <- tissue$paired_id[match(common_samples, tissue$sample_id)]
sample_order <- order(group != "Adjacent", pair)
common_samples <- common_samples[sample_order]
group <- group[sample_order]
pair <- pair[sample_order]
protein_matrix <- protein_matrix[, common_samples, drop = FALSE]
metabolite_matrix <- metabolite_matrix[
  , common_samples, drop = FALSE
]
hcc_index <- which(group == "HCC")
adjacent_index <- which(group == "Adjacent")
stopifnot(
  identical(
    as.character(pair[hcc_index]),
    as.character(pair[adjacent_index])
  )
)

msigdb <- as.data.table(readRDS(result_path(
  "references", "msigdb_tissue_sets.rds"
)))
gene_set <- function(name) {
  unique(toupper(msigdb$gene_symbol[msigdb$gs_name == name]))
}
pathways <- list(
  "NF-kB (TNFA)" = gene_set("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  "MAPK-ERK" = gene_set("REACTOME_MAPK1_MAPK3_SIGNALING"),
  "mTOR (PI3K-AKT)" = gene_set("HALLMARK_PI3K_AKT_MTOR_SIGNALING"),
  "Hypoxia" = gene_set("HALLMARK_HYPOXIA")
)

protein_assignment <- rep(NA_character_, nrow(protein_matrix))
protein_upper <- toupper(protein_genes)
for (pathway in names(pathways)) {
  hit <- (
    protein_upper %in% pathways[[pathway]] &
      is.na(protein_assignment)
  )
  protein_assignment[hit] <- pathway
}
protein_delta <- (
  rowMeans(protein_matrix[, hcc_index]) -
    rowMeans(protein_matrix[, adjacent_index])
)
protein_p <- vapply(seq_len(nrow(protein_matrix)), function(index) {
  tryCatch(
    t.test(
      protein_matrix[index, hcc_index],
      protein_matrix[index, adjacent_index],
      paired = TRUE
    )$p.value,
    error = function(error) NA_real_
  )
}, numeric(1))
protein_fdr <- p.adjust(protein_p, method = "BH")

negative_regulators <- c(
  "TSC1", "TSC2", "DEPTOR", "DEPDC5", "NPRL2", "NPRL3", "DDIT4",
  "PTEN", "PRKAG1", "PRKAG2", "PRKAA1", "PRKAB1", "PRKAB2",
  "STK11", "STRADA"
)
delta_minimum <- 0.75
top_n <- 10L
rows <- list()
row_metadata <- data.table()
for (index in seq_len(nrow(protein_matrix))) {
  pathway <- protein_assignment[[index]]
  if (is.na(pathway) || is.na(protein_fdr[[index]])) {
    next
  }
  if (
    protein_fdr[[index]] >= 0.05 ||
      protein_delta[[index]] < delta_minimum
  ) {
    next
  }
  if (
    protein_upper[[index]] %in% negative_regulators &&
      identical(pathway, "mTOR (PI3K-AKT)")
  ) {
    next
  }
  values <- protein_matrix[index, , drop = FALSE]
  rownames(values) <- protein_genes[[index]]
  rows[[length(rows) + 1L]] <- values
  row_metadata <- rbind(
    row_metadata,
    data.table(
      feature = protein_genes[[index]],
      layer = "Protein",
      module = pathway,
      delta = protein_delta[[index]],
      p_adj = protein_fdr[[index]]
    )
  )
}

internal_standards <- c(
  "20-HETE-d6", "Cholesterol(d7)", "7-dehydrocholesterol-d7",
  "Latanoprost (free acid)-d4", "EPA (d5)",
  "8-iso Prostaglandin E2-d4", "Arachidonoyl-EA(d8)",
  "tetranor-PGDM-d6", "Docosahexaenoic Acid-d5",
  "Stearic acid(d3)", "Pentadecylic acid(d3)", "Lignoceric acid(d3)"
)
metabolite_modules <- list(
  "Glycolytic shift" = list(up = c("Lactic Acid")),
  "PLA2 -> LysoPC/LysoPE" = list(up = c(
    "LysoPC(22:1/0:0)", "LysoPC(15:0/0:0)", "LysoPC(17:1/0:0)",
    "LysoPC(20:3/0:0)", "LysoPC(20:2/0:0)", "LysoPE(20:1/0:0)",
    "LysoPE(16:1/0:0)", "LysoPE(0:0/15:0)"
  )),
  "Kynurenine pathway" = list(up = c("Quinolinic acid"))
)

metabolite_statistic <- function(feature) {
  index <- which(metabolite_names == feature)
  if (length(index) == 0L) {
    return(NULL)
  }
  if (length(index) > 1L) {
    stop("Metabolite name is not unique: ", feature)
  }
  values <- metabolite_matrix[index[[1L]], , drop = FALSE]
  list(
    values = values,
    delta = mean(values[hcc_index]) - mean(values[adjacent_index]),
    p_value = tryCatch(
      t.test(
        as.numeric(values)[hcc_index],
        as.numeric(values)[adjacent_index],
        paired = TRUE
      )$p.value,
      error = function(error) NA_real_
    ),
    platform = metabolite_platform[index[[1L]]]
  )
}

standard_qc <- rbindlist(lapply(internal_standards, function(feature) {
  statistic <- metabolite_statistic(feature)
  if (is.null(statistic)) {
    return(NULL)
  }
  data.table(
    standard = feature,
    platform = statistic$platform,
    delta = round(statistic$delta, 3),
    p_value = signif(statistic$p_value, 3),
    flat = abs(statistic$delta) < 0.5
  )
}))
fwrite(
  standard_qc,
  result_path("fig4c", "QC_internal_standard_bias.csv")
)

candidates <- rbindlist(lapply(names(metabolite_modules), function(module) {
  definitions <- metabolite_modules[[module]]
  rbindlist(lapply(
    c(definitions$down, definitions$up),
    function(feature) {
      if (feature %in% internal_standards) {
        stop("Internal standard used as analyte: ", feature)
      }
      statistic <- metabolite_statistic(feature)
      if (is.null(statistic)) {
        return(NULL)
      }
      data.table(
        module = module,
        feature = feature,
        platform = statistic$platform,
        expected = if (feature %in% definitions$up) "up" else "down",
        delta = statistic$delta,
        p_value = statistic$p_value
      )
    }
  ))
}))
candidates[, p_adj := p.adjust(p_value, method = "BH")]
candidates[, retained := (
  p_adj < 0.05 &
    (
      (delta > 0 & expected == "up") |
        (delta < 0 & expected == "down")
    )
)]
fwrite(
  candidates,
  result_path("fig4c", "metabolite_module_screen.csv")
)

for (index in which(candidates$retained)) {
  statistic <- metabolite_statistic(candidates$feature[[index]])
  values <- statistic$values
  rownames(values) <- substr(candidates$feature[[index]], 1, 40)
  rows[[length(rows) + 1L]] <- values
  row_metadata <- rbind(
    row_metadata,
    data.table(
      feature = substr(candidates$feature[[index]], 1, 40),
      layer = "Metab",
      module = candidates$module[[index]],
      delta = candidates$delta[[index]],
      p_adj = candidates$p_adj[[index]]
    )
  )
}

row_metadata[, row_index := seq_len(.N)]
protein_top <- row_metadata[
  layer == "Protein",
  .SD[order(-delta)][seq_len(min(top_n, .N))],
  by = module
]$row_index
keep_rows <- sort(unique(c(
  protein_top,
  row_metadata[layer == "Metab", row_index]
)))
rows <- rows[keep_rows]
row_metadata <- row_metadata[keep_rows]
row_metadata[, row_index := NULL]
row_metadata[, significance := fifelse(
  p_adj < 0.001,
  "***",
  fifelse(p_adj < 0.01, "**", fifelse(p_adj < 0.05, "*", ""))
)]
fwrite(
  row_metadata,
  result_path("fig4c", "selected_features.csv")
)

matrix <- do.call(rbind, rows)
z_matrix <- t(scale(t(matrix)))

saveRDS(
  list(
    values = matrix,
    z_scores = z_matrix,
    samples = common_samples,
    group = group,
    pair = pair,
    features = row_metadata
  ),
  result_path("fig4c", "selected_feature_matrix.rds")
)
fwrite(
  as.data.table(z_matrix, keep.rownames = "feature"),
  result_path("fig4c", "selected_feature_z_scores.csv")
)
