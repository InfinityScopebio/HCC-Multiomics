#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(matrixStats)
})

script_path <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

metadata <- paired_tissue_metadata()
patients <- metadata$patients

paired_differential <- function(object, log2_transform, layer) {
  matrix <- object$mat
  if (log2_transform) {
    matrix <- log2(matrix + 1)
  }
  paired_patients <- patients[
    paste0("Can", patients) %in% colnames(matrix) &
      paste0("Para", patients) %in% colnames(matrix)
  ]
  hcc <- matrix[, paste0("Can", paired_patients), drop = FALSE]
  adjacent <- matrix[, paste0("Para", paired_patients), drop = FALSE]
  keep <- !(rowAnys(is.na(hcc)) & rowAnys(is.na(adjacent)))
  hcc <- hcc[keep, , drop = FALSE]
  adjacent <- adjacent[keep, , drop = FALSE]
  difference <- hcc - adjacent
  design <- matrix(1, ncol(difference), 1)
  colnames(design) <- "delta"
  fit <- eBayes(lmFit(difference, design))
  results <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
  results$feature <- rownames(difference)
  results$layer <- layer
  if ("SYMBOL" %in% colnames(object$meta)) {
    results$SYMBOL <- object$meta$SYMBOL[
      match(rownames(difference), object$feature_id)
    ]
  }
  if (identical(layer, "rna")) {
    results$SYMBOL <- rownames(difference)
  }
  results$contrast <- "HCC_vs_Adj"
  as.data.table(results)
}

rna <- readRDS(data_path("tissue_rna.rds"))
protein <- readRDS(data_path("tissue_protein.rds"))
metabolite <- readRDS(data_path("tissue_metabolite.rds"))

rna_de <- paired_differential(rna, TRUE, "rna")
protein_de <- paired_differential(protein, TRUE, "protein")
metabolite_de <- paired_differential(metabolite, TRUE, "metabolite")
differential <- rbindlist(
  list(rna_de, protein_de, metabolite_de),
  fill = TRUE
)
fwrite(
  differential,
  result_path("cross_layer", "paired_DE_per_layer.csv")
)

score_genes <- function(results, layer) {
  results <- results[
    !is.na(SYMBOL) & !is.na(logFC) & !is.na(P.Value)
  ]
  results[, signed_score := logFC * (-log10(P.Value))]
  results[, .(layer = layer, gene = SYMBOL, score = signed_score)]
}

rna_scores <- score_genes(rna_de, "rna")
protein_scores <- score_genes(protein_de, "protein")

gene_sets <- as.data.table(readRDS(result_path(
  "references", "msigdb_tissue_sets.rds"
)))
gene_sets <- gene_sets[
  gs_collection == "C2" &
    gs_subcollection %in% c("CP:REACTOME", "CP:KEGG_LEGACY")
]
pathway_genes <- gene_sets[
  ,
  .(path = gs_id, SYMBOL = gene_symbol, name = gs_name)
]
pathway_names <- unique(pathway_genes[, .(path, name)])

calculate_pathway_scores <- function(gene_scores, layer) {
  paths <- unique(pathway_genes$path)
  global_mean <- mean(gene_scores$score, na.rm = TRUE)
  global_sd <- sd(gene_scores$score, na.rm = TRUE)
  output <- rbindlist(lapply(paths, function(path_id) {
    genes <- pathway_genes[path == path_id, SYMBOL]
    values <- gene_scores[gene %in% genes, score]
    values <- values[!is.na(values)]
    if (length(values) < 3L) {
      return(NULL)
    }
    data.table(
      layer = layer,
      path = path_id,
      n_genes = length(values),
      z = (mean(values) - global_mean) / (
        global_sd / sqrt(length(values))
      )
    )
  }))
  merge(output, pathway_names, by = "path", all.x = TRUE)
}

rna_pathways <- calculate_pathway_scores(rna_scores, "rna")
protein_pathways <- calculate_pathway_scores(protein_scores, "protein")
pathway_scores <- rbindlist(list(rna_pathways, protein_pathways))
fwrite(
  pathway_scores,
  result_path("cross_layer", "pathway_zscore_per_layer.csv")
)

shared <- intersect(rna_pathways$path, protein_pathways$path)
x <- rna_pathways[match(shared, path), z]
y <- protein_pathways[match(shared, path), z]
complete <- complete.cases(x, y)
x <- x[complete]
y <- y[complete]
shared <- shared[complete]
rho <- cor(x, y, method = "spearman")
stopifnot(is.finite(rho), sd(x) > 0, sd(y) > 0)

set.seed(42)
bootstrap <- replicate(1000, {
  index <- sample(seq_along(x), replace = TRUE)
  cor(x[index], y[index], method = "spearman")
})
confidence_interval <- quantile(
  bootstrap, c(0.025, 0.975), na.rm = TRUE
)
correlation <- data.table(
  layer_x = "rna",
  layer_y = "protein",
  rho_point = rho,
  ci_low = confidence_interval[[1L]],
  ci_high = confidence_interval[[2L]],
  n_pathways = length(x)
)
fwrite(
  correlation,
  result_path("cross_layer", "cross_layer_rho.csv")
)

pathway_scores[, p := 2 * pnorm(-abs(z))]
pathway_scores[, fdr := p.adjust(p, "BH"), by = layer]
significant_rna <- rna_pathways[
  ,
  .(path, z_rna = z, n_rna = n_genes)
][pathway_scores[layer == "rna" & fdr < 0.05], on = "path", nomatch = 0]
significant_protein <- protein_pathways[
  ,
  .(path, z_prot = z, n_prot = n_genes)
][
  pathway_scores[layer == "protein" & fdr < 0.05],
  on = "path",
  nomatch = 0
]
triangle <- merge(significant_rna, significant_protein, by = "path")
triangle[, same_direction := sign(z_rna) == sign(z_prot)]
triangle <- merge(triangle, pathway_names, by = "path", all.x = TRUE)
fwrite(
  triangle[, .(path, name, z_rna, z_prot, same_direction)],
  result_path("cross_layer", "triangle_concordance.csv")
)

wide <- dcast(pathway_scores, path + name ~ layer, value.var = "z")
wide <- wide[complete.cases(rna, protein)]
triangle_status <- triangle[, .(
  path,
  dual_significant = TRUE,
  same_direction
)]
wide <- merge(wide, triangle_status, by = "path", all.x = TRUE)
wide[is.na(dual_significant), dual_significant := FALSE]
wide[
  dual_significant == TRUE,
  status := ifelse(same_direction, "Dual-significant, concordant", NA_character_)
]
