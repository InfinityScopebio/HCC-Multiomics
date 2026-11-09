#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(matrixStats)
  library(proDA)
})

script_path <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

phospho <- readRDS(data_path("tissue_phospho.rds"))
metadata <- paired_tissue_metadata()
patients <- metadata$patients
hcc_columns <- paste0("Can", patients)
adjacent_columns <- paste0("Para", patients)

matrix <- phospho$mat[, c(hcc_columns, adjacent_columns), drop = FALSE]
if (any(duplicated(rownames(matrix)))) {
  rownames(matrix) <- make.unique(rownames(matrix), sep = "__")
}
detected_hcc <- rowMeans(!is.na(matrix[, hcc_columns, drop = FALSE]))
detected_adjacent <- rowMeans(
  !is.na(matrix[, adjacent_columns, drop = FALSE])
)
keep <- detected_hcc >= 0.30 | detected_adjacent >= 0.30
filtered <- matrix[keep, , drop = FALSE]
if (nrow(filtered) > 2000L) {
  variance <- rowVars(log2(filtered + 1), na.rm = TRUE)
  filtered <- filtered[order(-variance)[seq_len(2000L)], , drop = FALSE]
}
filtered <- log2(filtered + 1)

group <- factor(
  ifelse(colnames(filtered) %in% hcc_columns, "HCC", "Adj"),
  levels = c("Adj", "HCC")
)
patient <- factor(sub("^(Can|Para)", "", colnames(filtered)))
design <- model.matrix(~ patient + group)
fit <- proDA(filtered, design = design)
proda <- as.data.table(test_diff(fit, "groupHCC"))
proda[, site := name]
fwrite(
  proda,
  result_path("ksea", "site_paired_DE_proDA.csv")
)

phospho_meta <- as.data.table(phospho$meta)
phospho_meta[, site_id_full := make.unique(site_id, sep = "__")]
proda <- merge(
  proda,
  phospho_meta[, .(
    site_id_full,
    SYMBOL,
    SiteAA,
    SiteLocation,
    FlankingRegion,
    Gene
  )],
  by.x = "site",
  by.y = "site_id_full",
  all.x = TRUE
)
proda[, sub_site := paste0(SYMBOL, "_", SiteAA, SiteLocation)]
proda <- proda[!is.na(SYMBOL) & SYMBOL != "" & !is.na(diff)]

prior <- as.data.table(readRDS(reference_path(
  "OmniPath", "omnipath_kinase_substrate.rds"
)))
if ("modification" %in% colnames(prior)) {
  prior <- prior[tolower(modification) == "phosphorylation"]
}
prior[, sub_site := paste0(
  substrate_genesymbol, "_", residue_type, residue_offset
)]
prior[, kinase := enzyme_genesymbol]
prior <- unique(prior[
  !is.na(kinase) & kinase != "" & !is.na(sub_site),
  .(kinase, sub_site)
])

global_logfc <- proda$diff
global_mean <- mean(global_logfc)
global_sd <- sd(global_logfc)
lookup <- setNames(proda$diff, proda$sub_site)
ksea <- prior[, {
  values <- lookup[intersect(sub_site, names(lookup))]
  values <- values[!is.na(values)]
  if (length(values) >= 3L) {
    z_score <- (
      (mean(values) - global_mean) * sqrt(length(values)) / global_sd
    )
    list(
      n_sub_data = length(values),
      mean_logFC = mean(values),
      ksea_z = z_score,
      ksea_p = 2 * pnorm(-abs(z_score))
    )
  } else {
    list(
      n_sub_data = length(values),
      mean_logFC = NA_real_,
      ksea_z = NA_real_,
      ksea_p = NA_real_
    )
  }
}, by = kinase]
ksea[, ksea_fdr := p.adjust(ksea_p, method = "BH")]
ksea[, abs_z := abs(ksea_z)]
setorder(ksea, -abs_z, na.last = TRUE)
fwrite(
  ksea,
  result_path("ksea", "kinase_activity_ksea.csv")
)

differential <- fread(result_path(
  "cross_layer", "paired_DE_per_layer.csv"
))
rna_de <- differential[layer == "rna"]
protein_de <- differential[layer == "protein"]
rna_de[, gene := SYMBOL]
protein_de[, gene := SYMBOL]

three_layer <- data.table(kinase = ksea[!is.na(ksea_z), kinase])
three_layer <- merge(
  three_layer,
  ksea[, .(
    kinase,
    n_sub_data,
    mean_logFC_phospho = mean_logFC,
    ksea_z,
    ksea_fdr
  )],
  by = "kinase"
)
three_layer <- merge(
  three_layer,
  rna_de[, .(
    kinase = gene,
    mRNA_logFC = logFC,
    mRNA_p = adj.P.Val
  )],
  by = "kinase",
  all.x = TRUE
)
three_layer <- merge(
  three_layer,
  protein_de[, .(
    kinase = gene,
    prot_logFC = logFC,
    prot_p = adj.P.Val
  )],
  by = "kinase",
  all.x = TRUE
)
three_layer[, sign_phospho := sign(ksea_z)]
three_layer[, sign_mRNA := sign(mRNA_logFC)]
three_layer[, sign_prot := sign(prot_logFC)]
three_layer[, three_layer_consistent := (
  !is.na(sign_mRNA) &
    !is.na(sign_prot) &
    sign_phospho == sign_mRNA &
    sign_phospho == sign_prot
)]
fwrite(
  three_layer,
  result_path("ksea", "kinase_three_layer.csv")
)

message(
  "KSEA complete: ",
  sum(!is.na(ksea$ksea_z)),
  " tested kinases; ",
  sum(ksea$ksea_fdr < 0.05, na.rm = TRUE),
  " at FDR < 0.05."
)
