#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(data.table)
  library(GSVA)
  library(MOFA2)
  library(msigdbr)
  library(survival)
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

input_path <- cfg$raw_input_path
expression_path <- file.path(cfg$tcga_raw_dir, "LIHC_HiSeqV2.gz")
clinical_path <- file.path(cfg$tcga_raw_dir, "LIHC_clinicalMatrix.tsv")
compound_pathway_path <- cfg$local_kegg_compound_pathway
require_analysis_mofa_boundary(cfg)
require_local_files(
  c(
    input_path,
    cfg$analysis_model_path,
    expression_path,
    clinical_path,
    compound_pathway_path,
    file.path(cfg$m0_dir, "serum_lcms.rds"),
    file.path(cfg$m0_dir, "serum_gcms.rds")
  )
)

clean_feature_names <- function(names) {
  names <- iconv(names, "UTF-8", "ASCII//TRANSLIT", sub = "")
  names <- gsub("[^A-Za-z0-9_().:/+-]", "_", names)
  make.unique(gsub("_+", "_", names))
}

empty_enrichment <- data.table(ID = character(), padj = numeric())

input <- readRDS(input_path)
group <- input$grp
stage <- input$stage
group_levels <- c("Healthy", "HepFibrosis", "Cirrhosis", "HCC")
model <- load_model(cfg$analysis_model_path)
scores <- get_factors(model)[[1]]
weights <- get_weights(model)
variance <- get_variance_explained(model)$r2_per_factor[[1]]
if (is.null(rownames(scores)) || !setequal(rownames(scores), input$common)) {
  stop("MOFA score sample IDs do not match the rebuilt serum input.")
}
scores <- scores[input$common, , drop = FALSE]
if (!identical(rownames(scores), input$common)) {
  stop("Failed to reorder MOFA scores to the rebuilt serum sample order.")
}
invisible(lapply(weights, function(view_weights) {
  missing_factors <- setdiff(colnames(scores), colnames(view_weights))
  if (length(missing_factors) > 0L) {
    stop(
      "MOFA weight view lacks factors: ",
      paste(missing_factors, collapse = ", ")
    )
  }
  TRUE
}))

serum_lcms <- readRDS(file.path(cfg$m0_dir, "serum_lcms.rds"))
serum_gcms <- readRDS(file.path(cfg$m0_dir, "serum_gcms.rds"))
lcms_meta <- as.data.table(serum_lcms$meta)
lcms_matrix <- log2(serum_lcms$mat + 1)
lcms_names <- lcms_meta$Metabolites
lcms_keep <- !is.na(lcms_names) & lcms_names != "" &
  rowMeans(!is.na(lcms_matrix)) >= 0.5
lcms_kegg <- lcms_meta$kegg[lcms_keep]
names(lcms_kegg) <- clean_feature_names(
  make.unique(as.character(lcms_names[lcms_keep]))
)

gcms_names <- serum_gcms$meta$Metabolites
gcms_matrix <- log2(serum_gcms$mat + 1)
gcms_keep <- !is.na(gcms_names) & gcms_names != ""
gcms_matrix <- gcms_matrix[gcms_keep, , drop = FALSE]
gcms_names <- gcms_names[gcms_keep]
gcms_kegg <- serum_gcms$meta$kegg[gcms_keep]
gcms_detected <- rowMeans(!is.na(gcms_matrix)) >= 0.5
gcms_names <- gcms_names[gcms_detected]
gcms_kegg <- gcms_kegg[gcms_detected]
names(gcms_kegg) <- clean_feature_names(
  make.unique(as.character(gcms_names))
)

metabolite_kegg <- c(lcms_kegg, gcms_kegg)
metabolite_kegg <- metabolite_kegg[
  !is.na(metabolite_kegg) & metabolite_kegg != ""
]
metabolite_kegg <- vapply(
  strsplit(metabolite_kegg, "[,; ]+"),
  `[`,
  character(1),
  1
)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
kegg <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:KEGG_LEGACY"
)
reactome <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
)
gene_set_table <- rbind(
  hallmark[, c("gs_name", "gene_symbol")],
  kegg[, c("gs_name", "gene_symbol")],
  reactome[, c("gs_name", "gene_symbol")]
)
term_to_gene <- as.data.frame(gene_set_table)
gene_sets <- split(gene_set_table$gene_symbol, gene_set_table$gs_name)

kegg_source <- unique(as.data.frame(kegg[, c("gs_name", "gs_exact_source")]))
source_to_set <- setNames(kegg_source$gs_name, kegg_source$gs_exact_source)
compound_pathway <- fread(
  compound_pathway_path,
  header = FALSE,
  col.names = c("compound", "pathway")
)
compound_pathway[, compound := sub("cpd:", "", compound)]
compound_pathway[, hsa := sub("path:map", "hsa", pathway)]
compound_pathway <- compound_pathway[hsa %in% names(source_to_set)]
compound_pathway[, gene_set := source_to_set[hsa]]
compound_to_pathway <- as.data.frame(
  compound_pathway[, .(gene_set, compound)]
)

enrich_genes <- function(genes) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  if (length(genes) < 3L) {
    return(empty_enrichment)
  }
  result <- tryCatch(
    enricher(
      genes,
      TERM2GENE = term_to_gene,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      minGSSize = 5,
      maxGSSize = 500
    ),
    error = function(error) NULL
  )
  if (is.null(result)) {
    return(empty_enrichment)
  }
  data <- as.data.frame(result)
  data <- data[data$Count >= 3, , drop = FALSE]
  if (nrow(data) == 0L) {
    return(empty_enrichment)
  }
  data.table(ID = data$ID, padj = data$p.adjust)
}

enrich_metabolites <- function(compounds) {
  compounds <- unique(compounds[!is.na(compounds) & compounds != ""])
  if (length(compounds) < 3L) {
    return(empty_enrichment)
  }
  result <- tryCatch(
    enricher(
      compounds,
      TERM2GENE = compound_to_pathway,
      universe = unique(metabolite_kegg),
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      minGSSize = 3,
      maxGSSize = 500
    ),
    error = function(error) NULL
  )
  if (is.null(result)) {
    return(empty_enrichment)
  }
  data <- as.data.frame(result)
  data <- data[data$Count >= 2, , drop = FALSE]
  if (nrow(data) == 0L) {
    return(empty_enrichment)
  }
  data.table(ID = data$ID, padj = data$p.adjust)
}

factor_pathways <- function(native_weights, factor_name) {
  protein <- native_weights$Protein[, factor_name]
  names(protein) <- rownames(native_weights$Protein)
  phospho <- native_weights$Phospho[, factor_name]
  names(phospho) <- sub(
    "\\.[0-9]+$",
    "",
    rownames(native_weights$Phospho)
  )
  genes <- c(protein, phospho)
  gene_quantiles <- quantile(genes, c(0.20, 0.80))
  up_genes <- names(genes)[genes > gene_quantiles[2] & genes > 0]
  down_genes <- names(genes)[genes < gene_quantiles[1] & genes < 0]

  lcms <- native_weights$LCMS[, factor_name]
  gcms <- native_weights$GCMS[, factor_name]
  metabolites <- c(lcms, gcms)
  compounds <- metabolite_kegg[names(metabolites)]
  metabolite_quantiles <- quantile(metabolites, c(0.20, 0.80))
  up_compounds <- compounds[
    metabolites > metabolite_quantiles[2] & metabolites > 0
  ]
  down_compounds <- compounds[
    metabolites < metabolite_quantiles[1] & metabolites < 0
  ]

  up_gene_enrichment <- enrich_genes(up_genes)
  down_gene_enrichment <- enrich_genes(down_genes)
  up_metabolite_enrichment <- enrich_metabolites(up_compounds)
  down_metabolite_enrichment <- enrich_metabolites(down_compounds)
  up_pathways <- unique(c(
    up_gene_enrichment$ID,
    up_metabolite_enrichment$ID
  ))
  down_pathways <- unique(c(
    down_gene_enrichment$ID,
    down_metabolite_enrichment$ID
  ))
  overlap <- intersect(up_pathways, down_pathways)
  combined <- rbind(
    up_gene_enrichment,
    down_gene_enrichment,
    up_metabolite_enrichment,
    down_metabolite_enrichment
  )
  p_map <- if (nrow(combined) > 0L) {
    minimum <- combined[, .(padj = min(padj)), by = ID]
    setNames(minimum$padj, minimum$ID)
  } else {
    numeric()
  }
  list(
    up = setdiff(up_pathways, overlap),
    down = setdiff(down_pathways, overlap),
    padj = p_map,
    up_gene = up_gene_enrichment$ID,
    up_metabolite = up_metabolite_enrichment$ID,
    down_gene = down_gene_enrichment$ID,
    down_metabolite = down_metabolite_enrichment$ID
  )
}

expression <- fread(
  cmd = paste("gunzip -c", shQuote(expression_path)),
  header = TRUE
)
genes <- expression[[1]]
expression_matrix <- as.matrix(expression[, -1])
rownames(expression_matrix) <- genes
clinical <- fread(clinical_path)
clinical[, OS_time := ifelse(
  vital_status == "DECEASED",
  as.numeric(days_to_death),
  as.numeric(days_to_last_followup)
)]
clinical[, OS := as.integer(vital_status == "DECEASED")]
clinical[, stage_num := as.numeric(factor(
  pathologic_stage,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage IIIA",
    "Stage IIIB",
    "Stage IIIC",
    "Stage IV",
    "Stage IVA",
    "Stage IVB"
  )
))]
clinical[, age_num := as.numeric(age_at_initial_pathologic_diagnosis)]
clinical[, sex := factor(gender)]
clinical[, grade_num := as.numeric(factor(
  neoplasm_histologic_grade,
  levels = c("G1", "G2", "G3", "G4")
))]
clinical[, stage_rom := gsub(
  "^Stage ([IVX]+).*",
  "\\1",
  pathologic_stage
)]
clinical[!stage_rom %in% c("I", "II", "III", "IV"), stage_rom := NA]
clinical[stage_rom %in% c("III", "IV"), stage_rom := "III-IV"]
survival <- clinical[
  sample_type == "Primary Tumor" &
    !is.na(OS_time) &
    OS_time > 0
]

project_directional_score <- function(up_pathways, down_pathways) {
  pathways <- intersect(
    unique(c(up_pathways, down_pathways)),
    names(gene_sets)
  )
  if (length(pathways) < 2L) {
    return(NULL)
  }
  parameters <- ssgseaParam(
    expression_matrix,
    gene_sets[pathways],
    normalize = TRUE,
    minSize = 5
  )
  activity <- suppressWarnings(suppressMessages(
    gsva(parameters, verbose = FALSE)
  ))
  up_present <- intersect(up_pathways, rownames(activity))
  down_present <- intersect(down_pathways, rownames(activity))
  up_score <- if (length(up_present) > 0L) {
    colMeans(t(scale(t(activity[up_present, , drop = FALSE]))))
  } else {
    0
  }
  down_score <- if (length(down_present) > 0L) {
    colMeans(t(scale(t(activity[down_present, , drop = FALSE]))))
  } else {
    0
  }
  list(score = up_score - down_score, up = up_score, down = down_score)
}

axis_store <- list()
pathway_rows <- list()
factor_rows <- list()
native_factor_rows <- list()
for (factor_name in colnames(scores)) {
  native_factor_score <- scores[, factor_name]
  native_group_means <- tapply(
    native_factor_score,
    group,
    mean
  )[group_levels]
  native_rho <- cor(native_factor_score, stage, method = "spearman")
  factor_score <- native_factor_score
  group_means <- native_group_means
  pathways <- factor_pathways(weights, factor_name)
  projection <- project_directional_score(pathways$up, pathways$down)
  if (is.null(projection)) {
    next
  }
  projected <- merge(
    survival,
    data.table(
      sampleID = names(projection$score),
      score = as.numeric(scale(projection$score)),
      up = as.numeric(scale(projection$up)),
      dn = as.numeric(scale(projection$down))
    ),
    by = "sampleID"
  )
  univariable <- summary(coxph(
    Surv(OS_time, OS) ~ score,
    data = projected
  ))
  factor_rows[[factor_name]] <- data.table(
    factor = factor_name,
    VE_total = round(sum(variance[factor_name, ]), 2),
    H = group_means[1],
    HF = group_means[2],
    Cir = group_means[3],
    HCC = group_means[4],
    strict = all(diff(group_means) > 0),
    HR_uni = univariable$conf.int["score", "exp(coef)"],
    HR_lo = univariable$conf.int["score", "lower .95"],
    HR_hi = univariable$conf.int["score", "upper .95"],
    p_uni = univariable$coefficients["score", "Pr(>|z|)"]
  )
  native_factor_rows[[factor_name]] <- data.table(
    factor = factor_name,
    H = native_group_means[1],
    HF = native_group_means[2],
    Cir = native_group_means[3],
    HCC = native_group_means[4],
    rho_stage = native_rho
  )
  pathway_rows[[factor_name]] <- rbind(
    if (length(pathways$up) > 0L) {
      data.table(
        factor = factor_name,
        direction = "UP",
        pathway = pathways$up,
        source = fifelse(
          pathways$up %in% pathways$up_gene &
            pathways$up %in% pathways$up_metabolite,
          "gene+metab",
          fifelse(
            pathways$up %in% pathways$up_metabolite,
            "metab",
            "gene"
          )
        ),
        padj = as.numeric(pathways$padj[pathways$up])
      )
    },
    if (length(pathways$down) > 0L) {
      data.table(
        factor = factor_name,
        direction = "DOWN",
        pathway = pathways$down,
        source = fifelse(
          pathways$down %in% pathways$down_gene &
            pathways$down %in% pathways$down_metabolite,
          "gene+metab",
          fifelse(
            pathways$down %in% pathways$down_metabolite,
            "metab",
            "gene"
          )
        ),
        padj = as.numeric(pathways$padj[pathways$down])
      )
    }
  )
  axis_store[[factor_name]] <- list(
    native_factor_score = native_factor_score,
    factor_score = factor_score,
    score_orientation = "native",
    stage_rho = native_rho,
    W = weights,
    up_pathways = pathways$up,
    down_pathways = pathways$down,
    tcga_df = projected,
    stats = as.list(factor_rows[[factor_name]])
  )
  cat(
    factor_name,
    "pathways",
    length(pathways$up),
    "UP",
    length(pathways$down),
    "DOWN\n"
  )
}

factor_order <- paste0("Factor", seq_len(15))
axis_store <- axis_store[factor_order]
pathway_table <- rbindlist(pathway_rows[factor_order])
pathway_table[, factor_index := match(factor, factor_order)]
setorder(pathway_table, factor_index, direction, pathway)
pathway_table[, factor_index := NULL]
factor_table <- rbindlist(factor_rows[factor_order])
native_factor_table <- rbindlist(native_factor_rows[factor_order])
saveRDS(
  axis_store,
  file.path(cfg$result_dir, "incdec_K15_all_axes_rebuilt.rds")
)
fwrite(
  pathway_table,
  file.path(cfg$result_dir, "incdec_K15_allfactors_pathways_rebuilt.tsv"),
  sep = "\t"
)
fwrite(
  factor_table,
  file.path(cfg$result_dir, "incdec_allfactors_TCGA_rebuilt.tsv"),
  sep = "\t"
)
fwrite(
  native_factor_table,
  file.path(cfg$result_dir, "incdec_allfactors_serum_native_rebuilt.tsv"),
  sep = "\t"
)

# Validation-only targets are read after all production outputs are materialized.
# They never supply weights, loadings, pathways or scores to the computation above.
require_local_files(c(cfg$validation$axis_store, cfg$validation$pathways))
historical_axes <- readRDS(cfg$validation$axis_store)
score_magnitude_differences <- vapply(factor_order, function(factor_name) {
  expected <- as.data.table(
    historical_axes[[factor_name]]$tcga_df
  )[, .(sampleID, score)]
  observed <- as.data.table(
    axis_store[[factor_name]]$tcga_df
  )[, .(sampleID, score)]
  setorder(expected, sampleID)
  setorder(observed, sampleID)
  if (!identical(expected$sampleID, observed$sampleID)) {
    return(Inf)
  }
  max(abs(abs(expected$score) - abs(observed$score)))
}, numeric(1))

historical_pathways <- fread(cfg$validation$pathways)
historical_pathways[, factor_index := match(factor, factor_order)]
setorder(historical_pathways, factor_index, pathway, source)
historical_pathways[, factor_index := NULL]
native_pathways_for_validation <- copy(pathway_table)
native_pathways_for_validation[, factor_index := match(factor, factor_order)]
setorder(native_pathways_for_validation, factor_index, pathway, source)
native_pathways_for_validation[, factor_index := NULL]
pathway_membership_exact <- identical(
  historical_pathways[, .(factor, pathway, source)],
  native_pathways_for_validation[, .(factor, pathway, source)]
)
pathway_padj_difference <- if (pathway_membership_exact) {
  max(abs(
    historical_pathways$padj - native_pathways_for_validation$padj
  ))
} else {
  Inf
}

audit <- data.table(
  check_id = c(
    "tcga_native_score_magnitudes",
    "native_pathway_membership_undirected",
    "native_pathway_padj",
    "local_only_execution"
  ),
  status = c(
    ifelse(
      max(score_magnitude_differences) <= 1e-12,
      "numeric_exact",
      "numeric_mismatch"
    ),
    ifelse(pathway_membership_exact, "exact", "mismatch"),
    ifelse(pathway_padj_difference <= 1e-12, "numeric_exact", "numeric_mismatch"),
    "exact"
  ),
  observed = c(
    sprintf(
      "max_abs_magnitude_diff=%.3g",
      max(score_magnitude_differences)
    ),
    sprintf(
      "historical_rows=%d; rebuilt_rows=%d",
      nrow(historical_pathways),
      nrow(pathway_table)
    ),
    sprintf("max_abs_diff=%.3g", pathway_padj_difference),
    paste(
      "TCGA expression, clinical, serum M0, K15 HDF5, kegg cache and",
      "installed msigdbr data only"
    )
  ),
  note = c(
    paste(
      "Production scores use native MOFA signs without multiplying a stage or",
      "progression direction. The legacy aligned target is checked by score",
      "magnitude only."
    ),
    paste(
      "Native UP/DOWN 20% loading tails and overlap removal follow step49g;",
      "legacy membership is checked without importing its aligned direction."
    ),
    paste(
      "Direction-invariant pathway membership is identical before comparing",
      "clusterProfiler BH values at full precision."
    ),
    paste0(
      "No downloader or network API is called; production model SHA256=",
      sha256_file(cfg$analysis_model_path)
    )
  )
)
fwrite(
  audit,
  file.path(cfg$audit_dir, "mofa_tcga_projection_audit.tsv"),
  sep = "\t"
)

cat(
  sprintf(
    paste0(
      "TCGA projection completed from the selected K15 analysis model; ",
      "validation ",
      "native-score magnitude max_abs_diff=%.12g, ",
      "native pathway_membership_exact=%s.\n"
    ),
    max(score_magnitude_differences),
    pathway_membership_exact
  )
)
