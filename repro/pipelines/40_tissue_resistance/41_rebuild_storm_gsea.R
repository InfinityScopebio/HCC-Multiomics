#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(data.table)
  library(dplyr)
  library(GEOquery)
  library(limma)
  library(msigdbr)
})

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "repro", "config", "project.toml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the project root")
    }
    current <- parent
  }
}

write_tsv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(value, path, sep = "\t", quote = FALSE, na = "NA")
}

root <- find_project_root()
result_dir <- file.path(root, "repo", "results", "panels", "main", "4E")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
audit_dir <- file.path(root, "repo", "audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

series_path <- file.path(
  root,
  "repo", "data", "public", "raw", "GSE109211",
  "GSE109211_series_matrix.txt.gz"
)
platform_path <- file.path(
  root,
  "repo", "data", "public", "raw", "GSE109211",
  "GPL13938_GEO_platform_full.txt"
)
stopifnot(file.exists(series_path), file.exists(platform_path))

storm_seed <- as.integer(Sys.getenv("STORM_GSEA_SEED", unset = "1"))
if (is.na(storm_seed)) {
  stop("STORM_GSEA_SEED must be an integer.")
}
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(storm_seed)
gse <- GEOquery::getGEO(filename = series_path, getGPL = FALSE)
expression <- Biobase::exprs(gse)
platform <- GEOquery::Table(GEOquery::getGEO(filename = platform_path))

probe_index <- match(rownames(expression), platform$ID)
symbol_column <- grep(
  "symbol|gene.?symbol|ILMN_Gene|^Gene$",
  colnames(platform),
  value = TRUE,
  ignore.case = TRUE
)[1]
if (is.na(symbol_column)) {
  stop("No gene-symbol column was found in the local GPL13938 table")
}
probe_qc <- data.frame(
  metric = c(
    "expression_probes",
    "platform_matched",
    "platform_unmatched",
    "nonempty_symbols",
    "selected_symbol_column"
  ),
  value = c(
    nrow(expression),
    sum(!is.na(probe_index)),
    sum(is.na(probe_index)),
    sum(
      !is.na(probe_index) &
        !is.na(platform[[symbol_column]][probe_index]) &
        platform[[symbol_column]][probe_index] != ""
    ),
    symbol_column
  )
)
write_tsv(probe_qc, file.path(result_dir, "storm_probe_mapping_qc.tsv"))

if (anyNA(probe_index)) {
  stop("The local GPL13938 manifest does not cover every expression probe")
}
symbols <- as.character(platform[[symbol_column]][probe_index])
symbols[symbols %in% c("", "---")] <- NA_character_
keep <- !is.na(symbols)
expression <- expression[keep, , drop = FALSE]
symbols <- symbols[keep]

mean_order <- order(-rowMeans(expression, na.rm = TRUE))
expression <- expression[mean_order, , drop = FALSE]
symbols <- symbols[mean_order]
duplicate_symbol <- duplicated(symbols)
expression <- expression[!duplicate_symbol, , drop = FALSE]
rownames(expression) <- symbols[!duplicate_symbol]

if (max(expression, na.rm = TRUE) > 50) {
  expression <- log2(expression - min(expression, na.rm = TRUE) + 1)
}

phenotype <- Biobase::pData(gse)
treatment <- ifelse(
  grepl("Sor", phenotype$characteristics_ch1.1, ignore.case = TRUE),
  "Sor",
  "Plac"
)
outcome <- ifelse(
  grepl("non-resp", phenotype$characteristics_ch1.2, ignore.case = TRUE),
  "NR",
  "R"
)
selected <- which(treatment == "Sor")
expression <- expression[, selected, drop = FALSE]
outcome <- outcome[selected]
group <- factor(outcome, levels = c("R", "NR"))

fit <- limma::lmFit(expression, stats::model.matrix(~group))
fit <- limma::eBayes(fit)
differential <- limma::topTable(
  fit,
  coef = 2,
  number = Inf,
  sort.by = "none"
)
ranked <- sort(
  stats::setNames(differential$t, rownames(differential)),
  decreasing = TRUE
)
write_tsv(
  data.frame(gene_symbol = names(ranked), moderated_t = unname(ranked)),
  file.path(result_dir, "storm_limma_rank.tsv")
)

msigdb <- msigdbr::msigdbr(species = "Homo sapiens")
gene_set <- function(name) {
  unique(msigdb$gene_symbol[msigdb$gs_name == name])
}
sets <- list(
  MAPK_ERK = gene_set("REACTOME_MAPK1_MAPK3_SIGNALING"),
  mTOR = gene_set("HALLMARK_PI3K_AKT_MTOR_SIGNALING"),
  NFKB = gene_set("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  HYPOXIA = gene_set("HALLMARK_HYPOXIA")
)
if (any(lengths(sets) == 0L)) {
  stop("At least one required MSigDB gene set was not found")
}
term_to_gene <- dplyr::bind_rows(lapply(
  names(sets),
  function(name) data.frame(term = name, gene = sets[[name]])
))
write_tsv(
  term_to_gene,
  file.path(result_dir, "storm_msigdb_gene_sets.tsv")
)

BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
set.seed(storm_seed)
gsea <- clusterProfiler::GSEA(
  ranked,
  TERM2GENE = term_to_gene,
  pvalueCutoff = 1,
  seed = TRUE,
  verbose = FALSE,
  minGSSize = 5,
  maxGSSize = 1000
)
gsea_result <- as.data.frame(gsea@result)
write_tsv(
  gsea_result,
  file.path(result_dir, "storm_gsea_results.tsv")
)
saveRDS(gsea, file.path(result_dir, "storm_gsea_result.rds"))


seed_audit <- data.table::data.table(
  pipeline = "tissue_resistance",
  step = c("limma_fit", "clusterProfiler_GSEA"),
  script_path = "repro/pipelines/40_tissue_resistance/41_rebuild_storm_gsea.R",
  language = "R",
  algorithm = c("limma empirical Bayes", "clusterProfiler::GSEA"),
  seed = c(NA_integer_, storm_seed),
  seed_source = c("not_applicable", "historical"),
  rng_kind = c(
    "deterministic",
    "Mersenne-Twister;Inversion;Rejection"
  ),
  status = c("deterministic", "fixed"),
  notes = c(
    "No stochastic operation is used by the fitted linear model.",
    paste(
      "The historical seed is set immediately before GSEA;",
      "clusterProfiler seed=TRUE with BiocParallel SerialParam."
    )
  )
)
write_tsv(
  seed_audit,
  file.path(audit_dir, "Fig4E_STORM_seed_audit.tsv")
)

cat("STORM samples:", ncol(expression), "\n")
cat("Responder / non-responder:", paste(names(table(group)), table(group)), "\n")
print(gsea_result[, c("ID", "setSize", "NES", "pvalue", "p.adjust")])
