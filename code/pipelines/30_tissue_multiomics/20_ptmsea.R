#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(limma)
  library(org.Hs.eg.db)
})

script_path <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

gmt_path <- reference_path(
  "PTMsigDB", "ptm.sig.db.all.uniprot.human.v2.0.0.gmt"
)
protein_path <- data_path("tissue_protein.rds")
phospho_path <- data_path("tissue_phospho.rds")
assert_files(c(gmt_path, protein_path, phospho_path))

lines <- readLines(gmt_path)
signatures <- list()
for (line in lines) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
  sites <- parts[-c(1L, 2L)]
  signatures[[parts[[1L]]]] <- sites[nzchar(sites)]
}

protein <- readRDS(protein_path)
phospho <- readRDS(phospho_path)
uniprot_to_symbol <- c(
  setNames(protein$meta$Gene, protein$meta$Accession),
  setNames(phospho$meta$SYMBOL, phospho$meta$ProteinAccession)
)
uniprot_to_symbol <- uniprot_to_symbol[
  !is.na(uniprot_to_symbol) & uniprot_to_symbol != ""
]
uniprot_to_symbol <- uniprot_to_symbol[
  !duplicated(names(uniprot_to_symbol))
]

uniprot_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "UNIPROT")
organism_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = uniprot_keys,
  keytype = "UNIPROT",
  columns = c("UNIPROT", "SYMBOL")
))
organism_map <- organism_map[!is.na(organism_map$SYMBOL), ]
organism_lookup <- setNames(organism_map$SYMBOL, organism_map$UNIPROT)
missing_accessions <- setdiff(
  names(organism_lookup), names(uniprot_to_symbol)
)
uniprot_to_symbol <- c(
  uniprot_to_symbol, organism_lookup[missing_accessions]
)

parse_site <- function(site) {
  parts <- strsplit(site, ";", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) {
    return(NA_character_)
  }
  accession <- parts[[1L]]
  position <- sub("^[A-Z]([0-9]+).*", "\\1", parts[[2L]])
  symbol <- unname(uniprot_to_symbol[accession])
  if (length(symbol) == 0L || is.na(symbol) || symbol == "") {
    return(NA_character_)
  }
  paste0(symbol, "_", position)
}

term_to_gene <- rbindlist(lapply(names(signatures), function(name) {
  keys <- vapply(signatures[[name]], parse_site, character(1))
  keys <- unique(keys[!is.na(keys)])
  if (length(keys) == 0L) {
    return(NULL)
  }
  data.table(gs_name = name, gene_symbol = keys)
}))

sample_metadata <- paired_tissue_metadata()
tissue <- sample_metadata$tissue
sample_ids <- colnames(phospho$mat)
groups <- tissue$group_4[match(sample_ids, tissue$sample_id)]
pairs <- tissue$paired_id[match(sample_ids, tissue$sample_id)]
keep_samples <- (
  !is.na(groups) & !is.na(pairs) & groups %in% c("HCC", "Adjacent")
)
matrix <- phospho$mat[, keep_samples, drop = FALSE]
groups <- groups[keep_samples]
pairs <- pairs[keep_samples]
rownames(matrix) <- phospho$meta$site_id

keep_sites <- rowMeans(!is.na(matrix)) >= 0.5
matrix <- matrix[keep_sites, , drop = FALSE]
matrix <- log2(matrix + 1)
sample_medians <- apply(matrix, 2, median, na.rm = TRUE)
matrix <- sweep(
  matrix,
  2,
  sample_medians - median(sample_medians, na.rm = TRUE),
  "-"
)
for (index in seq_len(nrow(matrix))) {
  missing <- is.na(matrix[index, ])
  if (any(missing) && !all(missing)) {
    matrix[index, missing] <- min(matrix[index, ], na.rm = TRUE)
  }
}
matrix[is.na(matrix)] <- 0

group_factor <- factor(groups, levels = c("Adjacent", "HCC"))
pair_factor <- factor(pairs)
design <- model.matrix(~ pair_factor + group_factor)
fit <- eBayes(lmFit(matrix, design), robust = TRUE)
differential <- as.data.table(topTable(
  fit,
  coef = "group_factorHCC",
  number = Inf,
  sort.by = "none"
))
differential[, SYMBOL := phospho$meta$SYMBOL[keep_sites]]
differential[, site_pos := as.character(
  phospho$meta$SiteLocation[keep_sites]
)]
differential <- differential[
  !is.na(SYMBOL) & SYMBOL != "" & !is.na(site_pos)
]
differential[, site_key := paste0(SYMBOL, "_", site_pos)]
differential[, source_order := .I]
differential <- differential[
  order(-abs(logFC), source_order, method = "radix")
][!duplicated(site_key)]
differential[, source_order := NULL]

ranked <- differential$logFC
names(ranked) <- differential$site_key
ranked <- ranked[
  order(-ranked, seq_along(ranked), method = "radix")
]

BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
run_ptmsea <- function() {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(42)
  clusterProfiler::GSEA(
    ranked,
    TERM2GENE = term_to_gene,
    pvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 2000,
    verbose = FALSE,
    seed = FALSE,
    eps = 0
  )
}

results <- as.data.table(run_ptmsea()@result)
repeat_results <- as.data.table(run_ptmsea()@result)
setorder(results, ID)
setorder(repeat_results, ID)
if (!identical(results, repeat_results)) {
  stop("PTM-SEA fixed-seed serial repeat was not fully identical.")
}
results[, category := sub("-.*", "", ID)]
fwrite(
  results,
  result_path("ptmsea", "PTMSEA_official_v2.0.0_all.csv")
)
fwrite(
  differential,
  result_path("ptmsea", "phosphosite_paired_limma.csv")
)
fwrite(
  data.table(
    signatures_parsed = length(signatures),
    mapped_term_site_edges = nrow(term_to_gene),
    ranked_sites = length(ranked),
    tied_ranked_sites = sum(duplicated(unname(ranked))),
    tie_order = "stable_original_matrix_order",
    overlap_sites = length(intersect(names(ranked), term_to_gene$gene_symbol)),
    signatures_tested = nrow(results),
    signatures_fdr_0_05 = sum(results$p.adjust < 0.05, na.rm = TRUE),
    rng_seed = 42L,
    GSEA_seed_argument = FALSE,
    registered_backend = "SerialParam",
    repeat_exact_full_result = TRUE,
    clusterProfiler_version =
      as.character(packageVersion("clusterProfiler")),
    fgsea_version = as.character(packageVersion("fgsea"))
  ),
  result_path("ptmsea", "PTMSEA_run_summary.tsv"),
  sep = "\t"
)

message(
  "PTM-SEA complete: ", nrow(results), " signatures tested; ",
  sum(results$p.adjust < 0.05, na.rm = TRUE), " at FDR < 0.05."
)
