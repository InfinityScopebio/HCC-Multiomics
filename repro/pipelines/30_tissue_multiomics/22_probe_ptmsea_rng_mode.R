#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BiocParallel)
  library(clusterProfiler)
  library(data.table)
  library(org.Hs.eg.db)
})

script_path <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
project_root <- normalizePath(
  file.path(dirname(normalizePath(script_path)), "..", "..", ".."),
  mustWork = TRUE
)

parse_argument <- function(name) {
  prefix <- paste0("--", name, "=")
  match <- grep(
    paste0("^", prefix),
    commandArgs(trailingOnly = TRUE),
    value = TRUE
  )
  if (length(match) != 1L) {
    stop("Expected exactly one ", prefix, " argument.")
  }
  sub(prefix, "", match, fixed = TRUE)
}

mode <- parse_argument("mode")
replicate_id <- suppressWarnings(as.integer(parse_argument("replicate")))
valid_modes <- c(
  "historical_default",
  "historical_serial",
  "explicit_seed_default",
  "explicit_seed_serial"
)
if (!mode %in% valid_modes) {
  stop("Unknown mode: ", mode)
}
if (is.na(replicate_id) || replicate_id < 1L) {
  stop("replicate must be a positive integer.")
}

gmt_path <- file.path(
  project_root,
  "repo",
  "data",
  "reference",
  "PTMsigDB",
  "ptm.sig.db.all.uniprot.human.v2.0.0.gmt"
)
protein_path <- file.path(
  project_root,
  "repo",
  "results",
  "00_preprocess",
  "multiomics",
  "dat",
  "tissue_protein.rds"
)
phospho_path <- file.path(
  project_root,
  "repo",
  "results",
  "00_preprocess",
  "multiomics",
  "dat",
  "tissue_phospho.rds"
)
differential_path <- file.path(
  project_root,
  "repo",
  "results",
  "30_tissue_multiomics",
  "ptmsea",
  "phosphosite_paired_limma.csv"
)
required_paths <- c(
  gmt_path,
  protein_path,
  phospho_path,
  differential_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing PTM-SEA probe input(s): ", paste(missing_paths, collapse = ", "))
}

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
missing_accessions <- setdiff(names(organism_lookup), names(uniprot_to_symbol))
uniprot_to_symbol <- c(
  uniprot_to_symbol,
  organism_lookup[missing_accessions]
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

differential <- fread(differential_path)
required_columns <- c("logFC", "site_key")
if (!all(required_columns %in% colnames(differential))) {
  stop("Paired phosphosite table lacks logFC/site_key.")
}
ranked <- differential$logFC
names(ranked) <- differential$site_key
ranked <- ranked[
  order(-ranked, seq_along(ranked), method = "radix")
]

use_serial <- endsWith(mode, "_serial")
use_internal_seed <- startsWith(mode, "explicit_seed_")
if (use_serial) {
  register(SerialParam(), default = TRUE)
}
backend <- class(registered()[[1L]])[[1L]]

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42)
gsea <- clusterProfiler::GSEA(
  ranked,
  TERM2GENE = term_to_gene,
  pvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 2000,
  verbose = FALSE,
  seed = use_internal_seed,
  eps = 0
)
results <- as.data.table(gsea@result)
results[, category := sub("-.*", "", ID)]
setorder(results, ID)

sha256_file <- function(path) {
  unname(tools::sha256sum(path))
}

output_root <- file.path(
  project_root,
  "repo",
  "audit",
  "ptmsea_rng_modes"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
stem <- sprintf("%s_repeat%d", mode, replicate_id)
result_path <- file.path(output_root, paste0(stem, ".csv"))
worker_path <- normalizePath(script_path, mustWork = TRUE)
production_script_path <- file.path(
  project_root,
  "repro",
  "pipelines",
  "30_tissue_multiomics",
  "20_ptmsea.R"
)
orchestrator_script_path <- file.path(
  project_root,
  "repro",
  "pipelines",
  "30_tissue_multiomics",
  "22_run_ptmsea_rng_probe.R"
)
summarizer_script_path <- file.path(
  project_root,
  "repro",
  "pipelines",
  "30_tissue_multiomics",
  "23_summarize_ptmsea_rng_modes.R"
)
run_all_script_path <- file.path(
  project_root,
  "repro",
  "pipelines",
  "30_tissue_multiomics",
  "run_all.sh"
)
project_config_path <- file.path(
  project_root,
  "repro",
  "config",
  "project.toml"
)
input_manifest_path <- file.path(
  project_root,
  "repro",
  "manifests",
  "input_sources.tsv"
)
fwrite(results, result_path)
fwrite(
  data.table(
    mode = mode,
    replicate = replicate_id,
    external_seed = 42L,
    GSEA_seed_argument = use_internal_seed,
    registered_backend = backend,
    signatures_tested = nrow(results),
    signatures_fdr_0_05 = sum(results$p.adjust < 0.05, na.rm = TRUE),
    gmt_sha256 = sha256_file(gmt_path),
    protein_sha256 = sha256_file(protein_path),
    phospho_sha256 = sha256_file(phospho_path),
    differential_sha256 = sha256_file(differential_path),
    worker_script_sha256 = sha256_file(worker_path),
    production_script_sha256 = sha256_file(production_script_path),
    orchestrator_script_sha256 = sha256_file(orchestrator_script_path),
    summarizer_script_sha256 = sha256_file(summarizer_script_path),
    run_all_script_sha256 = sha256_file(run_all_script_path),
    project_config_sha256 = sha256_file(project_config_path),
    input_manifest_sha256 = sha256_file(input_manifest_path),
    result_sha256 = sha256_file(result_path),
    clusterProfiler_version = as.character(packageVersion("clusterProfiler")),
    fgsea_version = as.character(packageVersion("fgsea")),
    BiocParallel_version = as.character(packageVersion("BiocParallel")),
    R_version = R.version.string
  ),
  file.path(output_root, paste0(stem, ".metadata.tsv")),
  sep = "\t"
)

message(
  mode,
  " repeat ",
  replicate_id,
  ": ",
  nrow(results),
  " signatures, ",
  sum(results$p.adjust < 0.05, na.rm = TRUE),
  " FDR<0.05, backend=",
  backend
)
