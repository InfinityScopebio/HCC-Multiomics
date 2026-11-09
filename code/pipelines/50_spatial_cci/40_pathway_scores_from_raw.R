#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(msigdbr)
  library(Seurat)
  library(UCell)
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
raw_root <- assert_production_path(
  file.path(root, "data", "public", "raw", "HCC-ST")
)
default_output <- file.path(
  root, "results", "50_spatial_cci", "pathways"
)
output_root <- normalizePath(
  parse_named_argument("output-root", default_output),
  mustWork = FALSE
)
section_argument <- parse_named_argument(
  "sections",
  paste(all_sections(), collapse = ",")
)
sections <- strsplit(section_argument, ",", fixed = TRUE)[[1L]]
invalid_sections <- setdiff(sections, all_sections())
if (length(invalid_sections) > 0L) {
  stop("Unknown sections: ", paste(invalid_sections, collapse = ", "))
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

msigdb <- msigdbr(species = "Homo sapiens")
gene_set <- function(name) {
  unique(msigdb$gene_symbol[msigdb$gs_name == name])
}
structural_noise <- grep(
  "^RP[SL]|^MT-|^MRP",
  unique(msigdb$gene_symbol),
  value = TRUE
)
signatures <- list(
  NFKB = gene_set("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  MAPK = c(
    "DUSP4", "DUSP5", "DUSP6", "ETV4", "ETV5",
    "SPRY2", "SPRY4", "PHLDA1", "SPRED1", "DUSP10"
  ),
  mTOR_signal = gene_set("HALLMARK_PI3K_AKT_MTOR_SIGNALING"),
  HYPOXIA = gene_set("HALLMARK_HYPOXIA"),
  FIBROSIS = c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
    "COL6A2", "COL6A3", "ACTA2", "TAGLN", "PDGFRB",
    "PDGFRA", "TIMP1", "LUM", "DCN", "POSTN", "SPARC",
    "BGN", "THBS2", "FBN1", "COL4A1"
  ),
  MAPK_KEGG = gene_set("KEGG_MAPK_SIGNALING_PATHWAY"),
  MAPK_ERK = gene_set("REACTOME_MAPK1_MAPK3_SIGNALING")
)
signatures <- lapply(
  signatures,
  function(genes) sort(unique(setdiff(genes, structural_noise)))
)
empty_signatures <- names(signatures)[lengths(signatures) == 0L]
if (length(empty_signatures) > 0L) {
  stop(
    "Current msigdbr lacks required terms: ",
    paste(empty_signatures, collapse = ", ")
  )
}

signature_table <- rbindlist(
  lapply(
    names(signatures),
    function(name) {
      data.table(
        program = name,
        gene = signatures[[name]],
        source = if (name %in% c("MAPK", "FIBROSIS")) {
          "fixed legacy gene list"
        } else {
          "msigdbr"
        }
      )
    }
  )
)
if ("db_version" %in% colnames(msigdb)) {
  signature_table[, msigdb_version := paste(unique(msigdb$db_version), collapse = ";")]
} else {
  signature_table[, msigdb_version := as.character(packageVersion("msigdbr"))]
}
fwrite(
  signature_table,
  file.path(output_root, "pathway_gene_sets.tsv"),
  sep = "\t"
)

score_section <- function(section) {
  section_index <- match(section, all_sections())
  if (is.na(section_index)) {
    stop("Unable to derive canonical seed index for ", section)
  }
  section_seed <- derive_r_seed(seed, section_index)
  set.seed(section_seed)
  section_dir <- assert_production_path(file.path(raw_root, section))
  spatial <- Load10X_Spatial(
    section_dir,
    filename = "filtered_feature_bc_matrix.h5"
  )
  spatial <- subset(spatial, subset = nCount_Spatial >= 200)
  spatial <- NormalizeData(spatial, verbose = FALSE)
  spatial <- AddModuleScore_UCell(
    spatial,
    features = signatures,
    name = "",
    ncores = 1L
  )
  metadata <- spatial[[]][, names(signatures), drop = FALSE]
  metadata$barcode <- rownames(metadata)
  metadata$section <- section
  metadata$seed <- section_seed
  metadata
}

scores <- bind_rows(
  lapply(sections, score_section)
)
saveRDS(scores, file.path(output_root, "pathway_scores_raw_visium.rds"))
fwrite(
  as.data.table(scores),
  file.path(output_root, "pathway_scores_raw_visium.csv.gz")
)
status <- as.data.table(scores)[
  ,
  .(n_scored_spots = .N, seed = unique(seed)),
  by = section
]
fwrite(status, file.path(output_root, "pathway_section_status.tsv"), sep = "\t")
cat(
  "Pathway scoring completed from raw Visium: ",
  nrow(scores),
  " spots across ",
  length(unique(scores$section)),
  " sections.\n",
  sep = ""
)
