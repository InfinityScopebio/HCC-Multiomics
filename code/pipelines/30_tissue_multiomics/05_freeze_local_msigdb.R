#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_path <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
source(file.path(dirname(normalizePath(script_path)), "common.R"))

cache_path <- Sys.getenv(
  "MSIGDB_LOCAL_CACHE",
  unset = path.expand(
    "~/Library/Caches/org.R-project.R/R/msigdbr/msigdb.2025.1.Hs.rds"
  )
)
assert_files(cache_path)
database <- as.data.table(readRDS(cache_path))
required_columns <- c(
  "db_gene_symbol", "gs_id", "gs_name", "gs_collection",
  "gs_subcollection", "db_version"
)
stopifnot(all(required_columns %in% colnames(database)))

selected <- database[
  gs_collection == "H" |
    (
      gs_collection == "C2" &
        gs_subcollection %in% c("CP:REACTOME", "CP:KEGG_LEGACY")
    ),
  .(
    gs_id,
    gs_name,
    gs_collection,
    gs_subcollection,
    gene_symbol = db_gene_symbol,
    db_version
  )
]
selected <- unique(selected)
stopifnot(
  nrow(selected) > 0L,
  uniqueN(selected$db_version) == 1L,
  all(c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_HYPOXIA",
    "REACTOME_MAPK1_MAPK3_SIGNALING"
  ) %in% selected$gs_name)
)

output <- result_path("references", "msigdb_tissue_sets.rds")
saveRDS(selected, output, version = 3)
fwrite(
  data.table(
    source_path = cache_path,
    source_md5 = unname(tools::md5sum(cache_path)),
    db_version = unique(selected$db_version),
    rows = nrow(selected),
    gene_sets = uniqueN(selected$gs_id),
    collections = paste(sort(unique(selected$gs_collection)), collapse = ";"),
    subcollections = paste(
      sort(unique(na.omit(selected$gs_subcollection))),
      collapse = ";"
    ),
    frozen_path = output,
    frozen_md5 = unname(tools::md5sum(output))
  ),
  result_path("references", "msigdb_tissue_sets_manifest.tsv"),
  sep = "\t"
)

message(
  "Frozen local MSigDB ", unique(selected$db_version), ": ",
  uniqueN(selected$gs_id), " gene sets."
)
