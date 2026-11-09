#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
source(file.path(
  dirname(sub("^--file=", "", grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )[[1L]])),
  "common.R"
))

required <- c(
  data_path("sample_sheet.csv"),
  data_path("tissue_rna.rds"),
  data_path("tissue_protein.rds"),
  data_path("tissue_phospho.rds"),
  data_path("tissue_metabolite.rds"),
  reference_path("PROGENy", "progeny_net.rds"),
  reference_path(
    "PTMsigDB", "ptm.sig.db.all.uniprot.human.v2.0.0.gmt"
  ),
  reference_path("OmniPath", "omnipath_kinase_substrate.rds")
)
assert_files(required)

msigdb_cache <- Sys.getenv(
  "MSIGDB_LOCAL_CACHE",
  unset = path.expand(
    "~/Library/Caches/org.R-project.R/R/msigdbr/msigdb.2025.1.Hs.rds"
  )
)
assert_files(msigdb_cache)

registry_path <- file.path(
  project_root, "repro", "manifests", "panel_registry.tsv"
)
assert_files(registry_path)
registry <- fread(registry_path)
scope <- registry[
  (collection == "main" & panel %in% c("4A", "4B", "4C")) |
    (
      collection == "supplemental" &
        panel %in% c("S3A", "S3B", "S3C")
    )
]
setorder(scope, collection, panel)
stopifnot(
  nrow(scope) == 6L,
  all(scope$stage == "40_tissue_resistance")
)
fwrite(
  scope,
  result_path("panel_registry_scope.tsv"),
  sep = "\t"
)

hashes <- data.table(
  role = c(
    "sample_sheet", "tissue_rna", "tissue_protein", "tissue_phospho",
    "tissue_metabolite", "progeny", "ptmsigdb", "omnipath"
  ),
  path = required,
  md5 = unname(tools::md5sum(required)),
  bytes = file.info(required)$size
)
hashes <- rbind(
  hashes,
  data.table(
    role = "msigdb_local_cache",
    path = msigdb_cache,
    md5 = unname(tools::md5sum(msigdb_cache)),
    bytes = file.info(msigdb_cache)$size
  )
)
fwrite(hashes, file.path(audit_root, "tissue_input_hashes.tsv"), sep = "\t")

metadata <- paired_tissue_metadata()
stopifnot(length(metadata$patients) == 17L)
stopifnot(!any(grepl("^CC[1-3]$", metadata$tissue$sample_id)))

objects <- lapply(required[2:5], readRDS)
names(objects) <- c("rna", "protein", "phospho", "metabolite")
dimensions <- rbindlist(lapply(names(objects), function(layer) {
  object <- objects[[layer]]
  data.table(
    layer = layer,
    features = nrow(object$mat),
    samples = ncol(object$mat),
    hcc_samples = sum(grepl("^Can", colnames(object$mat))),
    adjacent_samples = sum(grepl("^Para", colnames(object$mat)))
  )
}))
stopifnot(all(dimensions$samples == 34L))
stopifnot(all(dimensions$hcc_samples == 17L))
stopifnot(all(dimensions$adjacent_samples == 17L))
fwrite(
  dimensions,
  file.path(audit_root, "tissue_input_dimensions.tsv"),
  sep = "\t"
)

message("Validated local tissue inputs and frozen references.")
