## Phase M0.3-M0.5 — load all 8 matrices, ID standardize, missingness QC,
## save canonical RDS bundle.
##
## Output:
##   multiomics/dat/serum_{lcms,gcms,protein,phospho}.rds
##   multiomics/dat/tissue_{rna,protein,phospho,metabolite}.rds
##   multiomics/qc/missval_summary.csv
##   multiomics/qc/id_mapping_report.csv
##
## Each RDS = list(meta = feature-meta data.table, mat = numeric matrix
##                 [feature × sample], feature_id = char, na_pct = numeric).
## Sample columns subset to canonical sample_sheet.

suppressMessages({
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)
})
dir.create("multiomics/dat", recursive = TRUE, showWarnings = FALSE)
dir.create("multiomics/qc",  recursive = TRUE, showWarnings = FALSE)

ss <- fread("multiomics/dat/sample_sheet.csv")
serum_ids  <- ss[source == "serum",  sample_id]
tissue_ids <- ss[source == "tissue", sample_id]

map_uniprot_to_symbol <- function(acc_vec) {
  acc_clean <- sub(";.*","", acc_vec)        # drop multi-acc remainder
  acc_clean <- sub("-\\d+$","", acc_clean)   # drop isoform suffix
  ans <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db, keys = unique(acc_clean),
      columns = c("SYMBOL","ENSEMBL"), keytype = "UNIPROT")),
    error = function(e) data.frame(UNIPROT = unique(acc_clean),
                                   SYMBOL = NA, ENSEMBL = NA))
  ans <- as.data.table(ans)
  ans <- ans[, .SD[1], by = UNIPROT]
  ans[match(acc_clean, UNIPROT), .(SYMBOL, ENSEMBL)]
}

map_symbol_to_uni <- function(sym) {
  ans <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db, keys = unique(sym),
      columns = c("UNIPROT","ENSEMBL"), keytype = "SYMBOL")),
    error = function(e) data.frame(SYMBOL = unique(sym),
                                   UNIPROT = NA, ENSEMBL = NA))
  ans <- as.data.table(ans)[, .SD[1], by = SYMBOL]
  ans[match(sym, SYMBOL), .(UNIPROT, ENSEMBL)]
}

id_report <- list()

## ============== SERUM LCMS ==============
cat("\n[SERUM LCMS]\n")
d <- fread("LiuF_data/serum/代谢_LCMS数据矩阵.csv")
meta_cols <- 1:21
mat_cols  <- intersect(serum_ids, colnames(d))
mat <- as.matrix(d[, ..mat_cols])
rownames(mat) <- d$ID
meta <- d[, ..meta_cols]
na_pct <- rowMeans(is.na(mat))
cat("dim", paste(dim(mat),collapse="x"), " sample matched", length(mat_cols),
    " na median", round(median(na_pct),3),"\n")
id_report$serum_lcms <- list(n_feat = nrow(mat),
                             with_kegg = sum(!is.na(meta$kegg) & meta$kegg != ""),
                             with_hmdb = sum(!is.na(meta$HMDB) & meta$HMDB != ""))
saveRDS(list(meta = meta, mat = mat, feature_id = d$ID, na_pct = na_pct,
             layer = "serum_lcms"),
        "multiomics/dat/serum_lcms.rds")

## ============== SERUM GCMS ==============
cat("\n[SERUM GCMS]\n")
d <- fread("LiuF_data/serum/代谢_GCMS数据矩阵.csv")
meta_cols <- 1:16
mat_cols  <- intersect(serum_ids, colnames(d))
mat <- as.matrix(d[, ..mat_cols])
rownames(mat) <- as.character(d$`Alignment ID`)
meta <- d[, ..meta_cols]
na_pct <- rowMeans(is.na(mat))
cat("dim", paste(dim(mat),collapse="x"), " sample matched", length(mat_cols),
    " na median", round(median(na_pct),3),"\n")
id_report$serum_gcms <- list(n_feat = nrow(mat),
                             with_kegg = sum(!is.na(meta$kegg) & meta$kegg != ""),
                             with_hmdb = sum(!is.na(meta$HMDB) & meta$HMDB != ""))
saveRDS(list(meta = meta, mat = mat,
             feature_id = as.character(d$`Alignment ID`),
             na_pct = na_pct, layer = "serum_gcms"),
        "multiomics/dat/serum_gcms.rds")

## ============== SERUM PROTEIN (peptide-level → protein roll-up) ==============
cat("\n[SERUM PROTEIN]\n")
d <- fread("LiuF_data/serum/DIA_Peptide groups.csv")
meta_cols <- 1:6
mat_cols  <- intersect(serum_ids, colnames(d))
## Roll up peptide → ProteinGroups by sum (TMT/DIA standard)
prot_map  <- map_uniprot_to_symbol(d$ProteinGroups)
d[, SYMBOL := prot_map$SYMBOL][, ENSEMBL := prot_map$ENSEMBL]
mat_pep <- as.matrix(d[, ..mat_cols])
mat_pep[is.na(mat_pep)] <- 0
agg <- rowsum(mat_pep, group = d$ProteinGroups, na.rm = TRUE)
agg[agg == 0] <- NA
prot_meta <- unique(d[, .(ProteinGroups, SYMBOL, ENSEMBL)])
prot_meta <- prot_meta[match(rownames(agg), ProteinGroups)]
na_pct <- rowMeans(is.na(agg))
cat("peptides", nrow(d), " → proteins", nrow(agg),
    " sample matched", length(mat_cols),
    " na median", round(median(na_pct),3), "\n")
cat("UniProt → SYMBOL mapped:", sum(!is.na(prot_meta$SYMBOL)),"/", nrow(prot_meta),"\n")
id_report$serum_prot <- list(n_pep = nrow(d), n_prot = nrow(agg),
                             sym_map = sum(!is.na(prot_meta$SYMBOL)))
saveRDS(list(meta = prot_meta, mat = agg,
             feature_id = rownames(agg), na_pct = na_pct,
             layer = "serum_protein"),
        "multiomics/dat/serum_protein.rds")

## ============== SERUM PHOSPHO ==============
cat("\n[SERUM PHOSPHO]\n")
d <- fread("LiuF_data/serum/DIAP_Phospho(STY)Sites.csv")
mat_cols <- intersect(serum_ids, colnames(d))
mat <- as.matrix(d[, ..mat_cols])
site_id <- paste(d$ProteinAccessions, d$SiteAA, d$SiteLocation, sep = "_")
rownames(mat) <- site_id
prot_map <- map_uniprot_to_symbol(d$ProteinAccessions)
meta <- data.table(site_id, ProteinAccession = d$ProteinAccessions,
                   Gene = d$`Gene Name`, SYMBOL = prot_map$SYMBOL,
                   ENSEMBL = prot_map$ENSEMBL,
                   SiteAA = d$SiteAA, SiteLocation = d$SiteLocation,
                   FlankingRegion = d$FlankingRegion)
na_pct <- rowMeans(is.na(mat))
cat("sites", nrow(mat), " sample matched", length(mat_cols),
    " na median", round(median(na_pct),3), "\n")
cat("UniProt → SYMBOL mapped:", sum(!is.na(meta$SYMBOL)),"/", nrow(meta),"\n")
id_report$serum_phos <- list(n_site = nrow(mat),
                             sym_map = sum(!is.na(meta$SYMBOL)),
                             na_median = median(na_pct))
saveRDS(list(meta = meta, mat = mat, feature_id = site_id,
             na_pct = na_pct, layer = "serum_phospho"),
        "multiomics/dat/serum_phospho.rds")

## ============== TISSUE RNA ==============
cat("\n[TISSUE RNA]\n")
d <- fread("LiuF_data/tissue/RNAseq_counts_anno.xls")
mat_cols <- intersect(tissue_ids, colnames(d))
mat <- as.matrix(d[, ..mat_cols])
rownames(mat) <- d$id
meta_cols <- setdiff(colnames(d), c(mat_cols, "id"))
meta <- d[, c("id", meta_cols), with = FALSE]
## Map id (looks like SYMBOL) → ENSEMBL
sym_map <- map_symbol_to_uni(d$id)
meta[, SYMBOL := id]
meta[, ENSEMBL := sym_map$ENSEMBL]
meta[, UNIPROT := sym_map$UNIPROT]
na_pct <- rowMeans(mat == 0)
cat("genes", nrow(mat), " sample matched", length(mat_cols),
    " genes with all-zero", sum(na_pct == 1), "\n")
cat("SYMBOL → ENSEMBL mapped:", sum(!is.na(meta$ENSEMBL)),"/", nrow(meta),"\n")
id_report$tissue_rna <- list(n_gene = nrow(mat),
                             ens_map = sum(!is.na(meta$ENSEMBL)),
                             zero_genes = sum(na_pct == 1))
saveRDS(list(meta = meta, mat = mat, feature_id = d$id,
             na_pct = na_pct, layer = "tissue_rna"),
        "multiomics/dat/tissue_rna.rds")

## ============== TISSUE PROTEIN ==============
cat("\n[TISSUE PROTEIN]\n")
d <- fread("LiuF_data/tissue/DIA_Protein quantitation.csv")
mat_cols <- intersect(tissue_ids, colnames(d))
mat <- as.matrix(d[, ..mat_cols])
prot_map <- map_uniprot_to_symbol(d$Accession)
rownames(mat) <- d$Accession
meta <- data.table(Accession = d$Accession, Gene = d$`Gene Name`,
                   SYMBOL = prot_map$SYMBOL, ENSEMBL = prot_map$ENSEMBL,
                   Description = d$Description, MW_kDa = d$`MW [kDa]`,
                   Qvalue = d$Qvalue, Cscore = d$Cscore)
na_pct <- rowMeans(is.na(mat))
cat("proteins", nrow(mat), " sample matched", length(mat_cols),
    " na median", round(median(na_pct),3),"\n")
cat("UniProt → SYMBOL mapped:", sum(!is.na(meta$SYMBOL)),"/", nrow(meta),"\n")
id_report$tissue_prot <- list(n_prot = nrow(mat),
                              sym_map = sum(!is.na(meta$SYMBOL)),
                              na_median = median(na_pct))
saveRDS(list(meta = meta, mat = mat, feature_id = d$Accession,
             na_pct = na_pct, layer = "tissue_protein"),
        "multiomics/dat/tissue_protein.rds")

## ============== TISSUE PHOSPHO ==============
cat("\n[TISSUE PHOSPHO]\n")
d <- fread("LiuF_data/tissue/DIAP_Phospho(STY)Sites.csv")
mat_cols <- intersect(tissue_ids, colnames(d))
mat <- as.matrix(d[, ..mat_cols])
site_id <- paste(d$ProteinAccessions, d$SiteAA, d$SiteLocation, sep = "_")
rownames(mat) <- site_id
prot_map <- map_uniprot_to_symbol(d$ProteinAccessions)
meta <- data.table(site_id, ProteinAccession = d$ProteinAccessions,
                   Gene = d$`Gene Name`, SYMBOL = prot_map$SYMBOL,
                   ENSEMBL = prot_map$ENSEMBL,
                   SiteAA = d$SiteAA, SiteLocation = d$SiteLocation,
                   FlankingRegion = d$FlankingRegion)
na_pct <- rowMeans(is.na(mat))
cat("sites", nrow(mat), " sample matched", length(mat_cols),
    " na median", round(median(na_pct),3),"\n")
cat("UniProt → SYMBOL mapped:", sum(!is.na(meta$SYMBOL)),"/", nrow(meta),"\n")
id_report$tissue_phos <- list(n_site = nrow(mat),
                              sym_map = sum(!is.na(meta$SYMBOL)),
                              na_median = median(na_pct))
saveRDS(list(meta = meta, mat = mat, feature_id = site_id,
             na_pct = na_pct, layer = "tissue_phospho"),
        "multiomics/dat/tissue_phospho.rds")

## ============== TISSUE METABOLITE ==============
cat("\n[TISSUE METABOLITE]\n")
d <- fread("LiuF_data/tissue/代谢.csv")
meta_cols <- 1:25
mat_cols  <- intersect(tissue_ids, colnames(d))
mat <- as.matrix(d[, ..mat_cols])
rownames(mat) <- d$ID
meta <- d[, ..meta_cols]
na_pct <- rowMeans(is.na(mat))
cat("dim", paste(dim(mat),collapse="x"), " sample matched", length(mat_cols),
    " na median", round(median(na_pct),3),"\n")
id_report$tissue_meta <- list(n_feat = nrow(mat),
                              with_kegg = sum(!is.na(meta$KEGG) & meta$KEGG != ""),
                              with_hmdb = sum(!is.na(meta$HMDB) & meta$HMDB != ""))
saveRDS(list(meta = meta, mat = mat, feature_id = d$ID,
             na_pct = na_pct, layer = "tissue_metabolite"),
        "multiomics/dat/tissue_metabolite.rds")

## ============== Missingness summary ==============
ms_files <- list.files("multiomics/dat", "\\.rds$", full.names = TRUE)
ms <- lapply(ms_files, function(f) {
  o <- readRDS(f)
  data.table(layer = o$layer, n_feat = nrow(o$mat), n_sample = ncol(o$mat),
             na_median = round(median(o$na_pct), 3),
             na_q90    = round(quantile(o$na_pct, 0.9), 3),
             feat_lt30pct_NA = sum(o$na_pct < 0.3),
             feat_gt70pct_NA = sum(o$na_pct > 0.7))
})
ms <- rbindlist(ms)
fwrite(ms, "multiomics/qc/missval_summary.csv")
cat("\n=== missingness summary ===\n"); print(ms)

## ============== ID mapping summary ==============
ir <- rbindlist(lapply(names(id_report), function(n) {
  x <- id_report[[n]]
  # data.table >= 1.18 reserves `key` as a constructor argument.
  data.table(layer = n, metric = names(x), value = unlist(x))
}))
fwrite(ir, "multiomics/qc/id_mapping_report.csv")
cat("\n=== id mapping ===\n"); print(ir)

cat("\n[", format(Sys.time()), "] M0 done\n")
