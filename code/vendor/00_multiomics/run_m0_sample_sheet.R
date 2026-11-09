## Phase M0.2 — canonical sample sheet
## Output: multiomics/dat/sample_sheet.csv  (one row per physical sample,
##         + presence flags for each of 8 omics layers)

suppressMessages({library(data.table); library(readxl)})
dir.create("multiomics/dat", recursive = TRUE, showWarnings = FALSE)
dir.create("multiomics/qc",  recursive = TRUE, showWarnings = FALSE)

## --- Pull sample columns from each matrix header ---
read_cols <- function(path, skip_meta) {
  d <- fread(path, nrows = 1)
  setdiff(colnames(d)[(skip_meta + 1):ncol(d)], "")
}
serum_lcms <- read_cols("LiuF_data/serum/代谢_LCMS数据矩阵.csv", 21)
serum_gcms <- read_cols("LiuF_data/serum/代谢_GCMS数据矩阵.csv", 16)
serum_prot <- read_cols("LiuF_data/serum/DIA_Peptide groups.csv", 6)
serum_phos <- read_cols("LiuF_data/serum/DIAP_Phospho(STY)Sites.csv", 14)
tis_rna    <- {
  d <- fread("LiuF_data/tissue/RNAseq_counts_anno.xls", nrows = 1)
  cn <- colnames(d); cn[grepl("^(Can|Para)\\d+$", cn)]
}
tis_prot   <- {
  d <- fread("LiuF_data/tissue/DIA_Protein quantitation.csv", nrows = 1)
  cn <- colnames(d); cn[grepl("^(Can|Para)\\d+$", cn)]
}
tis_phos   <- {
  d <- fread("LiuF_data/tissue/DIAP_Phospho(STY)Sites.csv", nrows = 1)
  cn <- colnames(d); cn[grepl("^(Can|Para)\\d+$", cn)]
}
tis_meta   <- {
  d <- fread("LiuF_data/tissue/代谢.csv", nrows = 1)
  cn <- colnames(d); cn[grepl("^(Can|Para)\\d+$", cn)]
}

## --- Strip QCs ---
strip_qc <- function(x) x[!grepl("^QC", x)]
serum_lcms <- strip_qc(serum_lcms)
serum_gcms <- strip_qc(serum_gcms)

## --- Build canonical sample sheet ---
all_serum  <- sort(unique(c(serum_lcms, serum_gcms, serum_prot, serum_phos)))
all_tissue <- sort(unique(c(tis_rna, tis_prot, tis_phos, tis_meta)))

assign_group <- function(x) {
  fcase(
    grepl("^H\\d+$", x),     "Healthy",
    grepl("^HF\\d+$", x),    "HepFibrosis",
    grepl("^Cir\\d+$", x),   "Cirrhosis",
    grepl("^Lcan\\d+$", x),  "HCC_serum",
    grepl("^Can\\d+$", x),   "HCC_tumor",
    grepl("^Para\\d+$", x),  "Adjacent",
    default = NA_character_
  )
}
patient_id <- function(x) {
  sub("^(H|HF|Cir|Lcan|Can|Para)", "", x)
}

ss <- data.table(
  sample_id = c(all_serum, all_tissue),
  source    = c(rep("serum",  length(all_serum)),
                rep("tissue", length(all_tissue)))
)
ss[, group       := assign_group(sample_id)]
ss[, patient_num := patient_id(sample_id)]
ss[, paired_id   := fcase(
        group == "HCC_serum",  paste0("LiuF_HCC_", patient_num),
        group == "HCC_tumor",  paste0("LiuF_HCC_", patient_num),
        group == "Adjacent",   paste0("LiuF_HCC_", patient_num),
        default = NA_character_)]
ss[, group_4 := fcase(
        group == "Healthy",     "Healthy",
        group == "HepFibrosis", "HepFibrosis",
        group == "Cirrhosis",   "Cirrhosis",
        group %in% c("HCC_serum","HCC_tumor"), "HCC",
        group == "Adjacent",    "Adjacent")]

## --- Presence flags per layer ---
ss[, serum_lcms := sample_id %in% serum_lcms]
ss[, serum_gcms := sample_id %in% serum_gcms]
ss[, serum_prot := sample_id %in% serum_prot]
ss[, serum_phos := sample_id %in% serum_phos]
ss[, tissue_rna  := sample_id %in% tis_rna]
ss[, tissue_prot := sample_id %in% tis_prot]
ss[, tissue_phos := sample_id %in% tis_phos]
ss[, tissue_meta := sample_id %in% tis_meta]
ss[, n_layers := serum_lcms + serum_gcms + serum_prot + serum_phos +
                 tissue_rna + tissue_prot + tissue_phos + tissue_meta]

setcolorder(ss, c("sample_id","source","group","group_4","patient_num",
                  "paired_id","n_layers"))

fwrite(ss, "multiomics/dat/sample_sheet.csv")

## --- Summaries ---
cat("\n=== group × source ===\n")
print(dcast(ss, group ~ source, value.var = "sample_id", fun.aggregate = length))

cat("\n=== layer coverage (serum) ===\n")
serum_only <- ss[source == "serum"]
print(serum_only[, .(n = .N,
                     lcms = sum(serum_lcms), gcms = sum(serum_gcms),
                     prot = sum(serum_prot), phos = sum(serum_phos)),
                 by = group_4])

cat("\n=== layer coverage (tissue) ===\n")
tis_only <- ss[source == "tissue"]
print(tis_only[, .(n = .N,
                   rna = sum(tissue_rna), prot = sum(tissue_prot),
                   phos = sum(tissue_phos), meta = sum(tissue_meta)),
               by = group_4])

cat("\n=== tissue paired check ===\n")
paired_check <- tis_only[, .N, by = paired_id][order(paired_id)]
cat("paired_id with both Can/Para:", sum(paired_check$N == 2),
    "/ unique paired:", nrow(paired_check), "\n")
miss_paired <- paired_check[N != 2]
if (nrow(miss_paired)) print(miss_paired)

cat("\n=== Lcan ↔ Can patient_num cross-source check ===\n")
lc <- ss[group == "HCC_serum", patient_num]
cn <- ss[group == "HCC_tumor", patient_num]
cat("serum HCC patient_num:", paste(sort(as.integer(lc)), collapse=","), "\n")
cat("tissue HCC patient_num:", paste(sort(as.integer(cn)), collapse=","), "\n")
cat("intersect:", length(intersect(lc, cn)), "/ union:", length(union(lc, cn)), "\n")

cat("\n=== n samples per layer ===\n")
print(c(lcms=length(serum_lcms), gcms=length(serum_gcms),
        s_prot=length(serum_prot), s_phos=length(serum_phos),
        rna=length(tis_rna), t_prot=length(tis_prot),
        t_phos=length(tis_phos), t_meta=length(tis_meta)))

cat("\n=== output ===\n")
cat("multiomics/dat/sample_sheet.csv  (rows =", nrow(ss), ")\n")
