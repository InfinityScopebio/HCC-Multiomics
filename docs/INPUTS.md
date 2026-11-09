# Local input contract

The repository is code-first. Place inputs at the following relative paths
before running the corresponding stages. Missing files must produce an explicit
error; do not replace them with final result tables.

## Internal scRNA-seq

Six 10x filtered count directories are required for the optional raw atlas
stage. `CC1–CC3` are prohibited and are not part of the analysis set.

```text
data/private/raw/data/scRNAData/{HCC1,HCC2,HCC3,NC1,NC2,NC3}/filtered_feature_bc_matrix/
  matrix.mtx.gz
  features.tsv.gz
  barcodes.tsv.gz
```

Public scRNA inputs:

```text
data/public/raw/GSE115469/GSE115469_Data.csv.gz
data/public/raw/GSE136103/GSE136103_extracted/*_matrix.mtx.gz
data/public/raw/GSE136103/GSE136103_extracted/*_genes.tsv.gz
data/public/raw/GSE136103/GSE136103_extracted/*_barcodes.tsv.gz
data/public/raw/GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz
data/public/raw/GSE149614/GSE149614_HCC.metadata.updated.txt.gz
```

The GSE136103 directory contains 26 files of each type. Accessions are
GSE115469, GSE136103, and GSE149614.

## Curated single-cell boundary

The default CoVarNet, projection, and spatial profile calculations use these
checksummed expert-curated objects:

```text
data/curation/scRNA/combined_annotated_v5.rds
data/curation/scRNA/hep_annotated.rds
data/curation/scRNA/t_annotated.rds
data/curation/scRNA/mp_annotated.rds
data/curation/scRNA/fib_annotated.rds
data/curation/scRNA/b_annotated.rds
```

`config/manifests/scrna_frozen_curation.tsv` records their expected dimensions
and SHA-256 values. These are annotation/curation boundaries, not raw matrices.

## Tissue and serum matrices

```text
data/private/raw/data/tissue/RNAseq_counts_anno.xls
data/private/raw/data/tissue/DIA_Protein quantitation.csv
data/private/raw/data/tissue/DIAP_Phospho(STY)Sites.csv
data/private/raw/data/tissue/代谢.csv
data/private/raw/data/serum/DIA_Peptide groups.csv
data/private/raw/data/serum/DIAP_Phospho(STY)Sites.csv
data/private/raw/data/serum/代谢_GCMS数据矩阵.csv
data/private/raw/data/serum/代谢_LCMS数据矩阵.csv
data/private/raw/data/样本信息.xlsx
```

## Bulk and drug-response cohorts

```text
data/public/raw/TCGA-LIHC/LIHC_HiSeqV2.gz
data/public/raw/TCGA-LIHC/LIHC_clinicalMatrix.tsv
data/public/raw/GSE63898/GSE63898_series_matrix.txt.gz
data/public/raw/GSE63898/GPL13667.soft.gz
data/public/raw/GTEx-Liver/GTEx_LIVER_recount3.rds
data/public/raw/GTEx-Liver/GTEx_LIVER_recount3_rse.rds
data/public/raw/DepMap/Model.csv
data/public/raw/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
data/public/raw/DepMap/Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv
data/public/raw/GSE109211/GSE109211_series_matrix.txt.gz
data/public/raw/GSE109211/GPL13938_GEO_platform_full.txt
```

TCGA uses the locally retained UCSC Xena legacy HiSeqV2 matrix. Other
accessions are GSE63898, GTEx Liver, DepMap/PRISM, and BIOSTORM GSE109211.
Release-specific files should be retained rather than replaced with a current
release having different identifiers.

## Spatial transcriptomics

Sixteen Space Ranger sample directories are expected:

```text
data/public/raw/HCC-ST/*/filtered_feature_bc_matrix.h5
data/public/raw/HCC-ST/*/spatial/tissue_positions_list.csv
data/public/raw/HCC-ST/*/spatial/scalefactors_json.json
data/public/raw/HCC-ST/*/spatial/tissue_hires_image.png
```

The raw-RCTD branch additionally requires the curated scRNA reference above.
The tissue image is used for input auditing, not by this computation-only
snapshot for rendering.

## Knowledge resources

```text
data/reference/MSigDB/msigdb_v2023.2.Hs_files_to_download_locally.zip
data/reference/PROGENy/progeny_net.rds
data/reference/CollecTRI/collectri_net.csv
data/reference/PTMsigDB/ptm.sig.db.all.flanking.human.v2.0.0.gmt
data/reference/PTMsigDB/ptm.sig.db.all.uniprot.human.v2.0.0.gmt
data/reference/OmniPath/omnipath_kinase_substrate.rds
```

The KEGG compound-to-pathway table is tracked at
`data/reference/KEGG/kegg_cpd_pathway.tsv`. COMMOT's bundled human
CellChat ligand/receptor table is read from the installed `commot` package and
its hash is recorded.

## Small inputs included with the code

Only fixed annotation/curation or knowledge inputs are included:

- cluster-to-state annotation CSVs and raw-derived overlay maps;
- `data/curation/resistance/HepC_signatures_main.gmt`;
- the KEGG compound-to-pathway TSV.

No patient-level clinical workbook, raw expression matrix, Seurat RDS, RCTD
weight, MOFA model, COMMOT cache, or final result table is included.
