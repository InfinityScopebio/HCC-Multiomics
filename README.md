# HCC Multi-omics: core computational workflows

Computation-only code snapshot assembled for **Reviewer 2, minor comment 10
(R2-m10)**. It covers the principal numerical workflows behind the HCC
multi-omics manuscript: CoVarNet, MuSiC, RCTD, COMMOT, CORNETO, MOFA, pathway
scoring, and custom cell-state/module projections.

This repository intentionally contains **no plotting, figure-rendering,
pixel-comparison, or page-assembly code**. Numerical results are written as
CSV, TSV, JSON, RDS, HDF5, or H5AD files. Raw/private matrices and frozen final
result objects are not published here. Animal-specific statistical workflows
and their dedicated inputs are intentionally excluded from this snapshot.

## Workflow map

| Module | Main computations | Entry point or scripts |
|---|---|---|
| Multi-omics preprocessing | Sample sheet and eight tissue/serum M0 matrices | `repro/scripts/10_run_raw_multiomics.sh` |
| scRNA atlas | Raw-count ingestion, RPCA/Harmony integration, subtype reconstruction | `repro/pipelines/20_scrna_cm02/03_*.R`, `04_*.R` |
| CoVarNet | Frequency filtering, K=5 nsNMF, CM02/network/membership | `repro/pipelines/20_scrna_cm02/05_rebuild_covarnet_cm02.R` |
| MuSiC and custom projections | Reference construction, four-cohort deconvolution, CM02 trajectory and survival | `repro/pipelines/20_scrna_cm02/27_*.R` through `31_*.R` |
| Hepatocyte activities | PROGENy MLM and CollecTRI ULM | `repro/pipelines/20_scrna_cm02/31_rebuild_hepc_activities.R` |
| Tissue multi-omics | PROGENy, PTM-SEA, RNA/protein concordance, KSEA, paired feature statistics | `repro/pipelines/30_tissue_multiomics/` |
| Drug response | DepMap/PRISM state signatures, TKI/chemotherapy composites, STORM limma/GSEA | `repro/pipelines/40_tissue_resistance/` |
| Spatial analysis | RCTD, CM recovery, radial profiles, pathway scores, COMMOT, LR-aware CCI | `repro/pipelines/50_spatial_cci/` |
| CORNETO | Curated MAPK prior network, expression-gated routes, solver repeatability | `repro/pipelines/50_spatial_cci/80_corneto_mapk.py` |
| Serum MOFA | Strict-monotonic input, K grid/final model training, TCGA axes and classifier | `repro/pipelines/60_serum_mofa/` |

## Requirements

- R 4.5 or newer; the validated article runtime was R 4.5.2.
- Python for general workflows and a separate Python environment containing
  scanpy, COMMOT and CORNETO for spatial analysis.
- Scientific package versions are inventoried in `repro/env/R_packages.tsv`,
  `repro/env/R_sessionInfo.txt`, `repro/env/python_base_freeze.txt`, and
  `repro/env/python_c2l_freeze.txt`.
- Inputs must be supplied locally under `repo/` as described in
  `docs/INPUTS.md`. No script silently downloads private data or substitutes a
  final result for a requested recomputation.

Runtime overrides:

```bash
export LIUF_RSCRIPT=/path/to/Rscript
export LIUF_PYTHON=/path/to/python3
export LIUF_C2L_PYTHON=/path/to/spatial/python
export MOFA_PYTHON=/path/to/mofa/python
```

## Run a module

List the available stages:

```bash
bash repro/run_core.sh --list
```

Examples:

```bash
# Curated-atlas CoVarNet rebuild
bash repro/run_core.sh covarnet

# MuSiC, CM02 trajectory, survival, and Hep_C pathway activities
bash repro/run_core.sh bulk_projection

# Tissue and resistance statistics
bash repro/run_core.sh tissue
bash repro/run_core.sh resistance

```

High-cost branches are explicit and are never launched by `--list`:

```bash
# Raw-count scRNA atlas reconstruction
bash repro/run_core.sh scrna_atlas

# Raw multi-omics preprocessing, required before tissue/MOFA stages
bash repro/run_core.sh preprocess

# Article-boundary spatial workflow
bash repro/run_core.sh spatial

# Diagnostic branch: refit all 16 RCTD sections before downstream spatial work
LIUF_SPATIAL_RUN_VARIANT=raw_rctd_diagnostic bash repro/run_core.sh spatial

# Retrain MOFA K=2..20 across six seeds and the final/replay K=15 models
MOFA_RUN_RAW_DIAGNOSTIC=1 MOFA_FORCE_RETRAIN=1 MOFA_WORKERS=4 \
  bash repro/run_core.sh serum_mofa
```

The default MOFA downstream branch uses the declared historical K=15 article
model. Set `MOFA_ANALYSIS_SOURCE=raw_derived` after raw retraining to compute
downstream results from the new model. These alternatives are diagnostic and
must not be described as byte-identical article replacement without separate
numeric equivalence checks.

## Random seeds and execution contract

The project master seed is `20260807`, exported as `LIUF_MASTER_SEED` and
`PYTHONHASHSEED`; BLAS/OpenMP thread counts default to one. Historical
step-specific seeds are preserved when they affect accepted results.

| Analysis | Seed contract |
|---|---|
| CoVarNet K=5 nsNMF | 77, sequential backend |
| Bulk CM02 UMAP | 42 |
| Hep_C PROGENy/CollecTRI | 42 |
| PTM-SEA | Historical 42 semantics; explicit serial repeat probe |
| STORM GSEA | 1, `BiocParallel::SerialParam` |
| COMMOT | Master seed plus canonical section index |
| CORNETO/HiGHS | Solver seed 0, one thread, two isolated complete runs |
| MOFA | Model seed 1; K-scan seeds 1–6 |

Additional per-stage seeds and RNG kinds are recorded in the scripts and the
generated `repo/audit` tables. A seed declaration alone is not presented as
proof of repeatability.

## Article-equivalence boundaries

Exact article-equivalent analysis is not synonymous with an exclusive raw-only
recomputation. The following boundaries are explicit:

- Expert-curated Seurat atlas and subtype objects are checksum-locked inputs;
  raw atlas reconstruction is a separate stage.
- The article spatial workflow validates 18 historical RCTD/CoVarNet objects
  before downstream calculations. Raw RCTD refitting is isolated as a
  diagnostic branch.
- Historical positive-flow CORNETO topology and current signal-state solutions
  are not treated as interchangeable. The current solver branch recomputes its
  own routes twice.
- Historical MOFA input/model objects support article-boundary downstream
  calculations; the raw branch retrains rather than copying a frozen model.
- Validation-only targets used for numeric comparisons are never accepted as
  silent substitutes for production computation.

See `docs/BOUNDARIES.md` for the expected local boundary files and their
manifests.

## Validation

```bash
python3 -m unittest discover -s repro/tests -p 'test_*.py' -v
python3 -m unittest repro.tests.test_public_core -v

# Parse every R script without executing analyses
Rscript -e 'files <- list.files("repro", pattern="[.]R$", recursive=TRUE, full.names=TRUE); for (f in files) parse(f); cat(length(files), "R scripts parsed\n")'

# Check shell entry points
find repro -name '*.sh' -print0 | xargs -0 -n1 bash -n
```

The public-core test checks method coverage, absence of plotting calls,
absence of binary figure/raw-data assets in the code snapshot, Python syntax,
and the SHA-256 values in `SOURCE_MANIFEST.tsv`.

## Provenance

`SOURCE_MANIFEST.tsv` maps every copied source file to its original and public
SHA-256. Files marked `computation-only cleanup` have plotting segments removed
while retaining numerical outputs and execution/seed audits. The repository
contains no manuscript DOCX, figure PDF, slide deck, patient workbook, raw
omics matrix, or private project history.
