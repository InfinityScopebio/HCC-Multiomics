# HCC Multi-omics Analysis Workflows

Reproducible R and Python workflows for hepatocellular carcinoma multi-omics
research. The repository covers CoVarNet, MuSiC, RCTD, COMMOT, CORNETO, MOFA,
pathway scoring, and custom cell-state/module projections.

## Repository layout

| Path | Contents |
|---|---|
| `code/` | Pipelines, fixed upstream sources, shared libraries, utilities, tests, and environment inventories |
| `config/` | Project configuration and static input/boundary manifests |
| `data/` | Small annotation/curation inputs and the local data root |
| `docs/` | Input contracts and analysis boundaries |
| `run_core.sh` | Unified stage runner |
| `SOURCE_MANIFEST.tsv` | Source/public paths, hashes, sizes, and transformations |

The main subdivisions under `code/` are `pipelines/`, `vendor/`, `lib/`,
`scripts/`, `tests/`, and `env/`. Local runs write numerical outputs to
`results/`, execution records to `audit/`, logs to `logs/`, and generated
inventories to `manifests/`. These output directories are excluded from version
control.

## Analysis modules

| Stage | Main computations |
|---|---|
| `preprocess` | Sample sheet and eight tissue/serum M0 matrices |
| `scrna_atlas` | Raw-count ingestion, RPCA/Harmony integration, subtype reconstruction |
| `covarnet` | Frequency filtering, K=5 nsNMF, CM02 network and membership |
| `bulk_projection` | MuSiC, CM02 trajectory, survival, and Hep_C activities |
| `atlas_tables` | Composition and ratio-of-enrichment tables |
| `tissue` | PROGENy, PTM-SEA, cross-layer concordance, KSEA, paired feature statistics |
| `ptmsea_rng` | PTM-SEA seed-mode comparison |
| `resistance` | DepMap/PRISM composites and STORM limma/GSEA |
| `spatial` | RCTD, CM recovery, radial profiles, pathways, COMMOT/CCI, CORNETO |
| `serum_mofa` | MOFA input/model training, TCGA axes, survival, and classifier |

## Requirements

- R 4.5 or newer; the validated article runtime was R 4.5.2.
- Python for general workflows and a separate environment containing scanpy,
  COMMOT, and CORNETO for spatial analysis.
- Local inputs arranged as described in `docs/INPUTS.md` and
  `docs/BOUNDARIES.md`.

Package versions are recorded in:

- `code/env/R_packages.tsv` and `code/env/R_sessionInfo.txt`
- `code/env/python_base_freeze.txt`
- `code/env/python_c2l_freeze.txt`

Runtime overrides:

```bash
export LIUF_RSCRIPT=/path/to/Rscript
export LIUF_PYTHON=/path/to/python3
export LIUF_C2L_PYTHON=/path/to/spatial/python
export MOFA_PYTHON=/path/to/mofa/python
```

## Run

List stages without starting a computation:

```bash
bash run_core.sh --list
```

Examples:

```bash
bash run_core.sh covarnet
bash run_core.sh bulk_projection
bash run_core.sh tissue
bash run_core.sh resistance
```

Higher-cost branches are explicit:

```bash
# Raw-count scRNA atlas reconstruction
bash run_core.sh scrna_atlas

# Raw multi-omics preprocessing
bash run_core.sh preprocess

# Article-boundary spatial workflow
bash run_core.sh spatial

# Refit all 16 RCTD sections before downstream spatial calculations
LIUF_SPATIAL_RUN_VARIANT=raw_rctd_diagnostic bash run_core.sh spatial

# Retrain the MOFA K grid and final/replay K=15 models
MOFA_RUN_RAW_DIAGNOSTIC=1 MOFA_FORCE_RETRAIN=1 MOFA_WORKERS=4 \
  bash run_core.sh serum_mofa
```

The default MOFA downstream branch uses the declared historical K=15 article
model. After raw retraining, set `MOFA_ANALYSIS_SOURCE=raw_derived` to calculate
downstream results from the new model. Diagnostic results require a separate
numerical equivalence assessment before being compared with article values.

## Random seeds

The master seed is `20260807`, exported as `LIUF_MASTER_SEED` and
`PYTHONHASHSEED`; BLAS/OpenMP thread counts default to one. Historical
step-specific seeds are retained where they affect accepted results.

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

Per-stage RNG settings and runtimes are recorded in the generated audit tables.

## Analysis boundaries

Exact article-equivalent analysis is not synonymous with an exclusive raw-only
recomputation. The following boundaries are explicit:

- Expert-curated Seurat objects are checksum-locked inputs; raw atlas
  reconstruction is a separate stage.
- The article spatial workflow validates 18 historical RCTD/CoVarNet objects
  before downstream calculations. Raw RCTD refitting is a diagnostic branch.
- Historical positive-flow CORNETO topology and current signal-state solutions
  are not treated as interchangeable. The current solver branch runs twice.
- Historical MOFA objects support article-boundary downstream calculations;
  raw retraining does not silently substitute a frozen model.
- Validation targets are used only for numerical comparison and never as
  silent substitutes for production computation.

See `docs/BOUNDARIES.md` for paths and manifests.

## Validation

```bash
python3 -m unittest discover -s code/tests -p 'test_*.py' -v

Rscript -e 'files <- list.files("code", pattern="[.]R$", recursive=TRUE, full.names=TRUE); for (f in files) parse(f); cat(length(files), "R scripts parsed\n")'

bash -n run_core.sh
for file in code/scripts/*.sh code/lib/shell/*.sh code/pipelines/*/*.sh; do bash -n "$file"; done
```

The tests check repository layout, configuration and key paths, method
coverage, production-code patterns, Python syntax, and every public hash in
`SOURCE_MANIFEST.tsv`.

## Provenance

`SOURCE_MANIFEST.tsv` preserves the original source path and original-byte
SHA-256 alongside the relocated public path, normalized public SHA-256, and
file size. `source_cleanup` marks scripts with numerical-output-preserving
cleanup; `layout_normalized` marks relocated sources; `documented_template`
marks repository-level entry points and documentation.
