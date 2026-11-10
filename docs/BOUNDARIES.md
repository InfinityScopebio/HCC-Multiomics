# Historical and validation boundaries

## Spatial article boundary

For article-equivalent downstream spatial analysis, these local inputs are
validated by `config/manifests/spatial_historical_intermediates.tsv`:

```text
data/curation/spatial/RCTD/ST3_rctd_reference.rds
data/curation/spatial/RCTD/weights/*_weights.rds  # exactly 16
data/curation/spatial/CoVarNet/cm_network.rds
```

The default spatial stage does not refit those 18 objects. The explicit
`LIUF_SPATIAL_RUN_VARIANT=raw_rctd_diagnostic` branch recomputes RCTD and writes
to the diagnostic result root. Article-equivalent provenance does not assert
byte-level identity of newly computed COMMOT caches or CCI tables with historical
validation targets; numerical comparisons are evaluated separately.

The article-used positive-flow CORNETO topology has a separate manifest,
`config/manifests/spatial_historical_corneto.tsv`:

```text
data/curation/spatial/CORNETO/step5_corneto_per_receptor.csv
data/curation/spatial/CORNETO/step5_corneto_measurements.csv
```

These objects document the historical solver selection; they are not presented
as a new CORNETO optimization. The current signal-state solver branch in
`80_corneto_mapk.py` recomputes and compares two independent runs.

## Serum MOFA article boundary

`config/manifests/serum_mofa_historical_boundary.tsv` registers:

```text
data/curation/serum_mofa/incdec_input.rds
data/curation/serum_mofa/incdec_Kselect_summary.csv
data/curation/serum_mofa/incdec_mofa_K15.hdf5
```

The default downstream branch uses this historical native K=15 model. Raw
retraining requires the M0 serum matrices and never falls back to this model
when a requested raw model is missing or has invalid provenance.

Raw-derived execution skips historical materialization, article replay, and the
article finalizer. `61_finalize_raw_audit.R` instead validates the current
model grid and final/replay metadata, binds downstream audits to the same input
and model hashes, and writes `mofa_raw_verdict.tsv`, `mofa_raw_check_status.tsv`,
`mofa_raw_output_checksums.tsv`, and `MOFA_RAW_AUDIT.md`. Article numerical
equivalence is established only when every designated numerical check is exact;
other completed results are explicitly diagnostic differences. Running with
`MOFA_RUN_RAW_DIAGNOSTIC=0` audits existing raw outputs and is not a fresh
rebuild.

## CM02 downstream manifests

`config/manifests/cm02_historical_downstream.tsv` and
`config/manifests/supplemental2_rank_survey_historical.tsv` document historical
CM02 inputs. `code/pipelines/20_scrna_cm02/10_materialize_historical_downstream.R`
and `20_materialize_supplemental2_rank.R` validate and materialize those declared
boundaries when historical replay is explicitly requested; they are not silently
used by the raw-derived runner.

## Validation-only targets

Some MOFA downstream scripts compare newly calculated results with historical
validation targets under:

```text
data/validation/frozen_targets/MOFA/
```

These files are needed only for full article-equivalence audits. They are not
production substitutes for model fitting. Numerical comparisons and execution
records are retained as tables.

## Generated output locations

- `results/00_preprocess`: M0 matrices and sample sheet.
- `results/20_scrna_cm02`: CoVarNet, MuSiC, trajectory, survival, Hep_C
  activities, composition and RoE tables.
- `results/30_tissue_multiomics`: pathway, cross-layer, KSEA and paired
  feature tables.
- `results/40_tissue_resistance`: DepMap/PRISM scores and composite
  associations.
- `results/panels/main/4E`: STORM limma ranks and GSEA result tables.
- `results/50_spatial_cci`: RCTD/CM recovery, radial/pathway tables,
  COMMOT caches, patient-balanced CCI and CORNETO routes.
- `results/60_serum_mofa`: MOFA models, TCGA axes, survival, classifier
  and supplemental tables.
- `audit/`: input hashes, seed declarations, runtimes and repeatability
  receipts.

These directories are ignored by Git so local measurements and generated
results cannot be published accidentally.
