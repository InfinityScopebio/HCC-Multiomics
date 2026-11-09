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
to the diagnostic result root.

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

## Validation-only targets

Some MOFA downstream scripts compare newly calculated results with historical
validation targets under:

```text
data/validation/frozen_targets/MOFA/
```

These files are needed only for full article-equivalence audits. They are not
production substitutes for model fitting. Numeric comparisons are retained in
this workflow repository; figure/pixel comparisons are not included.

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
- `repo/audit`: input hashes, seed declarations, runtimes and repeatability
  receipts.

These directories are ignored by Git so local measurements and generated
results cannot be published accidentally.
