# Data root

Place local inputs at the paths documented in `docs/INPUTS.md` and
`docs/BOUNDARIES.md`. Small, explicitly inventoried annotation and curation
inputs are tracked under `data/reference` and `data/curation/resistance`.

Private matrices, public-cohort files, validation objects, curated Seurat or
spatial objects, and serum MOFA objects are supplied locally and excluded from
version control. Numerical outputs are written to the sibling `results/`
directory; execution records, logs, and generated inventories are written to
`audit/`, `logs/`, and `manifests/`.
