# Local data and output root

This directory is reserved for local inputs and generated outputs. It is not a
source of published raw data or frozen result objects in this repository.

Expected local input families are documented in `docs/INPUTS.md`. Generated
analysis tables and models are written to `repo/results`; execution and seed
audits are written to `repo/audit`. Those paths are ignored by Git, except for
the small, explicitly inventoried annotation and curation inputs already
tracked under `repo/data/reference` and `repo/data/curation/resistance`.
