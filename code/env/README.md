# Environment inventory

## Validated runtimes

- R 4.5.2 (`/usr/local/bin/Rscript` in the validated local environment).
- General Python 3.12.2. Python 3.11 or newer is required for `tomllib` in
  `code/lib/python/project_config.py`.
- Spatial Python 3.10 or newer with the package versions in
  `code/env/python_c2l_freeze.txt`. The main runner checks scanpy 1.11.5,
  COMMOT 0.0.3, and CORNETO 0.9.1a6 before spatial analysis.

`R_sessionInfo.txt` records the article runtime. `R_packages.tsv` is a combined
inventory from the two R library locations listed in `R_libpaths.txt`; it is not
a single-library lockfile. Do not select package versions from only one of the
two rows without considering the complete R library search path.

## Installation

Use a compatible operating system and R/Python versions. Existing validated
environments are preferred. To create new environments from the inventories:

```bash
python3.12 -m venv .venv
.venv/bin/pip install -r code/env/python_base_freeze.txt

python3.10 -m venv .venv-c2l
.venv-c2l/bin/pip install -r code/env/python_c2l_freeze.txt
```

For R, use R 4.5.2 and install the packages recorded in `R_packages.tsv` from
the appropriate CRAN/Bioconductor or upstream sources. `R_packages.tsv` and
`R_sessionInfo.txt` document versions; they do not constitute an automated R
environment installer.

## Runtime overrides

```bash
export LIUF_RSCRIPT=/path/to/Rscript
export LIUF_PYTHON=/path/to/python3.12
export LIUF_C2L_PYTHON=/path/to/spatial/python
export MOFA_PYTHON=/path/to/mofa/python
```

The STORM and serum MOFA analyses use MSigDB release `2025.1.Hs`. Their scripts
require the validated local `msigdbr` cache and its recorded checksum, then
verify the returned release and gene-set hash. They do not select or download a
different release at runtime. `STORM_MSIGDB_VERSION` and `MOFA_MSIGDB_VERSION`
declare the expected release; any mismatch is fatal.

An explicitly invalid `LIUF_C2L_PYTHON` is fatal. A fallback interpreter is
never selected silently after an explicit override.

## Numerical execution

The master seed is `20260807`; BLAS/OpenMP thread counts default to one.
Package versions and per-stage seeds are necessary but not sufficient for
reproducibility. Provenance checks, immutable-input hashes, and repeatability
audits are also required.
