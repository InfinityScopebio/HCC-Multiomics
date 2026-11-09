#!/usr/bin/env bash
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${PIPELINE_DIR}/../../.." && pwd)"
source "${PROJECT_ROOT}/code/lib/shell/repro_env.sh"
liuf_init_repro_env "${PROJECT_ROOT}"

ROOT="${LIUFANG_ROOT}"
PIPELINE_DIR="${ROOT}/code/pipelines/30_tissue_multiomics"
RESULT_DIR="${ROOT}/results/30_tissue_multiomics"
LOG_DIR="${RESULT_DIR}/logs"
R="${LIUF_RSCRIPT}"
PYTHON="${LIUF_PYTHON}"

mkdir -p "${LOG_DIR}"


run_r() {
  local script="$1"
  local name
  name="$(basename "${script}" .R)"
  "${R}" "${PIPELINE_DIR}/${script}" 2>&1 | tee "${LOG_DIR}/${name}.log"
}

run_r 00_validate_inputs.R
run_r 05_freeze_local_msigdb.R
run_r 10_progeny.R
run_r 20_ptmsea.R
run_r 22_run_ptmsea_rng_probe.R
run_r 23_summarize_ptmsea_rng_modes.R
run_r 30_cross_layer.R
run_r 40_ksea_three_layer.R
run_r 60_fig4c.R
run_r 92_seed_audit.R

printf 'Tissue multi-omics pipeline completed: %s\n' "${RESULT_DIR}"
