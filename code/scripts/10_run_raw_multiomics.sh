#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/code/lib/shell/repro_env.sh"
liuf_init_repro_env "${PROJECT_ROOT}"

R_SCRIPT="${LIUF_RSCRIPT}"
STAGE_DIR="${PROJECT_ROOT}/results/00_preprocess"
INPUT_LINK="${STAGE_DIR}/LiuF_data"
PRIVATE_INPUT="${PROJECT_ROOT}/data/private/raw/data"
LOG_DIR="${PROJECT_ROOT}/logs/00_preprocess"

if [[ ! -d "${PRIVATE_INPUT}" ]]; then
  echo "Materialized private input is missing: ${PRIVATE_INPUT}" >&2
  exit 1
fi

mkdir -p "${STAGE_DIR}" "${LOG_DIR}"
if [[ -e "${INPUT_LINK}" && ! -L "${INPUT_LINK}" ]]; then
  echo "Refusing to replace non-symlink input path: ${INPUT_LINK}" >&2
  exit 1
fi
ln -sfn "${PRIVATE_INPUT}" "${INPUT_LINK}"

"${R_SCRIPT}" -e '
required <- c("data.table", "readxl", "org.Hs.eg.db", "AnnotationDbi")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "))
}
' >"${LOG_DIR}/package_check.log" 2>&1

(
  cd "${STAGE_DIR}"
  "${R_SCRIPT}" \
    "${PROJECT_ROOT}/code/vendor/00_multiomics/run_m0_sample_sheet.R" \
    >"${LOG_DIR}/10_sample_sheet.stdout.log" \
    2>"${LOG_DIR}/10_sample_sheet.stderr.log"
  "${R_SCRIPT}" \
    "${PROJECT_ROOT}/code/vendor/00_multiomics/run_m0_load_qc.R" \
    >"${LOG_DIR}/11_load_qc.stdout.log" \
    2>"${LOG_DIR}/11_load_qc.stderr.log"
)

printf 'Raw multi-omics preprocessing complete: %s\n' \
  "${STAGE_DIR}/multiomics"
