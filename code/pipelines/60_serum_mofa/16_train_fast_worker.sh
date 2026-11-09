#!/usr/bin/env bash

set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: 16_train_fast_worker.sh K seed" >&2
  exit 2
fi

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "${PROJECT_ROOT}/code/lib/shell/repro_env.sh"
liuf_init_repro_env "${PROJECT_ROOT}"

PIPELINE_DIR="${PROJECT_ROOT}/code/pipelines/60_serum_mofa"
RESULT_DIR="${PROJECT_ROOT}/results/60_serum_mofa"
LOG_DIR="${RESULT_DIR}/logs"
RSCRIPT="${LIUF_RSCRIPT}"
k="$1"
seed="$2"

liuf_require_positive_integer "${k}" "K"
liuf_require_r_seed "${seed}" "MOFA worker seed"

mkdir -p "${LOG_DIR}"
export MOFA_WORKER_SEED="${seed}"
export PYTHONHASHSEED="${seed}"
export SKLEARN_RANDOM_STATE="${seed}"

exec "${RSCRIPT}" "${PIPELINE_DIR}/15_train_mofa_model.R" \
  "--profile=fast1000" "--k=${k}" "--seed=${seed}" \
  >"${LOG_DIR}/train_fast1000_K${k}_seed${seed}.log" 2>&1
