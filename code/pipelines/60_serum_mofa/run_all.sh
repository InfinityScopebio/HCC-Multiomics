#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "${PROJECT_ROOT}/code/lib/shell/repro_env.sh"
liuf_init_repro_env "${PROJECT_ROOT}"
source "${PROJECT_ROOT}/code/lib/shell/bounded_worker_pool.sh"

PIPELINE_DIR="${PROJECT_ROOT}/code/pipelines/60_serum_mofa"
CANONICAL_RESULT_DIR="${PROJECT_ROOT}/results/60_serum_mofa"
RESULT_DIR="${CANONICAL_RESULT_DIR}"
if [[ "${MOFA_ANALYSIS_SOURCE:-article}" == "raw_derived" ]]; then
  RESULT_DIR="${CANONICAL_RESULT_DIR}/raw_to_figure"
fi
LOG_DIR="${RESULT_DIR}/logs"
RSCRIPT="${LIUF_RSCRIPT}"
MOFA_WORKERS="${MOFA_WORKERS:-4}"
MOFA_FORCE_RETRAIN="${MOFA_FORCE_RETRAIN:-0}"
MOFA_RUN_RAW_DIAGNOSTIC="${MOFA_RUN_RAW_DIAGNOSTIC:-0}"
MOFA_ANALYSIS_SOURCE="${MOFA_ANALYSIS_SOURCE:-article}"
if [[ "${MOFA_ANALYSIS_SOURCE}" != "article" && "${MOFA_ANALYSIS_SOURCE}" != "raw_derived" ]]; then
  liuf_fail "MOFA_ANALYSIS_SOURCE must be article or raw_derived."
fi
for variable in MOFA_RUN_RAW_DIAGNOSTIC MOFA_FORCE_RETRAIN; do
  if [[ "${!variable}" != "0" && "${!variable}" != "1" ]]; then
    liuf_fail "${variable} must be 0 or 1."
  fi
done
if [[ "${MOFA_ANALYSIS_SOURCE}" == "article" && "${MOFA_RUN_RAW_DIAGNOSTIC}" == "1" ]]; then
  liuf_fail "Raw MOFA retraining requires MOFA_ANALYSIS_SOURCE=raw_derived."
fi
export MOFA_ANALYSIS_SOURCE
MOFA_PERM_B="${MOFA_PERM_B:-1000}"
RUNTIME_FILE="${RESULT_DIR}/pipeline_stage_runtime.tsv"

export MOFA_FORCE_RETRAIN
export MOFA_RUN_RAW_DIAGNOSTIC
if [[ "${MOFA_PERM_B}" != "1000" ]]; then
  liuf_fail "MOFA_PERM_B is fixed at 1000 for the accepted analysis."
fi
export MOFA_PERM_B

mkdir -p "${LOG_DIR}"
liuf_require_positive_integer "${MOFA_WORKERS}" "MOFA_WORKERS"

printf 'stage\tstarted_epoch\tfinished_epoch\telapsed_seconds\n' >"${RUNTIME_FILE}"

run_stage() {
  local script="$1"
  local started
  local finished
  started="$(date +%s)"
  echo "[serum-mofa] ${script}"
  "${RSCRIPT}" "${PIPELINE_DIR}/${script}" \
    2>&1 | tee "${LOG_DIR}/${script%.R}.log"
  finished="$(date +%s)"
  printf '%s\t%s\t%s\t%s\n' \
    "${script}" "${started}" "${finished}" "$((finished - started))" \
    >>"${RUNTIME_FILE}"
}


train_model() {
  local profile="$1"
  local k="$2"
  local seed="$3"
  local log_file="${LOG_DIR}/train_${profile}_K${k}_seed${seed}.log"
  MOFA_WORKER_SEED="${seed}" \
    PYTHONHASHSEED="${seed}" \
    SKLEARN_RANDOM_STATE="${seed}" \
    "${RSCRIPT}" "${PIPELINE_DIR}/15_train_mofa_model.R" \
    "--profile=${profile}" "--k=${k}" "--seed=${seed}" \
    >"${log_file}" 2>&1
}

export PROJECT_ROOT PIPELINE_DIR RESULT_DIR LOG_DIR

if [[ "${MOFA_ANALYSIS_SOURCE}" != "raw_derived" ]]; then
  run_stage "05_materialize_historical_boundary.R"
fi
run_stage "10_rebuild_incdec_input.R"
if [[ "${MOFA_ANALYSIS_SOURCE}" != "raw_derived" ]]; then
  run_stage "18_article_k_scan_and_variance.R"
fi

if [[ "${MOFA_RUN_RAW_DIAGNOSTIC}" == "1" ]]; then
  training_started="$(date +%s)"
  echo "[serum-mofa] diagnostic final slow/2000 K15 seed 1"
  train_model "slow2000" 15 1

  echo "[serum-mofa] diagnostic fast/1000 K=2..20 x seeds=1..6 (${MOFA_WORKERS} workers)"
  job_file="${RESULT_DIR}/k_scan_training_jobs.tsv"
  : >"${job_file}"
  for k in $(seq 2 20); do
    for seed in $(seq 1 6); do
      printf '%s\t%s\n' "${k}" "${seed}" >>"${job_file}"
    done
  done
  liuf_worker_pool_init "${MOFA_WORKERS}"
  while IFS=$'\t' read -r k seed; do
    liuf_worker_pool_submit \
      "K${k}_seed${seed}" \
      env LIUF_RSCRIPT="${LIUF_RSCRIPT}" \
      "${PIPELINE_DIR}/16_train_fast_worker.sh" "${k}" "${seed}"
  done <"${job_file}"
  liuf_worker_pool_finish

  echo "[serum-mofa] diagnostic deterministic replay models"
  train_model "fast1000_replay" 15 1
  train_model "slow2000_replay" 15 1
  training_finished="$(date +%s)"
  printf '%s\t%s\t%s\t%s\n' \
    "15_train_mofa_models_diagnostic" \
    "${training_started}" \
    "${training_finished}" \
    "$((training_finished - training_started))" \
    >>"${RUNTIME_FILE}"
  run_stage "20_k_scan_and_variance.R"
else
  echo "[serum-mofa] raw MOFA retraining diagnostic skipped (set MOFA_RUN_RAW_DIAGNOSTIC=1 to run)"
fi
run_stage "25_rebuild_tcga_axes.R"
run_stage "30_tcga_survival.R"

classifier_started="$(date +%s)"
echo "[serum-mofa] 40_factor2_classifier.R (B=1000)"
"${RSCRIPT}" \
  "${PIPELINE_DIR}/40_factor2_classifier.R" \
  2>&1 | tee "${LOG_DIR}/40_factor2_classifier.log"
classifier_finished="$(date +%s)"
printf '%s\t%s\t%s\t%s\n' \
  "40_factor2_classifier.R" \
  "${classifier_started}" \
  "${classifier_finished}" \
  "$((classifier_finished - classifier_started))" \
  >>"${RUNTIME_FILE}"

run_stage "50_supplemental_tables.R"
if [[ "${MOFA_ANALYSIS_SOURCE}" == "raw_derived" ]]; then
  run_stage "61_finalize_raw_audit.R"
else
  run_stage "61_finalize_core_audit.R"
fi

echo "[serum-mofa] complete"
