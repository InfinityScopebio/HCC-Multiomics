#!/usr/bin/env bash

# Bash 3.2-compatible bounded process pool with batch-level fail-fast behavior.

LIUF_WORKER_POOL_SHELL_DIR="$(
  cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P
)"
LIUF_WORKER_POOL_GROUP_EXEC="$(
  cd "${LIUF_WORKER_POOL_SHELL_DIR}/../python" &&
    pwd -P
)/process_group_exec.py"

liuf_worker_pool_clear_state() {
  LIUF_WORKER_POOL_PIDS=()
  LIUF_WORKER_POOL_LABELS=()
  LIUF_WORKER_POOL_READY_FILES=()
  LIUF_WORKER_POOL_ACK_FILES=()
  LIUF_WORKER_POOL_GROUP_FILES=()
  LIUF_WORKER_POOL_GROUP_IDS=()
  LIUF_WORKER_POOL_GROUP_READY=()
}

liuf_worker_pool_group_running() {
  local group_id="$1"

  kill -0 "-${group_id}" 2>/dev/null
}

liuf_worker_pool_process_running() {
  local process_id="$1"

  liuf_worker_pool_group_running "${process_id}" ||
    kill -0 "${process_id}" 2>/dev/null
}

liuf_worker_pool_terminate_group() {
  local group_id="$1"
  local group_ready="${2:-1}"
  local tick=0

  if [[ "${group_ready}" == "1" ]]; then
    kill -TERM "-${group_id}" 2>/dev/null || return 0
  else
    kill -TERM "${group_id}" 2>/dev/null || true
  fi
  while [[ "${tick}" -lt "${LIUF_WORKER_POOL_TERM_GRACE_TICKS}" ]]; do
    if [[ "${group_ready}" == "1" ]]; then
      liuf_worker_pool_group_running "${group_id}" || return 0
    elif ! kill -0 "${group_id}" 2>/dev/null; then
      return
    fi
    sleep 0.1
    tick=$((tick + 1))
  done
  if [[ "${group_ready}" == "1" ]]; then
    kill -KILL "-${group_id}" 2>/dev/null || true
  else
    kill -KILL "${group_id}" 2>/dev/null || true
  fi
}

liuf_worker_pool_remove_traps() {
  trap - EXIT INT TERM HUP
}

liuf_worker_pool_terminate_all() {
  local pid
  local index
  local group_ready
  local group_id
  local ready_pid
  local ready_file
  local acknowledgement_file
  local group_file

  if [[ "${#LIUF_WORKER_POOL_PIDS[@]}" -gt 0 ]]; then
    for ((index = 0; index < ${#LIUF_WORKER_POOL_PIDS[@]}; index++)); do
      pid="${LIUF_WORKER_POOL_PIDS[${index}]}"
      [[ -z "${pid}" ]] && continue
      group_ready="${LIUF_WORKER_POOL_GROUP_READY[${index}]:-0}"
      group_id="${LIUF_WORKER_POOL_GROUP_IDS[${index}]:-}"
      group_file="${LIUF_WORKER_POOL_GROUP_FILES[${index}]:-}"
      if [[ "${group_ready}" != "1" && -n "${group_file}" && -s "${group_file}" ]]; then
        group_id="$(tr -d '\r\n' < "${group_file}")"
        if [[ "${group_id}" =~ ^[1-9][0-9]*$ ]]; then
          group_ready="1"
          LIUF_WORKER_POOL_GROUP_IDS[${index}]="${group_id}"
          LIUF_WORKER_POOL_GROUP_READY[${index}]="1"
        fi
      fi
      if [[ "${group_ready}" == "1" ]]; then
        liuf_worker_pool_terminate_group "${group_id}" 1
      fi
      liuf_worker_pool_terminate_group "${pid}" 0
      wait "${pid}" 2>/dev/null || true
      LIUF_WORKER_POOL_PIDS[${index}]=""
      LIUF_WORKER_POOL_GROUP_IDS[${index}]=""
    done
  fi
  for ready_file in "${LIUF_WORKER_POOL_READY_FILES[@]}"; do
    rm -f "${ready_file}"
  done
  for acknowledgement_file in "${LIUF_WORKER_POOL_ACK_FILES[@]}"; do
    rm -f "${acknowledgement_file}"
  done
  for group_file in "${LIUF_WORKER_POOL_GROUP_FILES[@]}"; do
    rm -f "${group_file}"
  done
  liuf_worker_pool_clear_state
  if [[ -n "${LIUF_WORKER_POOL_READY_DIR:-}" ]]; then
    rmdir "${LIUF_WORKER_POOL_READY_DIR}" 2>/dev/null || true
    LIUF_WORKER_POOL_READY_DIR=""
  fi
}

liuf_worker_pool_handle_exit() {
  local exit_code=$?

  if [[ "${LIUF_WORKER_POOL_ACTIVE:-0}" == "1" ]]; then
    liuf_worker_pool_terminate_all
  fi
  return "${exit_code}"
}

liuf_worker_pool_handle_signal() {
  local signal_name="$1"
  local exit_code="$2"

  liuf_worker_pool_terminate_all
  LIUF_WORKER_POOL_ACTIVE=0
  liuf_worker_pool_remove_traps
  printf 'worker pool interrupted by %s\n' "${signal_name}" >&2
  exit "${exit_code}"
}

liuf_worker_pool_init() {
  local limit="$1"

  liuf_require_positive_integer "${limit}" "worker pool limit" || return
  if [[ "${LIUF_WORKER_POOL_ACTIVE:-0}" == "1" ]]; then
    liuf_fail "worker pool is already active"
    return
  fi
  LIUF_WORKER_POOL_LIMIT="${limit}"
  LIUF_WORKER_POOL_TERM_GRACE_TICKS="${LIUF_WORKER_POOL_TERM_GRACE_TICKS:-50}"
  LIUF_WORKER_POOL_READY_TIMEOUT_TICKS="${LIUF_WORKER_POOL_READY_TIMEOUT_TICKS:-500}"
  liuf_require_positive_integer \
    "${LIUF_WORKER_POOL_TERM_GRACE_TICKS}" \
    "LIUF_WORKER_POOL_TERM_GRACE_TICKS" || return
  liuf_require_positive_integer \
    "${LIUF_WORKER_POOL_READY_TIMEOUT_TICKS}" \
    "LIUF_WORKER_POOL_READY_TIMEOUT_TICKS" || return
  if [[ ! -f "${LIUF_WORKER_POOL_GROUP_EXEC}" ]]; then
    liuf_fail "worker process-group launcher is missing"
    return
  fi
  liuf_require_command "${LIUF_PYTHON:-python3}" "LIUF_PYTHON" || return
  LIUF_WORKER_POOL_READY_DIR="$(
    mktemp -d "${TMPDIR:-/tmp}/liuf-worker-pool.XXXXXX"
  )"
  liuf_worker_pool_clear_state
  LIUF_WORKER_POOL_ACTIVE=1
  trap 'liuf_worker_pool_handle_exit' EXIT
  trap 'liuf_worker_pool_handle_signal INT 130' INT
  trap 'liuf_worker_pool_handle_signal TERM 143' TERM
  trap 'liuf_worker_pool_handle_signal HUP 129' HUP
}

liuf_worker_pool_wait_batch() {
  local failed=0
  local index
  local process_id
  local return_code
  local ready_file
  local acknowledgement_file
  local group_file
  local group_id

  for ((index = 0; index < ${#LIUF_WORKER_POOL_PIDS[@]}; index++)); do
    process_id="${LIUF_WORKER_POOL_PIDS[${index}]}"
    [[ -z "${process_id}" ]] && continue
    if wait "${process_id}"; then
      return_code=0
    else
      return_code=$?
    fi
    group_id="${LIUF_WORKER_POOL_GROUP_IDS[${index}]:-}"
    group_file="${LIUF_WORKER_POOL_GROUP_FILES[${index}]:-}"
    if [[ -z "${group_id}" && -n "${group_file}" && -s "${group_file}" ]]; then
      group_id="$(tr -d '\r\n' < "${group_file}")"
    fi
    if (
      [[ "${group_id}" =~ ^[1-9][0-9]*$ ]] &&
      liuf_worker_pool_group_running "${group_id}"
    ); then
      liuf_worker_pool_terminate_group "${group_id}" 1
      if [[ "${return_code}" == "0" ]]; then
        return_code=125
      fi
    fi
    LIUF_WORKER_POOL_PIDS[${index}]=""
    LIUF_WORKER_POOL_GROUP_IDS[${index}]=""
    LIUF_WORKER_POOL_GROUP_READY[${index}]="0"
    if [[ "${return_code}" != "0" ]]; then
      printf 'worker failed (%s): exit %s\n' \
        "${LIUF_WORKER_POOL_LABELS[${index}]}" \
        "${return_code}" >&2
      failed=1
    fi
  done
  for ready_file in "${LIUF_WORKER_POOL_READY_FILES[@]}"; do
    rm -f "${ready_file}"
  done
  for acknowledgement_file in "${LIUF_WORKER_POOL_ACK_FILES[@]}"; do
    rm -f "${acknowledgement_file}"
  done
  for group_file in "${LIUF_WORKER_POOL_GROUP_FILES[@]}"; do
    rm -f "${group_file}"
  done
  liuf_worker_pool_clear_state
  return "${failed}"
}

liuf_worker_pool_submit() {
  local label="$1"
  local process_id
  local ready_file
  local acknowledgement_file
  local acknowledgement_tmp
  local group_file
  local group_id
  local ready_pid
  local tick=0
  shift

  ready_file="${LIUF_WORKER_POOL_READY_DIR}/ready-${#LIUF_WORKER_POOL_PIDS[@]}"
  acknowledgement_file="${LIUF_WORKER_POOL_READY_DIR}/ack-${#LIUF_WORKER_POOL_PIDS[@]}"
  group_file="${LIUF_WORKER_POOL_READY_DIR}/group-${#LIUF_WORKER_POOL_PIDS[@]}"
  "${LIUF_PYTHON:-python3}" "${LIUF_WORKER_POOL_GROUP_EXEC}" \
    --ready-file "${ready_file}" \
    --ack-file "${acknowledgement_file}" \
    --group-file "${group_file}" -- "$@" &
  process_id=$!
  LIUF_WORKER_POOL_PIDS[${#LIUF_WORKER_POOL_PIDS[@]}]="${process_id}"
  LIUF_WORKER_POOL_LABELS[${#LIUF_WORKER_POOL_LABELS[@]}]="${label}"
  LIUF_WORKER_POOL_READY_FILES[
    ${#LIUF_WORKER_POOL_READY_FILES[@]}
  ]="${ready_file}"
  LIUF_WORKER_POOL_ACK_FILES[
    ${#LIUF_WORKER_POOL_ACK_FILES[@]}
  ]="${acknowledgement_file}"
  LIUF_WORKER_POOL_GROUP_FILES[
    ${#LIUF_WORKER_POOL_GROUP_FILES[@]}
  ]="${group_file}"
  LIUF_WORKER_POOL_GROUP_IDS[
    ${#LIUF_WORKER_POOL_GROUP_IDS[@]}
  ]=""
  LIUF_WORKER_POOL_GROUP_READY[
    ${#LIUF_WORKER_POOL_GROUP_READY[@]}
  ]="0"
  while [[ ! -s "${ready_file}" ]]; do
    if ! kill -0 "${process_id}" 2>/dev/null; then
      printf 'worker failed before process-group readiness (%s)\n' \
        "${label}" >&2
      return 125
    fi
    if [[ "${tick}" -ge "${LIUF_WORKER_POOL_READY_TIMEOUT_TICKS}" ]]; then
      printf 'worker process-group readiness timed out (%s)\n' \
        "${label}" >&2
      liuf_worker_pool_terminate_group "${process_id}" 0
      return 124
    fi
    sleep 0.01
    tick=$((tick + 1))
  done
  ready_pid="$(tr -d '\r\n' < "${ready_file}")"
  if [[ "${ready_pid}" != "${process_id}" ]]; then
    printf 'worker process-group PID mismatch (%s): %s != %s\n' \
      "${label}" "${ready_pid}" "${process_id}" >&2
    liuf_worker_pool_terminate_group "${process_id}" 1
    return 125
  fi
  acknowledgement_tmp="${acknowledgement_file}.tmp-$$"
  printf '%s\n' "${process_id}" > "${acknowledgement_tmp}"
  mv -f "${acknowledgement_tmp}" "${acknowledgement_file}"
  tick=0
  while [[ ! -s "${group_file}" ]]; do
    if ! kill -0 "${process_id}" 2>/dev/null; then
      printf 'worker supervisor failed before command readiness (%s)\n' \
        "${label}" >&2
      return 125
    fi
    if [[ "${tick}" -ge "${LIUF_WORKER_POOL_READY_TIMEOUT_TICKS}" ]]; then
      printf 'worker command process-group readiness timed out (%s)\n' \
        "${label}" >&2
      liuf_worker_pool_terminate_group "${process_id}" 0
      return 124
    fi
    sleep 0.01
    tick=$((tick + 1))
  done
  group_id="$(tr -d '\r\n' < "${group_file}")"
  if ! [[ "${group_id}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'worker command process-group PID is invalid (%s): %s\n' \
      "${label}" "${group_id}" >&2
    liuf_worker_pool_terminate_group "${process_id}" 0
    return 125
  fi
  LIUF_WORKER_POOL_GROUP_IDS[
    $((${#LIUF_WORKER_POOL_GROUP_IDS[@]} - 1))
  ]="${group_id}"
  LIUF_WORKER_POOL_GROUP_READY[
    $((${#LIUF_WORKER_POOL_GROUP_READY[@]} - 1))
  ]="1"
  if [[ "${#LIUF_WORKER_POOL_PIDS[@]}" -ge "${LIUF_WORKER_POOL_LIMIT}" ]]; then
    liuf_worker_pool_wait_batch
  fi
}

liuf_worker_pool_finish() {
  if [[ "${#LIUF_WORKER_POOL_PIDS[@]}" -gt 0 ]]; then
    liuf_worker_pool_wait_batch || return
  fi
  LIUF_WORKER_POOL_ACTIVE=0
  liuf_worker_pool_remove_traps
  if [[ -n "${LIUF_WORKER_POOL_READY_DIR:-}" ]]; then
    rmdir "${LIUF_WORKER_POOL_READY_DIR}" 2>/dev/null || true
    LIUF_WORKER_POOL_READY_DIR=""
  fi
}
