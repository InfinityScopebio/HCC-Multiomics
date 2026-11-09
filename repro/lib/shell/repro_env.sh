#!/usr/bin/env bash

# Shared deterministic runtime contract for every reproducibility pipeline.

liuf_fail() {
  printf 'repro environment error: %s\n' "$1" >&2
  return 2
}

liuf_require_r_seed() {
  local value="$1"
  local name="$2"
  local canonical_value

  if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
    liuf_fail "${name} must be an integer in [0, 2147483647]."
    return
  fi
  if [[ "${value}" != "0" && "${value}" == 0* ]]; then
    liuf_fail "${name} must use canonical decimal notation without leading zeros."
    return
  fi
  canonical_value="${value}"
  if (
    ((${#canonical_value} > 10)) ||
      (
        ((${#canonical_value} == 10)) &&
          [[ "${canonical_value}" > "2147483647" ]]
      )
  ); then
    liuf_fail "${name} must be an integer in [0, 2147483647]."
    return
  fi
}

liuf_require_positive_integer() {
  local value="$1"
  local name="$2"

  if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
    liuf_fail "${name} must be a positive integer."
    return
  fi
}

liuf_require_command() {
  local command_name="$1"
  local variable_name="$2"

  if ! command -v "${command_name}" >/dev/null 2>&1; then
    liuf_fail "${variable_name} runtime is unavailable: ${command_name}"
    return
  fi
}

liuf_require_r_runtime() {
  local rscript="$1"
  local python="$2"
  local runtime_marker

  if ! runtime_marker="$(
    "${python}" - "${rscript}" <<'PY'
import subprocess
import sys

rscript = sys.argv[1]
expected = (
    "ba7816bf8f01cfea414140de5dae2223"
    "b00361a396177a9cb410ff61f20015ad"
)
probe = (
    'minimum <- package_version("4.5.0"); '
    'has_sha256 <- exists("sha256sum", '
    'envir = asNamespace("tools"), inherits = FALSE); '
    'probe_path <- tempfile("liuf-r-sha256-"); '
    'writeBin(charToRaw("abc"), probe_path); '
    "probe_sha256 <- if (has_sha256) "
    "unname(tools::sha256sum(probe_path)) else ''; "
    "unlink(probe_path); "
    "if (getRversion() < minimum || !has_sha256 || "
    f'probe_sha256 != "{expected}") quit(status = 2L); '
    'cat("LIUF_R_RUNTIME_OK")'
)
try:
    result = subprocess.run(
        [rscript, "--vanilla", "-e", probe],
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
except subprocess.TimeoutExpired:
    raise SystemExit(124)
if result.returncode != 0:
    raise SystemExit(result.returncode)
sys.stdout.write(result.stdout)
PY
  )"; then
    liuf_fail \
      "LIUF_RSCRIPT must provide R >= 4.5.0 and tools::sha256sum: ${rscript}"
    return
  fi
  if [[ "${runtime_marker}" != "LIUF_R_RUNTIME_OK" ]]; then
    liuf_fail \
      "LIUF_RSCRIPT failed the R >= 4.5.0 runtime probe: ${rscript}"
    return
  fi
}

liuf_init_repro_env() {
  local derived_root="$1"
  local canonical_root
  local configured_root

  if [[ ! -d "${derived_root}" ]]; then
    liuf_fail "project root does not exist: ${derived_root}"
    return
  fi
  canonical_root="$(cd "${derived_root}" && pwd -P)"
  if [[ ! -f "${canonical_root}/repro/config/project.toml" ]]; then
    liuf_fail "project root lacks repro/config/project.toml: ${canonical_root}"
    return
  fi

  if [[ -n "${LIUFANG_ROOT:-}" ]]; then
    if [[ ! -d "${LIUFANG_ROOT}" ]]; then
      liuf_fail "LIUFANG_ROOT does not exist: ${LIUFANG_ROOT}"
      return
    fi
    configured_root="$(cd "${LIUFANG_ROOT}" && pwd -P)"
    if [[ "${configured_root}" != "${canonical_root}" ]]; then
      liuf_fail \
        "LIUFANG_ROOT (${configured_root}) does not match ${canonical_root}."
      return
    fi
  fi

  LIUFANG_ROOT="${canonical_root}"
  LIUF_MASTER_SEED="${LIUF_MASTER_SEED:-20260807}"
  LIUF_THREADS="${LIUF_THREADS:-1}"
  if [[ -z "${LIUF_RSCRIPT:-}" ]]; then
    if [[ -x "/usr/local/bin/Rscript" ]]; then
      LIUF_RSCRIPT="/usr/local/bin/Rscript"
    else
      LIUF_RSCRIPT="Rscript"
    fi
  fi
  if [[ -z "${LIUF_PYTHON:-}" ]]; then
    if [[ -x "${HOME}/miniconda3/bin/python3" ]]; then
      LIUF_PYTHON="${HOME}/miniconda3/bin/python3"
    else
      LIUF_PYTHON="python3"
    fi
  fi

  liuf_require_r_seed "${LIUF_MASTER_SEED}" "LIUF_MASTER_SEED" || return
  liuf_require_positive_integer "${LIUF_THREADS}" "LIUF_THREADS" || return
  liuf_require_command "${LIUF_RSCRIPT}" "LIUF_RSCRIPT" || return
  liuf_require_command "${LIUF_PYTHON}" "LIUF_PYTHON" || return
  liuf_require_r_runtime "${LIUF_RSCRIPT}" "${LIUF_PYTHON}" || return

  export LIUFANG_ROOT
  export LIUF_MASTER_SEED
  export LIUF_THREADS
  export LIUF_RSCRIPT
  export LIUF_PYTHON
  export PYTHONHASHSEED="${LIUF_MASTER_SEED}"
  export OMP_NUM_THREADS="${LIUF_THREADS}"
  export OPENBLAS_NUM_THREADS="${LIUF_THREADS}"
  export MKL_NUM_THREADS="${LIUF_THREADS}"
  export VECLIB_MAXIMUM_THREADS="${LIUF_THREADS}"
  export NUMEXPR_NUM_THREADS="${LIUF_THREADS}"
  export BLIS_NUM_THREADS="${LIUF_THREADS}"
  export RCPP_PARALLEL_NUM_THREADS="${LIUF_THREADS}"
  export OMP_DYNAMIC="FALSE"
  export MKL_DYNAMIC="FALSE"
  export R_FUTURE_RNG_ONMISUSE="error"
}
