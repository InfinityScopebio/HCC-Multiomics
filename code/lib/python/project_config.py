from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import shutil
import subprocess
import tomllib
from typing import Any


R_SEED_MAX = 2**31 - 1
R_SEED_MODULUS = 2**31
R_RUNTIME_MARKER = "LIUF_R_RUNTIME_OK"
MINIMUM_R_VERSION = "4.5.0"
SHA256_ABC = (
    "ba7816bf8f01cfea414140de5dae2223"
    "b00361a396177a9cb410ff61f20015ad"
)


def derive_r_seed(base_seed: int, offset: int = 0) -> int:
    """Derive an R-compatible seed without integer overflow."""

    if isinstance(base_seed, bool) or not isinstance(base_seed, int):
        raise TypeError("base_seed must be an integer")
    if isinstance(offset, bool) or not isinstance(offset, int):
        raise TypeError("offset must be an integer")
    if not 0 <= base_seed <= R_SEED_MAX:
        raise ValueError(f"base_seed must be in [0, {R_SEED_MAX}]")
    return (base_seed + offset) % R_SEED_MODULUS


@dataclass(frozen=True)
class ProjectPaths:
    root: Path
    repo: Path
    repro: Path
    data: Path
    results: Path
    logs: Path
    audit: Path
    manifests: Path
    private_raw_source: Path

    def require_runtime_output(self, path: Path) -> Path:
        resolved = path.resolve()
        allowed = (self.results, self.logs, self.audit, self.manifests)
        if not any(resolved == base or base in resolved.parents for base in allowed):
            raise ValueError(f"Runtime output is outside results/, audit/, logs/ or manifests/: {resolved}")
        return resolved


@dataclass(frozen=True)
class ValidationConfig:
    absolute_tolerance: float
    relative_tolerance: float
    require_sample_order: bool
    require_feature_order: bool


@dataclass(frozen=True)
class ExecutionConfig:
    default_seed: int
    threads: int
    fail_fast: bool


@dataclass(frozen=True)
class ProjectConfig:
    paths: ProjectPaths
    validation: ValidationConfig
    execution: ExecutionConfig
    raw: dict[str, Any]


def find_project_root(start: Path | None = None) -> Path:
    current = (start or Path(__file__)).resolve()
    if current.is_file():
        current = current.parent
    for candidate in (current, *current.parents):
        if (candidate / "config/project.toml").is_file() and (candidate / ".git").exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate project root containing config/project.toml and .git"
    )


def load_project_config(config_path: Path | None = None) -> ProjectConfig:
    root = find_project_root()
    path = config_path or root / "config" / "project.toml"
    with path.open("rb") as handle:
        raw = tomllib.load(handle)

    project = raw["project"]
    paths_raw = raw["paths"]
    repo = (root / project["repo_dir"]).resolve()
    repro = (root / project["repro_dir"]).resolve()
    paths = ProjectPaths(
        root=root,
        repo=repo,
        repro=repro,
        data=(root / paths_raw["data_dir"]).resolve(),
        results=(root / paths_raw["result_dir"]).resolve(),
        logs=(root / paths_raw["log_dir"]).resolve(),
        audit=(root / paths_raw["audit_dir"]).resolve(),
        manifests=(root / paths_raw["manifest_dir"]).resolve(),
        private_raw_source=(root / paths_raw["private_raw_source"]).resolve(),
    )
    validation_raw = raw["validation"]
    validation = ValidationConfig(
        absolute_tolerance=float(validation_raw["numeric_absolute_tolerance"]),
        relative_tolerance=float(validation_raw["numeric_relative_tolerance"]),
        require_sample_order=bool(validation_raw["require_sample_order"]),
        require_feature_order=bool(validation_raw["require_feature_order"]),
    )
    execution_raw = raw["execution"]
    default_seed = execution_raw["default_seed"]
    threads = execution_raw["threads"]
    if type(default_seed) is not int:
        raise TypeError("execution.default_seed must be a TOML integer")
    if type(threads) is not int:
        raise TypeError("execution.threads must be a TOML integer")
    if not 0 <= default_seed <= R_SEED_MAX:
        raise ValueError(
            "execution.default_seed must be an integer in "
            f"[0, {R_SEED_MAX}]"
        )
    if threads < 1:
        raise ValueError("execution.threads must be a positive integer")
    execution = ExecutionConfig(
        default_seed=default_seed,
        threads=threads,
        fail_fast=bool(execution_raw["fail_fast"]),
    )
    return ProjectConfig(
        paths=paths,
        validation=validation,
        execution=execution,
        raw=raw,
    )


def validate_rscript_runtime(
    path: str,
    timeout_seconds: float = 10.0,
) -> str:
    probe = (
        f'minimum <- package_version("{MINIMUM_R_VERSION}"); '
        'has_sha256 <- exists("sha256sum", '
        'envir = asNamespace("tools"), inherits = FALSE); '
        'probe_path <- tempfile("liuf-r-sha256-"); '
        'writeBin(charToRaw("abc"), probe_path); '
        "probe_sha256 <- if (has_sha256) "
        "unname(tools::sha256sum(probe_path)) else ''; "
        "unlink(probe_path); "
        "if (getRversion() < minimum || !has_sha256 || "
        f'probe_sha256 != "{SHA256_ABC}") quit(status = 2L); '
        f'cat("{R_RUNTIME_MARKER}")'
    )
    try:
        result = subprocess.run(
            [path, "--vanilla", "-e", probe],
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"Rscript runtime probe timed out after {timeout_seconds:g}s: {path}"
        ) from error
    if result.returncode != 0 or result.stdout.strip() != R_RUNTIME_MARKER:
        detail = result.stderr.strip() or result.stdout.strip() or "probe failed"
        raise RuntimeError(
            "Rscript must provide "
            f"R >= {MINIMUM_R_VERSION} and tools::sha256sum: {path} ({detail})"
        )
    return path


def resolve_rscript() -> str:
    requested = os.environ.get("LIUF_RSCRIPT")
    if requested:
        path = shutil.which(requested) or requested
        if not Path(path).is_file():
            raise FileNotFoundError(f"LIUF_RSCRIPT does not exist: {path}")
        return validate_rscript_runtime(str(Path(path).resolve()))

    preferred = Path("/usr/local/bin/Rscript")
    if preferred.is_file():
        return validate_rscript_runtime(str(preferred))

    discovered = shutil.which("Rscript")
    if discovered is None:
        raise FileNotFoundError("Rscript was not found")
    return validate_rscript_runtime(str(Path(discovered).resolve()))
