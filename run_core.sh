#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"

usage() {
  cat <<'EOF'
Usage: bash run_core.sh <stage>

Stages:
  preprocess       Rebuild the eight tissue/serum M0 matrices
  scrna_atlas      Rebuild the raw-derived scRNA atlas and subtype objects
  covarnet         Validate curated inputs and rebuild K=5 CoVarNet/CM02
  bulk_projection  Run MuSiC, CM02 trajectory, survival and Hep_C activities
  atlas_tables     Compute composition and ratio-of-enrichment tables
  tissue           Run PROGENy, PTM-SEA, cross-layer and tissue statistics
  ptmsea_rng       Run the PTM-SEA seed-mode comparison
  resistance       Run DepMap/PRISM composites and STORM limma/GSEA
  spatial          Run spatial recovery, pathways, COMMOT/CCI and CORNETO
  serum_mofa       Run the serum MOFA downstream workflow

Set LIUF_SPATIAL_RUN_VARIANT=raw_rctd_diagnostic to refit all 16 RCTD sections.
Set MOFA_ANALYSIS_SOURCE=raw_derived, MOFA_RUN_RAW_DIAGNOSTIC=1, and
MOFA_FORCE_RETRAIN=1 to retrain the MOFA grid.
EOF
}

if [[ "${1:---list}" == "--list" || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

source "${PROJECT_ROOT}/code/lib/shell/repro_env.sh"
liuf_init_repro_env "${PROJECT_ROOT}"
cd "${PROJECT_ROOT}"

R="${LIUF_RSCRIPT}"
PYTHON="${LIUF_PYTHON}"
if [[ -n "${LIUF_C2L_PYTHON:-}" && ! -x "${LIUF_C2L_PYTHON}" ]]; then
  printf 'LIUF_C2L_PYTHON is not executable: %s\n' "${LIUF_C2L_PYTHON}" >&2
  exit 1
fi
SPATIAL_PYTHON="${LIUF_C2L_PYTHON:-${HOME}/miniconda3/envs/c2l/bin/python}"
if [[ ! -x "${SPATIAL_PYTHON}" ]]; then
  SPATIAL_PYTHON="python3"
fi

run_r() {
  local script="$1"
  shift
  printf '[core] %s %s\n' "${R}" "${script} $*"
  "${R}" "${PROJECT_ROOT}/${script}" "$@"
}

run_py() {
  local script="$1"
  shift
  printf '[core] %s %s\n' "${PYTHON}" "${script} $*"
  "${PYTHON}" "${PROJECT_ROOT}/${script}" "$@"
}

run_spatial_py() {
  local script="$1"
  shift
  printf '[core] %s %s\n' "${SPATIAL_PYTHON}" "${script} $*"
  "${SPATIAL_PYTHON}" "${PROJECT_ROOT}/${script}" "$@"
}

verify_spatial_runtime() {
  "${SPATIAL_PYTHON}" - <<'PY'
from importlib import import_module
from importlib.metadata import version

required = {
    "scanpy": "1.11.5",
    "commot": "0.0.3",
    "corneto": "0.9.1a6",
}
for name, expected in required.items():
    try:
        module = import_module(name)
    except ImportError as error:
        raise SystemExit(f"Missing spatial dependency {name}: {error}")
    observed = version(name)
    if observed != expected:
        raise SystemExit(f"{name} {observed} does not match required {expected}")
    print(f"{name}={observed}")
PY
}

case "$1" in
  preprocess)
    bash "${PROJECT_ROOT}/code/scripts/10_run_raw_multiomics.sh"
    ;;
  scrna_atlas)
    run_r code/pipelines/20_scrna_cm02/03_rebuild_atlas_from_local_raw.R
    run_r code/pipelines/20_scrna_cm02/04_rebuild_atlas_panel_data.R
    ;;
  covarnet)
    run_r code/pipelines/20_scrna_cm02/02_audit_frozen_curation.R
    run_r code/pipelines/20_scrna_cm02/05_rebuild_covarnet_cm02.R
    ;;
  bulk_projection)
    export FIG3_DOWNSTREAM_SOURCE="raw_derived"
    run_r code/pipelines/20_scrna_cm02/27_build_music_reference.R
    for cohort in local TCGA GTEx; do
      run_r code/pipelines/20_scrna_cm02/28_run_music_rnaseq.R "${cohort}"
    done
    run_r code/pipelines/20_scrna_cm02/29_run_music_gse63898.R
    run_r code/pipelines/20_scrna_cm02/30_build_bulk_cm02_trajectory.R
    run_r code/pipelines/20_scrna_cm02/31_rebuild_hepc_activities.R
    run_r code/pipelines/20_scrna_cm02/11_replay_figure3c_tcga.R
    run_r code/pipelines/20_scrna_cm02/17_replay_supplemental2_bulk.R
    ;;
  atlas_tables)
    run_r code/pipelines/20_scrna_cm02/23_replay_figure2_composition_roe.R
    run_r code/pipelines/20_scrna_cm02/25_replay_supplemental1_composition.R
    ;;
  tissue)
    run_r code/pipelines/30_tissue_multiomics/00_validate_inputs.R
    run_r code/pipelines/30_tissue_multiomics/05_freeze_local_msigdb.R
    run_r code/pipelines/30_tissue_multiomics/10_progeny.R
    run_r code/pipelines/30_tissue_multiomics/20_ptmsea.R
    run_r code/pipelines/30_tissue_multiomics/30_cross_layer.R
    run_r code/pipelines/30_tissue_multiomics/40_ksea_three_layer.R
    run_r code/pipelines/30_tissue_multiomics/60_fig4c.R
    ;;
  ptmsea_rng)
    run_r code/pipelines/30_tissue_multiomics/22_run_ptmsea_rng_probe.R
    run_r code/pipelines/30_tissue_multiomics/23_summarize_ptmsea_rng_modes.R
    ;;
  resistance)
    run_r code/pipelines/40_tissue_resistance/30_rebuild_depmap_prism.R
    run_r code/pipelines/40_tissue_resistance/41_rebuild_storm_gsea.R
    ;;
  spatial)
    verify_spatial_runtime
    export LIUF_SPATIAL_RUN_VARIANT="${LIUF_SPATIAL_RUN_VARIANT:-article_equivalent}"
    run_spatial_py code/pipelines/50_spatial_cci/01_raw_visium_qc.py
    if [[ "${LIUF_SPATIAL_RUN_VARIANT}" == "raw_rctd_diagnostic" ]]; then
      run_r code/pipelines/50_spatial_cci/10_rctd_from_raw.R --force
    else
      run_spatial_py code/pipelines/50_spatial_cci/10_validate_historical_intermediates.py
    fi
    run_r code/pipelines/50_spatial_cci/20_recover_cm02.R
    run_r code/pipelines/50_spatial_cci/30_radial_from_raw.R
    run_r code/pipelines/50_spatial_cci/40_pathway_scores_from_raw.R
    run_r code/pipelines/50_spatial_cci/41_merge_cm02_pathways.R
    run_spatial_py code/pipelines/50_spatial_cci/49_cellchat_database.py
    run_r code/pipelines/50_spatial_cci/50_export_scrna_profiles.R
    run_spatial_py code/pipelines/50_spatial_cci/51_commot_from_raw.py --force
    run_spatial_py code/pipelines/50_spatial_cci/60_pairwise_commot.py
    run_spatial_py code/pipelines/50_spatial_cci/70_lr_overlap_river.py
    run_spatial_py code/pipelines/50_spatial_cci/80_corneto_mapk.py
    ;;
  serum_mofa)
    bash "${PROJECT_ROOT}/code/pipelines/60_serum_mofa/run_all.sh"
    ;;
  *)
    printf 'Unknown stage: %s\n\n' "$1" >&2
    usage >&2
    exit 2
    ;;
esac
