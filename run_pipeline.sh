#!/bin/bash
# =============================================================================
# RNA-seq Analysis Pipeline
# =============================================================================
# Main entry point that runs all pipeline steps sequentially.
# Each step can also be run independently.
#
# Usage:
#   bash run_pipeline.sh                    # Run all steps
#   bash run_pipeline.sh --from 3           # Start from step 3
#   bash run_pipeline.sh --only 4           # Run only step 4
#   bash run_pipeline.sh --only 4,5         # Run only steps 4 and 5
#   bash run_pipeline.sh --config my.sh     # Use custom config
#   bash run_pipeline.sh --from 4 \
#       --merge-samples a/samples.tsv,b/samples.tsv \
#       --merge-counts a/counts/count.txt,b/counts/count.txt \
#       --merge-exon-counts a/counts/exon_count.txt,b/counts/exon_count.txt
#
# Steps:
#   1. STAR mapping
#   2. featureCounts (gene-level)
#   3. featureCounts (exon-level)
#   4. DESeq2 differential expression
#   5. fGSEA gene set enrichment
#   6. ssGSEA and selected-gene visualization
#   7. Exon-level DE (DEXSeq, splicing variants)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

# Load utility functions
source "${SCRIPT_DIR}/scripts/utils.sh"

CONFIG="${SCRIPT_DIR}/config.sh"
FROM_STEP=1
ONLY_STEP=0
ONLY_STEPS=()
FORCE="false"
MERGE_SAMPLES=""
MERGE_COUNTS=""
MERGE_EXON_COUNTS=""

parse_step_list() {
    local raw="${1}"
    local item
    IFS=',' read -ra step_items <<< "${raw}"
    for item in "${step_items[@]}"; do
        item="$(echo "${item}" | xargs)"
        if [[ ! "${item}" =~ ^[1-7]$ ]]; then
            log_error "Invalid step in --only: ${item} (allowed: 1-7)"
            exit 1
        fi
        ONLY_STEPS+=("${item}")
    done
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG="$2"
            shift 2
            ;;
        --from)
            FROM_STEP="$2"
            shift 2
            ;;
        --only)
            ONLY_STEP="$2"
            parse_step_list "$2"
            shift 2
            ;;
        --force)
            FORCE="true"
            shift
            ;;
        --merge-samples)
            MERGE_SAMPLES="$2"
            shift 2
            ;;
        --merge-counts)
            MERGE_COUNTS="$2"
            shift 2
            ;;
        --merge-exon-counts)
            MERGE_EXON_COUNTS="$2"
            shift 2
            ;;
        --version)
            echo "RNAseq_pipeline_takubo v$(get_pipeline_version)"
            exit 0
            ;;
        -h|--help)
            head -25 "$0" | tail -20
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# --- Conda 環境チェック ---
EXPECTED_ENV="rnaseq_takubo"
current_env="${CONDA_DEFAULT_ENV:-}"
if [[ "${current_env}" != "${EXPECTED_ENV}" ]]; then
    echo -e "\033[1;33m[WARN]\033[0m Conda環境 '${EXPECTED_ENV}' がアクティブではありません。"
    echo -e "       現在の環境: '${current_env:-none}'"
    echo ""
    echo "  セットアップ:   bash setup_env.sh"
    echo "  アクティベート: conda activate ${EXPECTED_ENV}"
    echo ""
    read -rp "このまま続行しますか？ [y/N]: " yn
    case "${yn}" in
        [Yy]*) echo "[INFO] 現在の環境で続行します。" ;;
        *)     echo "[INFO] 中断しました。"; exit 1 ;;
    esac
fi

source "${CONFIG}"

resolve_rscript() {
    if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/Rscript" ]]; then
        echo "${CONDA_PREFIX}/bin/Rscript"
        return 0
    fi

    if command -v Rscript >/dev/null 2>&1; then
        command -v Rscript
        return 0
    fi

    return 1
}

RSCRIPT_BIN="$(resolve_rscript || true)"
if [[ -z "${RSCRIPT_BIN}" ]]; then
    log_error "Rscript was not found. Activate the environment and try again."
    release_lock "${PROJECT_DIR}"
    exit 1
fi

# --- Pre-flight checks ---
log_info "Pipeline version: $(get_pipeline_version)"
log_info "Using Rscript: ${RSCRIPT_BIN}"
validate_samples "${CONFIG}"
acquire_lock "${PROJECT_DIR}"
record_provenance "${CONFIG}" "${PROJECT_DIR}"

# Export PIPELINE_ROOT for R scripts (plot_utils.R)
export PIPELINE_ROOT="${SCRIPT_DIR}"
export PIPELINE_MERGE_SAMPLES="${MERGE_SAMPLES}"
export PIPELINE_MERGE_COUNTS="${MERGE_COUNTS}"
export PIPELINE_MERGE_EXON_COUNTS="${MERGE_EXON_COUNTS}"

echo "####################################################"
echo "# RNA-seq Analysis Pipeline  v$(get_pipeline_version)"
echo "# Project : ${PROJECT_DIR}"
echo "# Genome  : ${GENOME}"
echo "# Config  : ${CONFIG}"
echo "# Conda   : ${CONDA_DEFAULT_ENV:-none}"
if [[ -n "${MERGE_SAMPLES}" ]]; then
    echo "# Merge samples : ${MERGE_SAMPLES}"
fi
if [[ -n "${MERGE_COUNTS}" ]]; then
    echo "# Merge counts  : ${MERGE_COUNTS}"
fi
if [[ -n "${MERGE_EXON_COUNTS}" ]]; then
    echo "# Merge exon    : ${MERGE_EXON_COUNTS}"
fi
echo "# Started : $(date)"
echo "####################################################"
echo ""

should_run() {
    local step=$1
    if [[ ${#ONLY_STEPS[@]} -gt 0 ]]; then
        local only_step
        for only_step in "${ONLY_STEPS[@]}"; do
            if [[ ${step} -eq ${only_step} ]]; then
                return 0
            fi
        done
        return 1
    elif [[ ${ONLY_STEP} -ne 0 ]]; then
        [[ ${step} -eq ${ONLY_STEP} ]]
    else
        [[ ${step} -ge ${FROM_STEP} ]]
    fi
}

# Step 1: STAR mapping
if should_run 1; then
    echo ">>> Step 1: STAR Mapping"
    bash "${SCRIPT_DIR}/scripts/01_mapping.sh" "${CONFIG}"
    echo ""
fi

# Step 2: featureCounts (gene-level)
if should_run 2; then
    echo ">>> Step 2: featureCounts (Gene-level)"
    bash "${SCRIPT_DIR}/scripts/02_featurecounts.sh" "${CONFIG}"
    echo ""
fi

# Step 3: featureCounts (exon-level)
if should_run 3; then
    echo ">>> Step 3: featureCounts (Exon-level)"
    bash "${SCRIPT_DIR}/scripts/03_featurecounts_exon.sh" "${CONFIG}"
    echo ""
fi

# Step 4: DESeq2
if should_run 4; then
    echo ">>> Step 4: DESeq2"
    cd "${PROJECT_DIR}"
    "${RSCRIPT_BIN}" "${SCRIPT_DIR}/scripts/04_deseq2.R" "${CONFIG}"
    echo ""
fi

# Step 5: fGSEA
if should_run 5; then
    echo ">>> Step 5: fGSEA"
    cd "${PROJECT_DIR}"
    "${RSCRIPT_BIN}" "${SCRIPT_DIR}/scripts/05_fgsea.R" "${CONFIG}"
    echo ""
fi

# Step 6: ssGSEA and selected-gene visualization
if should_run 6; then
    echo ">>> Step 6: ssGSEA and selected-gene visualization"
    cd "${PROJECT_DIR}"
    "${RSCRIPT_BIN}" "${SCRIPT_DIR}/scripts/06_ssgsea_heatmap.R" "${CONFIG}"
    echo ""
fi

# Step 7: Exon-level DE (DEXSeq)
if should_run 7; then
    echo ">>> Step 7: Exon-level DE (DEXSeq)"
    cd "${PROJECT_DIR}"
    "${RSCRIPT_BIN}" "${SCRIPT_DIR}/scripts/07_exon_de.R" "${CONFIG}"
    echo ""
fi

echo "####################################################"
echo "# Pipeline Complete: $(date)"
echo "####################################################"

release_lock "${PROJECT_DIR}"
log_info "All done."
