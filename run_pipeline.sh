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
#   bash run_pipeline.sh --config my.sh     # Use custom config
#
# Steps:
#   1. STAR mapping
#   2. featureCounts (gene-level)
#   3. featureCounts (exon-level)
#   4. DESeq2 differential expression
#   5. fGSEA gene set enrichment
#   6. Exon-level DE (DEXSeq, splicing variants)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

# Load utility functions
source "${SCRIPT_DIR}/scripts/utils.sh"

CONFIG="${SCRIPT_DIR}/config.sh"
FROM_STEP=1
ONLY_STEP=0

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
            shift 2
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

source "${CONFIG}"

# --- Pre-flight checks ---
log_info "Pipeline version: $(get_pipeline_version)"
validate_samples "${CONFIG}"
acquire_lock "${PROJECT_DIR}"
record_provenance "${CONFIG}" "${PROJECT_DIR}"

echo "####################################################"
echo "# RNA-seq Analysis Pipeline  v$(get_pipeline_version)"
echo "# Project : ${PROJECT_DIR}"
echo "# Genome  : ${GENOME}"
echo "# Config  : ${CONFIG}"
echo "# Started : $(date)"
echo "####################################################"
echo ""

should_run() {
    local step=$1
    if [[ ${ONLY_STEP} -ne 0 ]]; then
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
    Rscript "${SCRIPT_DIR}/scripts/04_deseq2.R" "${CONFIG}"
    echo ""
fi

# Step 5: fGSEA
if should_run 5; then
    echo ">>> Step 5: fGSEA"
    cd "${PROJECT_DIR}"
    Rscript "${SCRIPT_DIR}/scripts/05_fgsea.R" "${CONFIG}"
    echo ""
fi

# Step 6: Exon-level DE (DEXSeq)
if should_run 6; then
    echo ">>> Step 6: Exon-level DE (DEXSeq)"
    cd "${PROJECT_DIR}"
    Rscript "${SCRIPT_DIR}/scripts/06_exon_de.R" "${CONFIG}"
    echo ""
fi

echo "####################################################"
echo "# Pipeline Complete: $(date)"
echo "####################################################"

release_lock "${PROJECT_DIR}"
log_info "All done."
