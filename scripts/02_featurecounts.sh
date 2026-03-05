#!/bin/bash
# =============================================================================
# Step 2: featureCounts (Gene-level)
# =============================================================================
# Counts reads per gene using featureCounts from Subread package.
# Output: counts/count.txt (gene-level count matrix)
#
# Usage:
#   bash scripts/02_featurecounts.sh
#   bash scripts/02_featurecounts.sh my_config.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
CONFIG="${1:-${PROJECT_DIR}/config.sh}"
source "${CONFIG}"

echo "============================================"
echo " Step 2: featureCounts (Gene-level)"
echo " GTF: ${GTF}"
echo "============================================"

if [[ ! -f "${GTF}" ]]; then
    echo "ERROR: GTF file not found: ${GTF}" >&2
    exit 1
fi

mkdir -p "${COUNTS_DIR}"

# Collect all BAM files
BAM_FILES=()
while IFS=$'\t' read -r sample_name group fq_prefix lane_suffix; do
    [[ "${sample_name}" =~ ^#.*$ || -z "${sample_name}" ]] && continue
    BAM="${BAM_DIR}/${sample_name}_Aligned.sortedByCoord.out.bam"
    if [[ -f "${BAM}" ]]; then
        BAM_FILES+=("${BAM}")
    else
        echo "WARNING: BAM not found: ${BAM}" >&2
    fi
done < "${SAMPLES_TSV}"

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found." >&2
    exit 1
fi

# Build featureCounts options
FC_OPTS=()
if [[ "${FEATURE_COUNTS_PAIRED}" == "true" ]]; then
    FC_OPTS+=("-p")
fi
if [[ "${FEATURE_COUNTS_MULTI}" == "true" ]]; then
    FC_OPTS+=("-M")
fi

COUNT_LONG="${COUNTS_DIR}/count_long.txt"
COUNT_OUT="${COUNTS_DIR}/count.txt"

echo "Running featureCounts on ${#BAM_FILES[@]} BAM files..."

featureCounts \
    "${FC_OPTS[@]}" \
    -a "${GTF}" \
    -g "${GENE_ATTRIBUTE}" \
    -t exon \
    -o "${COUNT_LONG}" \
    "${BAM_FILES[@]}"

# Extract Geneid + count columns (skip first header comment line, remove columns 2-6)
# Column 1 = Geneid, Columns 2-6 = Chr/Start/End/Strand/Length, Columns 7+ = counts
awk -F '\t' 'BEGIN{OFS="\t"} NR==1{next} {printf $1; for(i=7;i<=NF;i++) printf "\t"$i; print ""}' \
    "${COUNT_LONG}" > "${COUNT_OUT}"

echo "Gene-level count matrix: ${COUNT_OUT}"
echo "============================================"
echo " Step 2: featureCounts (Gene-level) Complete"
echo "============================================"
