#!/bin/bash
# =============================================================================
# Step 3: featureCounts (Exon-level)
# =============================================================================
# Counts reads per exon for downstream splicing variant analysis (DEXSeq).
# Output: counts/exon_count.txt
#
# Usage:
#   bash scripts/03_featurecounts_exon.sh
#   bash scripts/03_featurecounts_exon.sh my_config.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
CONFIG="${1:-${PROJECT_DIR}/config.sh}"
source "${CONFIG}"

echo "============================================"
echo " Step 3: featureCounts (Exon-level)"
echo " GTF: ${GTF}"
echo "============================================"

if [[ ! -f "${GTF}" ]]; then
    echo "ERROR: GTF file not found: ${GTF}" >&2
    exit 1
fi

mkdir -p "${COUNTS_DIR}"

# Collect BAM files
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

FC_OPTS=()
if [[ "${FEATURE_COUNTS_PAIRED}" == "true" ]]; then
    FC_OPTS+=("-p")
fi
if [[ "${FEATURE_COUNTS_MULTI}" == "true" ]]; then
    FC_OPTS+=("-M")
fi

EXON_COUNT_LONG="${COUNTS_DIR}/exon_count_long.txt"
EXON_COUNT_OUT="${COUNTS_DIR}/exon_count.txt"

echo "Running featureCounts (exon-level, -f flat) on ${#BAM_FILES[@]} BAM files..."

# -f: count at feature (exon) level instead of meta-feature (gene)
# -O: allow reads overlapping multiple exons
# --fraction: fractional counting for multi-overlapping reads
featureCounts \
    "${FC_OPTS[@]}" \
    -f \
    -O \
    --fraction \
    -a "${GTF}" \
    -g gene_name \
    -t exon \
    -o "${EXON_COUNT_LONG}" \
    "${BAM_FILES[@]}"

# Keep all annotation columns + counts for downstream DEXSeq
# Columns: Geneid, Chr, Start, End, Strand, Length, count_1, count_2, ...
awk -F '\t' 'BEGIN{OFS="\t"} NR==1{next} {print}' \
    "${EXON_COUNT_LONG}" > "${EXON_COUNT_OUT}"

echo "Exon-level count matrix: ${EXON_COUNT_OUT}"
echo "============================================"
echo " Step 3: featureCounts (Exon-level) Complete"
echo "============================================"
