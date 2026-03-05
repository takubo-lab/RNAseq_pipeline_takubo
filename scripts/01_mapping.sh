#!/bin/bash
# =============================================================================
# Step 1: STAR Mapping
# =============================================================================
# Maps paired-end FASTQ files to the reference genome using STAR.
# Reads sample information from samples.tsv and parameters from config.sh.
#
# Usage:
#   bash scripts/01_mapping.sh          # Uses config.sh in project root
#   bash scripts/01_mapping.sh my.sh    # Uses custom config file
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
CONFIG="${1:-${PROJECT_DIR}/config.sh}"
source "${CONFIG}"

echo "============================================"
echo " Step 1: STAR Mapping"
echo " Genome : ${GENOME}"
echo " Index  : ${STAR_INDEX}"
echo "============================================"

# Validate
if [[ ! -d "${STAR_INDEX}" ]]; then
    echo "ERROR: STAR index not found: ${STAR_INDEX}" >&2
    exit 1
fi

mkdir -p "${BAM_DIR}"

# Read samples.tsv (skip header and comments)
while IFS=$'\t' read -r sample_name group fq_prefix lane_suffix; do
    [[ "${sample_name}" =~ ^#.*$ || -z "${sample_name}" ]] && continue

    R1="${FASTQ_DIR}/${fq_prefix}_${lane_suffix}_1.fq.gz"
    R2="${FASTQ_DIR}/${fq_prefix}_${lane_suffix}_2.fq.gz"

    # Check if FASTQ files exist
    if [[ ! -f "${R1}" ]]; then
        echo "WARNING: R1 not found: ${R1}, skipping ${sample_name}" >&2
        continue
    fi

    # Check if output BAM already exists
    BAM_OUT="${BAM_DIR}/${sample_name}_Aligned.sortedByCoord.out.bam"
    if [[ -f "${BAM_OUT}" ]]; then
        echo "SKIP: BAM already exists for ${sample_name}"
        continue
    fi

    echo "Mapping: ${sample_name} ..."

    if [[ "${FEATURE_COUNTS_PAIRED}" == "true" ]]; then
        # Paired-end
        if [[ ! -f "${R2}" ]]; then
            echo "WARNING: R2 not found: ${R2}, skipping ${sample_name}" >&2
            continue
        fi
        STAR \
            --runThreadN "${STAR_THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesCommand zcat \
            --readFilesIn "${R1}" "${R2}" \
            --readMapNumber "${READ_MAP_NUMBER}" \
            --sjdbOverhang "${SJDB_OVERHANG}" \
            --outFileNamePrefix "${BAM_DIR}/${sample_name}_" \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMsortingThreadN "${STAR_THREADS}"
    else
        # Single-end
        STAR \
            --runThreadN "${STAR_THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesCommand zcat \
            --readFilesIn "${R1}" \
            --readMapNumber "${READ_MAP_NUMBER}" \
            --sjdbOverhang "${SJDB_OVERHANG}" \
            --outFileNamePrefix "${BAM_DIR}/${sample_name}_" \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMsortingThreadN "${STAR_THREADS}"
    fi

    echo "Done: ${sample_name}"
done < "${SAMPLES_TSV}"

echo "============================================"
echo " Step 1: STAR Mapping Complete"
echo "============================================"
