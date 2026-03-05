#!/bin/bash
# =============================================================================
# Utility Functions for RNA-seq Pipeline
# =============================================================================
# Provides: logging, provenance recording, sample validation, run locking.
# Source this file at the beginning of run_pipeline.sh.
# =============================================================================

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*"; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  $*" >&2; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; }

# ---------------------------------------------------------------------------
# Pipeline version
# ---------------------------------------------------------------------------
get_pipeline_version() {
    local version_file
    version_file="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/VERSION"
    if [[ -f "${version_file}" ]]; then
        cat "${version_file}" | tr -d '[:space:]'
    else
        echo "unknown"
    fi
}

# ---------------------------------------------------------------------------
# Provenance recording
# ---------------------------------------------------------------------------
# Writes a provenance.yml file in the project directory for reproducibility.
# Call this at the start of each pipeline run.
#
# Usage: record_provenance <config_path> <project_dir>
# ---------------------------------------------------------------------------
record_provenance() {
    local config_path="${1}"
    local project_dir="${2}"
    local provenance_file="${project_dir}/provenance.yml"
    local pipeline_version
    pipeline_version="$(get_pipeline_version)"
    local pipeline_dir
    pipeline_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

    # Pipeline git info (if available)
    local pipeline_commit="N/A"
    local pipeline_branch="N/A"
    if command -v git &>/dev/null && [[ -d "${pipeline_dir}/.git" ]]; then
        pipeline_commit="$(git -C "${pipeline_dir}" rev-parse --short HEAD 2>/dev/null || echo 'N/A')"
        pipeline_branch="$(git -C "${pipeline_dir}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'N/A')"
    fi

    cat > "${provenance_file}" <<EOF
# ============================================================
# Pipeline Provenance — auto-generated, do not edit manually
# ============================================================
pipeline:
  version: "${pipeline_version}"
  commit: "${pipeline_commit}"
  branch: "${pipeline_branch}"
  path: "${pipeline_dir}"

run:
  date: "$(date '+%Y-%m-%d %H:%M:%S')"
  user: "$(whoami)"
  hostname: "$(hostname)"
  config: "${config_path}"

environment:
  os: "$(uname -srm)"
  bash: "${BASH_VERSION}"
  star: "$(STAR --version 2>/dev/null || echo 'not found')"
  featureCounts: "$(featureCounts -v 2>&1 | head -1 || echo 'not found')"
  R: "$(R --version 2>/dev/null | head -1 || echo 'not found')"
EOF
    log_info "Provenance recorded → ${provenance_file}"
}

# ---------------------------------------------------------------------------
# Sample validation
# ---------------------------------------------------------------------------
# Validates that samples.tsv exists and has the required columns.
# Also checks that FASTQ files exist for each sample.
#
# Usage: validate_samples <config_path>
# ---------------------------------------------------------------------------
validate_samples() {
    local config_path="${1}"
    source "${config_path}"

    if [[ ! -f "${SAMPLES_TSV}" ]]; then
        log_error "samples.tsv not found: ${SAMPLES_TSV}"
        return 1
    fi

    # Check header
    local header
    header="$(head -1 "${SAMPLES_TSV}")"
    for col in sample_name group fq_prefix lane_suffix; do
        if ! echo "${header}" | grep -qw "${col}"; then
            log_error "Missing column '${col}' in ${SAMPLES_TSV}"
            return 1
        fi
    done

    # Check for duplicate sample names
    local dup_count
    dup_count="$(tail -n +2 "${SAMPLES_TSV}" | awk -F'\t' '{print $1}' | sort | uniq -d | wc -l)"
    if [[ "${dup_count}" -gt 0 ]]; then
        log_error "Duplicate sample names found in ${SAMPLES_TSV}:"
        tail -n +2 "${SAMPLES_TSV}" | awk -F'\t' '{print $1}' | sort | uniq -d >&2
        return 1
    fi

    # Check that at least 2 groups exist
    local group_count
    group_count="$(tail -n +2 "${SAMPLES_TSV}" | awk -F'\t' '{print $2}' | sort -u | wc -l)"
    if [[ "${group_count}" -lt 2 ]]; then
        log_error "At least 2 groups required in ${SAMPLES_TSV}, found ${group_count}"
        return 1
    fi

    # Check FASTQ files exist (optional — warn only)
    local missing=0
    while IFS=$'\t' read -r sname group prefix suffix; do
        [[ "${sname}" == "sample_name" ]] && continue   # skip header
        local r1="${FASTQ_DIR}/${prefix}_${suffix}_1.fq.gz"
        if [[ ! -f "${r1}" ]]; then
            log_warn "FASTQ not found: ${r1} (sample: ${sname})"
            missing=$((missing + 1))
        fi
    done < "${SAMPLES_TSV}"

    if [[ ${missing} -gt 0 ]]; then
        log_warn "${missing} FASTQ file(s) missing — mapping step will fail for those samples"
    fi

    local sample_count
    sample_count="$(tail -n +2 "${SAMPLES_TSV}" | wc -l)"
    log_info "Validated ${sample_count} samples in ${group_count} groups"
    return 0
}

# ---------------------------------------------------------------------------
# Run locking
# ---------------------------------------------------------------------------
# Prevents concurrent pipeline runs in the same project directory.
#
# Usage:
#   acquire_lock <project_dir>     # Call at the start
#   release_lock <project_dir>     # Called automatically via trap
# ---------------------------------------------------------------------------
LOCK_FILE=""

acquire_lock() {
    local project_dir="${1}"
    LOCK_FILE="${project_dir}/.pipeline.lock"

    if [[ -f "${LOCK_FILE}" ]]; then
        local lock_pid lock_time
        lock_pid="$(head -1 "${LOCK_FILE}" 2>/dev/null || echo '?')"
        lock_time="$(tail -1 "${LOCK_FILE}" 2>/dev/null || echo '?')"

        # Check if the locking process is still running
        if kill -0 "${lock_pid}" 2>/dev/null; then
            log_error "Pipeline already running (PID ${lock_pid}, started ${lock_time})"
            log_error "If this is stale, remove ${LOCK_FILE} manually."
            return 1
        else
            log_warn "Stale lock found (PID ${lock_pid} not running). Removing."
            rm -f "${LOCK_FILE}"
        fi
    fi

    echo "$$" > "${LOCK_FILE}"
    echo "$(date '+%Y-%m-%d %H:%M:%S')" >> "${LOCK_FILE}"
    trap 'release_lock "${LOCK_FILE%/*}"' EXIT INT TERM
    log_info "Lock acquired (PID $$)"
}

release_lock() {
    local project_dir="${1}"
    local lock="${project_dir}/.pipeline.lock"
    if [[ -f "${lock}" ]]; then
        rm -f "${lock}"
        log_info "Lock released"
    fi
}
