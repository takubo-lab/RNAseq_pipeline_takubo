#!/bin/bash
# =============================================================================
# RNA-seq Pipeline Configuration
# =============================================================================
# Edit this file to configure the pipeline for your experiment.
# All paths should be absolute or relative to the project directory.
# =============================================================================

# --- Pipeline location ---
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Project ---
# When using init_project.sh, PROJECT_DIR points to the project directory.
# For standalone use, it defaults to the pipeline directory itself.
PROJECT_DIR="${PIPELINE_DIR}"
SAMPLES_TSV="${PROJECT_DIR}/samples.tsv"

# --- Reference Genome ---
# Supported: "hg38" or "mm10"
GENOME="hg38"

# STAR genome index directory
# Adjust these paths to match your local STAR index locations
STAR_INDEX_HG38="${HOME}/STAR2.7.3/hg38_sjdb149_RefSeq"
STAR_INDEX_MM10="${HOME}/STAR2.7.3/mm10_sjdb149_RefSeq_ver241129"

# GTF annotation file
GTF_HG38="${HOME}/Genome/hg38/hg38.ncbiRefSeq.gtf"
GTF_MM10="${HOME}/Genome/mm10/mm10.ncbiRefSeq.gtf"

# --- STAR mapping parameters ---
STAR_THREADS=6
SJDB_OVERHANG=149
READ_MAP_NUMBER=-1       # -1 = map all reads

# --- featureCounts parameters ---
FEATURE_COUNTS_PAIRED=true   # true for paired-end, false for single-end
FEATURE_COUNTS_MULTI=true    # true to count multi-mapping reads (-M)
GENE_ATTRIBUTE="gene_name"   # attribute used for gene-level counting
EXON_ATTRIBUTE="exon_id"     # attribute used for exon-level counting

# --- FASTQ settings ---
# Pattern for finding FASTQ files: {SAMPLE_PREFIX}{LANE_SUFFIX}_{1,2}.fq.gz
# LANE_SUFFIX is specified per-sample in samples.tsv
FASTQ_DIR="${PROJECT_DIR}/fastq"

# --- Output directories ---
BAM_DIR="${PROJECT_DIR}/BAM"
COUNTS_DIR="${PROJECT_DIR}/counts"
DESEQ2_DIR="${PROJECT_DIR}/DESeq2"
FGSEA_DIR="${PROJECT_DIR}/fGSEA"
EXON_DE_DIR="${PROJECT_DIR}/exon_DE"

# --- DESeq2 parameters ---
MIN_COUNT_FILTER=50          # minimum total count across samples to keep a gene
PCA_TOP_GENES=500            # number of top variable genes for PCA
VOLCANO_PADJ=0.001           # adjusted p-value threshold for volcano plots
VOLCANO_LOG2FC=1             # log2 fold change threshold for volcano plots

# --- fGSEA parameters ---
FGSEA_MIN_SIZE=15
FGSEA_MAX_SIZE=500
FGSEA_NPERM=10000
FGSEA_NPROC=12
# Species for MSigDB: "Homo sapiens" or "Mus musculus"
# Auto-set based on GENOME
FGSEA_PADJ=0.05

# --- Gene set files (optional, in addition to MSigDB) ---
CUSTOM_GMT="${PROJECT_DIR}/Gene_set/Human_old_HSC_set2.gmt"

# --- Single-sample GSEA / selected-gene visualization ---
SSGSEA_OUTPUT_DIR="${PROJECT_DIR}/single_sample"
SSGSEA_MSIG_CATEGORIES="H,C2,C3,C4,C5,C6,C7"
SSGSEA_MIN_SIZE=10
SSGSEA_MAX_SIZE=5000
SSGSEA_MIN_EXPR=1
SSGSEA_MIN_SAMPLES_FRACTION=0.5

# Groups for per-gene plots. Leave empty to auto-detect when exactly 2 groups exist.
GENE_PLOT_GROUPS=""
GENE_PLOT_PAIRED=true

# --- Genes to highlight in volcano plots (comma-separated) ---
HIGHLIGHT_GENES="ITGA3,PROCR,KIT,HLF,MECOM,MYCT1,MLLT3,THY1,ADGRG1,CD34"

# =============================================================================
# Derived values (do not edit below this line)
# =============================================================================
if [[ "${GENOME}" == "hg38" ]]; then
    STAR_INDEX="${STAR_INDEX_HG38}"
    GTF="${GTF_HG38}"
    SPECIES="Homo sapiens"
    ORG_DB="org.Hs.eg.db"
elif [[ "${GENOME}" == "mm10" ]]; then
    STAR_INDEX="${STAR_INDEX_MM10}"
    GTF="${GTF_MM10}"
    SPECIES="Mus musculus"
    ORG_DB="org.Mm.eg.db"
else
    echo "ERROR: Unsupported GENOME '${GENOME}'. Use 'hg38' or 'mm10'." >&2
    exit 1
fi
