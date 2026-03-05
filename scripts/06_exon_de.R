#!/usr/bin/env Rscript
# =============================================================================
# Step 6: Exon-level Differential Expression (Splicing Variant Analysis)
# =============================================================================
# Uses DEXSeq to detect differential exon usage between groups,
# which indicates splicing variant differences.
#
# Usage:
#   Rscript scripts/06_exon_de.R [config.sh]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config.sh")

# --- Parse config ---
parse_config <- function(config_path) {
  lines <- readLines(config_path)
  lines <- lines[!grepl("^\\s*#", lines) & !grepl("^\\s*$", lines)]
  lines <- lines[!grepl("^(if|elif|else|fi|source|set|SCRIPT_DIR|PROJECT_DIR)", lines)]
  env <- list()
  for (line in lines) {
    m <- regmatches(line, regexpr("^([A-Za-z_][A-Za-z0-9_]*)=[\"']?([^\"']*)[\"']?", line))
    if (length(m) == 1 && nchar(m) > 0) {
      parts <- strsplit(m, "=", fixed = TRUE)[[1]]
      key <- parts[1]
      val <- paste(parts[-1], collapse = "=")
      val <- gsub('^["\']|["\']$', '', val)
      env[[key]] <- val
    }
  }
  env
}

cfg <- parse_config(config_file)

project_dir <- getwd()
samples_tsv <- file.path(project_dir, "samples.tsv")
exon_count  <- file.path(project_dir, "counts", "exon_count.txt")
output_dir  <- file.path(project_dir, "exon_DE")

# --- Load packages ---
suppressPackageStartupMessages({
  library(DEXSeq)
  library(tidyverse)
  library(data.table)
  library(magrittr)
  library(ggrepel)
  library(BiocParallel)
})

set.seed(123)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================\n")
cat(" Step 6: Exon-level DE (DEXSeq)\n")
cat(" Exon counts: ", exon_count, "\n")
cat(" Output     : ", output_dir, "\n")
cat("============================================\n")

# --- Read sample info ---
sample_info <- fread(samples_tsv, sep = "\t", header = TRUE)
sample_info <- sample_info[!grepl("^#", sample_info$sample_name), ]
colnames(sample_info) <- c("sample_name", "group", "fq_prefix", "lane_suffix")

# --- Read exon count matrix ---
# Format from featureCounts -f: Geneid, Chr, Start, End, Strand, Length, counts...
exon_raw <- fread(exon_count, header = TRUE)

# Create exon ID: gene:chr:start-end
exon_annot <- exon_raw[, 1:6]
colnames(exon_annot) <- c("gene_id", "chr", "start", "end", "strand", "length")
exon_annot$exon_id <- paste0(exon_annot$gene_id, ":", exon_annot$chr, ":",
                             exon_annot$start, "-", exon_annot$end)

# Count data (columns 7+)
count_mat <- as.matrix(exon_raw[, 7:ncol(exon_raw)])
# Round fractional counts (from --fraction)
count_mat <- round(count_mat)

# Clean column names
clean_names <- colnames(count_mat)
clean_names <- gsub(".*[/\\\\]", "", clean_names)
clean_names <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", clean_names)
clean_names <- gsub("-", "_", clean_names)
colnames(count_mat) <- clean_names

# Match to samples
matched <- sample_info$sample_name[sample_info$sample_name %in% clean_names]
if (length(matched) == 0) {
  stop("ERROR: No sample names match exon count columns.")
}
count_mat   <- count_mat[, matched, drop = FALSE]
sample_info <- sample_info[match(matched, sample_info$sample_name), ]

# --- Create DEXSeqDataSet ---
# Prepare annotation
gene_ids <- exon_annot$gene_id
exon_ids <- paste0("E", seq_len(nrow(exon_annot)))  # simple exon numbering

sample_table <- data.frame(
  row.names = sample_info$sample_name,
  condition = factor(sample_info$group)
)

cat("Building DEXSeqDataSet...\n")

dxd <- DEXSeqDataSet(
  countData   = count_mat,
  sampleData  = sample_table,
  design      = ~ sample + exon + condition:exon,
  featureID   = exon_ids,
  groupID     = gene_ids
)

# --- Run DEXSeq ---
cat("Estimating size factors...\n")
dxd <- estimateSizeFactors(dxd)

cat("Estimating dispersions...\n")
BPPARAM <- MulticoreParam(workers = min(4, parallel::detectCores()))
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)

cat("Testing for differential exon usage...\n")
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)

cat("Estimating exon fold changes...\n")
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = BPPARAM)

# --- Extract results ---
dxr <- DEXSeqResults(dxd)
dxr_df <- as.data.frame(dxr)

# Add original annotation
dxr_df$gene_name   <- gene_ids
dxr_df$exon_coord  <- exon_annot$exon_id
dxr_df$chr         <- exon_annot$chr
dxr_df$exon_start  <- exon_annot$start
dxr_df$exon_end    <- exon_annot$end

# Save full results
fwrite(dxr_df, file.path(output_dir, "DEXSeq_all_results.tsv"), sep = "\t")

# Filter significant
sig_results <- dxr_df %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(padj)
fwrite(sig_results, file.path(output_dir, "DEXSeq_significant_exons.tsv"), sep = "\t")

cat("Total exons tested : ", nrow(dxr_df), "\n")
cat("Significant exons  : ", nrow(sig_results), "\n")

# --- Summary per gene ---
gene_summary <- dxr_df %>%
  group_by(gene_name) %>%
  summarise(
    n_exons_tested = n(),
    n_sig_exons    = sum(!is.na(padj) & padj < 0.05),
    min_padj       = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_sig_exons > 0) %>%
  arrange(min_padj)

fwrite(gene_summary, file.path(output_dir, "DEXSeq_gene_summary.tsv"), sep = "\t")
cat("Genes with differential exon usage: ", nrow(gene_summary), "\n")

# --- Visualization: top genes ---
top_genes <- head(gene_summary$gene_name, 20)
if (length(top_genes) > 0) {
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, showWarnings = FALSE)

  for (g in top_genes) {
    tryCatch({
      png(file.path(plot_dir, paste0(g, "_exon_usage.png")),
          width = 1200, height = 600, res = 150)
      plotDEXSeq(dxr, g, legend = TRUE, cex.axis = 1.2, cex = 1.3,
                 lwd = 2, displayTranscripts = FALSE)
      dev.off()
    }, error = function(e) {
      message("  Could not plot: ", g, " - ", e$message)
    })
  }
}

cat("============================================\n")
cat(" Step 6: Exon-level DE Complete\n")
cat("============================================\n")
