#!/usr/bin/env Rscript
# =============================================================================
# Step 6: Exon-level Differential Expression (Splicing Variant Analysis)
# =============================================================================
# Uses DEXSeq to detect differential exon usage between groups,
# which indicates splicing variant differences.
#
# Cross-platform:
#   Linux : bash run_pipeline.sh (full pipeline)
#   Win/Mac: Rscript scripts/06_exon_de.R [config.sh]
#            (exon_count.txt と samples.tsv が事前に用意されている前提)
#
# 可視化設定: plot_config.default.yml / plot_config.yml
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
    line <- sub("\\s+#.*$", "", line)
    m <- regmatches(line, regexpr("^([A-Za-z_][A-Za-z0-9_]*)=[\"']?([^\"']*)[\"']?", line))
    if (length(m) == 1 && nchar(m) > 0) {
      parts <- strsplit(m, "=", fixed = TRUE)[[1]]
      key <- parts[1]
      val <- paste(parts[-1], collapse = "=")
      val <- gsub('^["\047]|["\047]$', '', val)
      env[[key]] <- val
    }
  }
  env
}

`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_sample_name <- function(x) {
  gsub("-", "_", x)
}

resolve_input_path <- function(path_value, project_dir, cfg = list()) {
  if (is.null(path_value) || !nzchar(path_value)) {
    return("")
  }

  resolved <- trimws(path_value)
  replacements <- c(
    PROJECT_DIR = project_dir,
    PIPELINE_DIR = cfg$PIPELINE_DIR %||% project_dir,
    HOME = Sys.getenv("HOME", unset = "")
  )

  for (nm in names(replacements)) {
    resolved <- gsub(paste0("\\$\\{", nm, "\\}"), replacements[[nm]], resolved)
    resolved <- gsub(paste0("\\$", nm), replacements[[nm]], resolved)
  }

  path.expand(resolved)
}

parse_path_list <- function(x, project_dir, cfg = list()) {
  if (is.null(x) || !nzchar(x)) {
    return(character())
  }

  unique(vapply(strsplit(x, ",", fixed = TRUE)[[1]], function(p) {
    resolve_input_path(p, project_dir, cfg)
  }, character(1)))
}

read_sample_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  required_cols <- c("sample_name", "group", "fq_prefix", "lane_suffix")

  if (length(lines) == 0) {
    stop("ERROR: samples.tsv is empty: ", path)
  }

  normalized <- sub("^\\s*#\\s*", "", lines)
  header_idx <- which(vapply(normalized, function(line) {
    fields <- strsplit(trimws(line), "\\t")[[1]]
    length(fields) >= length(required_cols) && identical(fields[seq_along(required_cols)], required_cols)
  }, logical(1)))[1]

  if (!is.na(header_idx) && grepl("^\\s*#", lines[header_idx])) {
    lines[header_idx] <- normalized[header_idx]
  }

  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!grepl("^\\s*#", lines)]

  dat <- fread(
    text = paste(lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    blank.lines.skip = TRUE
  )

  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop("ERROR: Missing required columns in samples.tsv: ",
         paste(missing_cols, collapse = ", "),
         " (file: ", path, ")")
  }

  dat[, c("sample_name", "group", "fq_prefix", "lane_suffix"), with = FALSE]
}

read_sample_tables <- function(sample_files) {
  sample_list <- lapply(sample_files, function(path) {
    if (!file.exists(path)) {
      stop("ERROR: samples.tsv not found: ", path)
    }
    read_sample_table(path)
  })

  sample_info <- dplyr::bind_rows(sample_list)
  dup <- sample_info$sample_name[duplicated(sample_info$sample_name)]
  if (length(dup) > 0) {
    stop("ERROR: Duplicate sample names found across merged sample sheets: ",
         paste(unique(dup), collapse = ", "))
  }
  sample_info
}

read_exon_count_tables <- function(exon_files) {
  exon_list <- lapply(exon_files, function(path) {
    if (!file.exists(path)) {
      stop("ERROR: Exon count file not found: ", path)
    }
    fread(path, header = TRUE)
  })

  key_cols <- colnames(exon_list[[1]])[1:6]
  merged <- Reduce(function(x, y) dplyr::full_join(x, y, by = key_cols), exon_list)
  merged[is.na(merged)] <- 0

  dup_cols <- colnames(merged)[duplicated(colnames(merged))]
  if (length(dup_cols) > 0) {
    stop("ERROR: Duplicate exon count columns found across merged exon count files: ",
         paste(unique(dup_cols), collapse = ", "))
  }

  merged
}

cfg <- parse_config(config_file)

project_dir <- getwd()
output_dir  <- file.path(project_dir, "exon_DE")

merge_samples_env <- Sys.getenv("PIPELINE_MERGE_SAMPLES", unset = "")
merge_exon_env    <- Sys.getenv("PIPELINE_MERGE_EXON_COUNTS", unset = "")

sample_files <- parse_path_list(merge_samples_env, project_dir, cfg)
exon_files   <- parse_path_list(merge_exon_env, project_dir, cfg)

if (length(sample_files) == 0) {
  sample_files <- file.path(project_dir, "samples.tsv")
}
if (length(exon_files) == 0) {
  exon_files <- file.path(project_dir, "counts", "exon_count.txt")
}

samples_tsv <- paste(sample_files, collapse = ", ")
exon_count  <- paste(exon_files, collapse = ", ")

# --- Load packages ---
suppressPackageStartupMessages({
  library(DEXSeq)
  library(tidyverse)
  library(data.table)
  library(magrittr)
  library(ggrepel)
  library(BiocParallel)
})

# --- Load plot config ---
script_dir <- tryCatch({
  normalizePath(dirname(sys.frame(1)$ofile), mustWork = TRUE)
}, error = function(e) {
  candidates <- c(
    file.path(project_dir, "scripts"),
    file.path(cfg$PIPELINE_DIR %||% ".", "scripts")
  )
  for (d in candidates) {
    if (file.exists(file.path(d, "plot_utils.R"))) return(normalizePath(d))
  }
  "scripts"
})
source(file.path(script_dir, "plot_utils.R"))
pcfg <- load_plot_config()

set.seed(123)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================\n")
cat(" Step 6: Exon-level DE (DEXSeq)\n")
cat(" Exon counts: ", exon_count, "\n")
cat(" Output     : ", output_dir, "\n")
cat("============================================\n")

# --- Read sample info ---
sample_info <- read_sample_tables(sample_files)

# --- Read exon count matrix ---
# Format from featureCounts -f: Geneid, Chr, Start, End, Strand, Length, counts...
exon_raw <- read_exon_count_tables(exon_files)

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
colnames(count_mat) <- clean_names

# Match to samples
count_keys <- normalize_sample_name(clean_names)
sample_keys <- normalize_sample_name(sample_info$sample_name)
matched_idx <- match(sample_keys, count_keys)

if (all(is.na(matched_idx))) {
  stop("ERROR: No sample names match exon count columns.")
}
sample_info <- sample_info[!is.na(matched_idx), ]
matched_idx <- matched_idx[!is.na(matched_idx)]

count_mat <- count_mat[, matched_idx, drop = FALSE]
colnames(count_mat) <- sample_info$sample_name

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
# Use platform-appropriate parallelism
if (.Platform$OS.type == "windows") {
  BPPARAM <- SerialParam()
} else {
  BPPARAM <- MulticoreParam(workers = min(4, parallel::detectCores()))
}
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

# --- Visualization: top genes (using plot_config) ---
top_genes <- head(gene_summary$gene_name, 20)
if (length(top_genes) > 0) {
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, showWarnings = FALSE)

  exon_size <- pcfg$figure_size$exon_usage %||% list(width = 12, height = 6)
  fmt <- pcfg$output$format %||% "png"
  dpi <- pcfg$output$dpi %||% 600

  for (g in top_genes) {
    tryCatch({
      out_file <- file.path(plot_dir, paste0(g, "_exon_usage.", fmt))
      if (fmt == "png") {
        png(out_file,
            width = exon_size$width, height = exon_size$height,
            units = "in", res = dpi)
      } else if (fmt == "pdf") {
        pdf(out_file,
            width = exon_size$width, height = exon_size$height)
      } else if (fmt == "svg") {
        svg(out_file,
            width = exon_size$width, height = exon_size$height)
      } else {
        png(out_file,
            width = exon_size$width, height = exon_size$height,
            units = "in", res = dpi)
      }
      plotDEXSeq(dxr, g, legend = TRUE, cex.axis = 1.2, cex = 1.3,
                 lwd = 2, displayTranscripts = FALSE)
      dev.off()
      message("[plot_utils] 保存: ", out_file)
    }, error = function(e) {
      message("  Could not plot: ", g, " - ", e$message)
      tryCatch(dev.off(), error = function(e2) NULL)
    })
  }
}

cat("============================================\n")
cat(" Step 6: Exon-level DE Complete\n")
cat("============================================\n")
