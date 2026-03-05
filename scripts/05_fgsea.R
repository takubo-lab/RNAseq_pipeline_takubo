#!/usr/bin/env Rscript
# =============================================================================
# Step 5: fGSEA (fast Gene Set Enrichment Analysis)
# =============================================================================
# Runs fGSEA on DESeq2 results for MSigDB gene sets + custom gene sets.
#
# Usage:
#   Rscript scripts/05_fgsea.R [config.sh]
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

project_dir  <- getwd()
deseq2_dir   <- file.path(project_dir, "DESeq2")
fgsea_dir    <- file.path(project_dir, "fGSEA")
genome       <- cfg$GENOME %||% "hg38"

fgsea_min    <- as.integer(cfg$FGSEA_MIN_SIZE %||% 15)
fgsea_max    <- as.integer(cfg$FGSEA_MAX_SIZE %||% 500)
fgsea_nperm  <- as.integer(cfg$FGSEA_NPERM %||% 10000)
fgsea_nproc  <- as.integer(cfg$FGSEA_NPROC %||% 12)
fgsea_padj   <- as.numeric(cfg$FGSEA_PADJ %||% 0.05)
custom_gmt   <- cfg$CUSTOM_GMT %||% ""

species <- ifelse(genome == "hg38", "Homo sapiens", "Mus musculus")

# --- Load packages ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(magrittr)
  library(fgsea)
  library(msigdbr)
})

set.seed(123)
dir.create(fgsea_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================\n")
cat(" Step 5: fGSEA\n")
cat(" Species : ", species, "\n")
cat(" DESeq2  : ", deseq2_dir, "\n")
cat(" Output  : ", fgsea_dir, "\n")
cat("============================================\n")

# --- Build gene set collection ---
msig_categories <- c("H", "C2", "C3", "C4", "C5", "C6", "C7")
msigdb_list <- lapply(msig_categories, function(cat) {
  msigdb <- msigdbr::msigdbr(species = species, category = cat)
  split(msigdb$gene_symbol, msigdb$gs_name)
})
names(msigdb_list) <- msig_categories

# Add custom gene set if available
if (nchar(custom_gmt) > 0 && file.exists(custom_gmt)) {
  cat("Loading custom gene set: ", custom_gmt, "\n")
  gmt_raw <- fread(custom_gmt, sep = "\t", header = FALSE, fill = TRUE, quote = "")
  gmt_t   <- t(gmt_raw) %>% data.frame()
  colnames(gmt_t) <- gmt_t[1, ]
  gmt_t <- gmt_t[-c(1, 2), , drop = FALSE]
  gmt_list <- as.list(gmt_t)
  gmt_list <- lapply(gmt_list, function(x) {
    x <- toupper(x)
    x[x != ""]
  })
  custom_name <- tools::file_path_sans_ext(basename(custom_gmt))
  msigdb_list[[custom_name]] <- gmt_list
  cat("  Added ", length(gmt_list), " custom sets from ", custom_name, "\n")
}

# --- fGSEA function ---
run_fgsea <- function(de_file, gene_set_list) {
  de <- fread(de_file, sep = "\t")
  ranks <- de$log2FoldChange
  names(ranks) <- de$geneID
  ranks <- na.omit(ranks)

  # Get comparison name
  file_name <- basename(de_file)
  file_name <- gsub("_DE_Results\\.tsv$", "", file_name)
  comp_dir  <- file.path(fgsea_dir, paste0("fGSEA_", file_name))
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

  cat_names <- names(gene_set_list)

  for (i in seq_along(gene_set_list)) {
    cat("  ", file_name, " - ", cat_names[i], "\n")
    fgsea_res <- fgseaMultilevel(
      pathways      = gene_set_list[[i]],
      stats         = ranks,
      eps           = 0,
      minSize       = fgsea_min,
      maxSize       = fgsea_max,
      nPermSimple   = fgsea_nperm,
      nproc         = fgsea_nproc
    )
    fgsea_res %<>% mutate(direction = case_when(
      padj < fgsea_padj & NES > 0 ~ "upReg",
      padj < fgsea_padj & NES < 0 ~ "downReg",
      TRUE ~ "none"
    ))

    out_file <- file.path(comp_dir,
                          paste0("fgsea_Results_", file_name, "_", cat_names[i], ".tsv"))
    fwrite(fgsea_res, out_file, sep = "\t")
  }
}

# --- Run on all DE result files ---
de_files <- list.files(deseq2_dir, pattern = "_DE_Results\\.tsv$", full.names = TRUE)
# Exclude "_with_paired" files
de_files <- de_files[!grepl("_with_paired", de_files)]

if (length(de_files) == 0) {
  stop("ERROR: No DE result files found in ", deseq2_dir)
}

cat("Found ", length(de_files), " DE result files.\n")

for (de_file in de_files) {
  cat("Processing: ", basename(de_file), "\n")
  run_fgsea(de_file, msigdb_list)
}

cat("============================================\n")
cat(" Step 5: fGSEA Complete\n")
cat("============================================\n")
