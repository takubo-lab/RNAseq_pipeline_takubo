#!/usr/bin/env Rscript
# =============================================================================
# Step 4: DESeq2 Differential Expression Analysis
# =============================================================================
# Performs all-group pairwise comparisons using DESeq2.
# Reads configuration from config.sh and sample info from samples.tsv.
#
# Usage:
#   Rscript scripts/04_deseq2.R [config.sh]
# =============================================================================

# --- Parse config ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config.sh")

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

# Resolve project directory
project_dir <- getwd()
samples_tsv <- file.path(project_dir, "samples.tsv")
count_file  <- file.path(project_dir, "counts", "count.txt")
output_dir  <- file.path(project_dir, "DESeq2")

min_count   <- as.integer(cfg$MIN_COUNT_FILTER %||% 50)
pca_top     <- as.integer(cfg$PCA_TOP_GENES %||% 500)
padj_thr    <- as.numeric(cfg$VOLCANO_PADJ %||% 0.001)
lfc_thr     <- as.numeric(cfg$VOLCANO_LOG2FC %||% 1)
gene_list   <- unlist(strsplit(cfg$HIGHLIGHT_GENES %||% "", ","))

# --- Load packages ---
suppressPackageStartupMessages({
  library(DESeq2)
  library(umap)
  library(tidyverse)
  library(data.table)
  library(magrittr)
  library(ggrepel)
})

set.seed(123)

cat("============================================\n")
cat(" Step 4: DESeq2 Analysis\n")
cat(" Count file : ", count_file, "\n")
cat(" Samples    : ", samples_tsv, "\n")
cat(" Output     : ", output_dir, "\n")
cat("============================================\n")

# --- Read sample info ---
sample_info <- fread(samples_tsv, sep = "\t", header = TRUE)
sample_info <- sample_info[!grepl("^#", sample_info$sample_name), ]
colnames(sample_info) <- c("sample_name", "group", "fq_prefix", "lane_suffix")

# --- Read count matrix ---
raw <- fread(count_file, header = TRUE)
geneid_col <- colnames(raw)[1]
rawCounts  <- raw %>%
  tibble::column_to_rownames(geneid_col)

# Clean column names (remove BAM path prefixes if present)
clean_names <- colnames(rawCounts)
clean_names <- gsub(".*[/\\\\]", "", clean_names)
clean_names <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", clean_names)
clean_names <- gsub("-", "_", clean_names)
colnames(rawCounts) <- clean_names

# Match samples to columns
matched <- sample_info$sample_name[sample_info$sample_name %in% clean_names]
if (length(matched) == 0) {
  stop("ERROR: No sample names in samples.tsv match count matrix columns.\n",
       "  Count columns: ", paste(clean_names[1:min(5, length(clean_names))], collapse = ", "), "\n",
       "  Sample names : ", paste(sample_info$sample_name[1:min(5, nrow(sample_info))], collapse = ", "))
}

rawCounts <- rawCounts[, matched, drop = FALSE]
sample_info <- sample_info[match(matched, sample_info$sample_name), ]

samples <- factor(sample_info$sample_name, levels = sample_info$sample_name)
groups  <- factor(sample_info$group)

sampleData <- data.frame(
  samples = samples,
  class   = groups,
  row.names = as.character(samples)
)

# --- Filter low-count genes ---
rawCounts_filt <- rawCounts[rowSums(rawCounts) > min_count, ]
rawCounts_filt <- rawCounts_filt[!duplicated(rawCounts_filt), ]

cat("Genes after filtering: ", nrow(rawCounts_filt), "\n")

# --- Create DESeq2 object ---
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dds <- DESeqDataSetFromMatrix(
  countData = rawCounts_filt,
  colData   = sampleData,
  design    = ~class
)

# Normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
norm_df <- data.frame(normalized_counts) %>%
  rownames_to_column("GeneName") %>%
  tibble()
fwrite(norm_df, file.path(output_dir, "Counts.normalized.txt"), sep = "\t")

# --- Run DESeq2 ---
dds <- DESeq(dds)

# --- PCA ---
vst_out  <- vst(dds)
stdev    <- apply(assay(vst_out), 1, sd)
top_n    <- min(pca_top, nrow(assay(vst_out)))
vst_short <- assay(vst_out)[rev(order(stdev))[1:top_n], ]

pcaData    <- data.frame(prcomp(t(vst_short))$x)
percentVar <- round(100 * summary(prcomp(t(vst_short)))$importance[2, ])

p_pca <- ggplot(pcaData, aes(x = PC1, y = PC2, fill = sampleData$class)) +
  geom_point(size = 8, alpha = 0.8, pch = 21, color = "Black") +
  geom_text_repel(aes(label = rownames(pcaData)), max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  labs(fill = "Group") +
  ggtitle("PCA") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

ggsave(file.path(output_dir, "PCA_class.png"), p_pca, width = 8, height = 7, dpi = 600)

# --- UMAP ---
config_umap <- umap.defaults
config_umap$n_neighbors <- min(3, ncol(vst_short) - 1)
u <- umap(t(vst_short), config = config_umap)
u_df <- data.frame(u$layout)

p_umap <- ggplot(u_df, aes(x = X1, y = X2, fill = sampleData$class)) +
  geom_point(size = 5, alpha = 0.8, pch = 21, color = "Black") +
  geom_text_repel(aes(label = rownames(u_df)), max.overlaps = 20) +
  xlab("UMAP_1") + ylab("UMAP_2") +
  labs(fill = "Group") +
  ggtitle("UMAP") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

ggsave(file.path(output_dir, "UMAP_all.png"), p_umap, width = 9, height = 7, dpi = 600)

# --- PC Loadings ---
loadings <- prcomp(t(vst_short))$rotation
loadings_df <- data.frame(PC1 = loadings[, 1], PC2 = loadings[, 2]) %>%
  rownames_to_column("Gene")
fwrite(loadings_df, file.path(output_dir, "PC_loadings.tsv"), sep = "\t")

# --- Pairwise comparisons (all groups) ---
all_levels <- levels(groups)
all_pairs  <- combn(all_levels, 2, simplify = FALSE)

cat("Performing ", length(all_pairs), " pairwise comparisons...\n")

# Volcano plot function
plot_volcano <- function(deseq2Results, gene_list = c(""), file_out,
                         padj_thr = 0.001, lfc_thr = 1) {
  set.seed(42)
  df <- data.frame(deseq2Results) %>%
    mutate(sig = case_when(
      geneID %in% gene_list ~ "marked",
      padj < padj_thr & log2FoldChange >  lfc_thr & !(geneID %in% gene_list) ~ "up",
      padj < padj_thr & log2FoldChange < -lfc_thr & !(geneID %in% gene_list) ~ "down",
      TRUE ~ "non"
    )) %>%
    mutate(mark = geneID %in% gene_list) %>%
    arrange(mark)

  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), fill = sig)) +
    geom_point(pch = 21, aes(size = mark), alpha = 0.7) +
    scale_fill_manual(values = c(down = "#33DFD4", marked = "#FFFF00",
                                 non = "#CACACA", up = "#F73719")) +
    scale_size_manual(values = c(1.5, 3)) +
    geom_text_repel(data = filter(df, mark & log2FoldChange > 0),
                    aes(label = geneID), size = 4, max.overlaps = 10,
                    nudge_y = 20, nudge_x = 30, box.padding = 0.5, color = "#000066") +
    geom_text_repel(data = filter(df, mark & log2FoldChange < 0),
                    aes(label = geneID), size = 4, max.overlaps = 10,
                    nudge_y = 20, nudge_x = -30, box.padding = 0.5, color = "#000066") +
    geom_vline(xintercept = lfc_thr, linetype = "dashed") +
    geom_vline(xintercept = -lfc_thr, linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    theme_bw()
  ggsave(paste0(file_out, "_volcano.png"), width = 6, height = 6)
}

# Paired test function (for paired experimental designs)
my_paired_test <- function(norm_counts, sampleData, group1, group2) {
  sampleData <- as.data.frame(sampleData)
  colnames(sampleData) <- c("samples", "class")
  sampleData$patient <- gsub(paste0("_", group1, "$|_", group2, "$"), "",
                             as.character(sampleData$samples))

  df_g1 <- sampleData %>% filter(class == group1) %>% arrange(patient)
  df_g2 <- sampleData %>% filter(class == group2) %>% arrange(patient)
  common_patients <- intersect(df_g1$patient, df_g2$patient)

  if (length(common_patients) < 2) {
    message("  Not enough paired samples, skipping paired test.")
    return(NULL)
  }

  df_g1 <- df_g1 %>% filter(patient %in% common_patients) %>% arrange(patient)
  df_g2 <- df_g2 %>% filter(patient %in% common_patients) %>% arrange(patient)
  message(sprintf("  Paired test: %d patients (%s vs %s)",
                  length(common_patients), group1, group2))

  mat_g1 <- log2(norm_counts[, as.character(df_g1$samples), drop = FALSE] + 1)
  mat_g2 <- log2(norm_counts[, as.character(df_g2$samples), drop = FALSE] + 1)

  n_genes     <- nrow(norm_counts)
  ttest_p     <- numeric(n_genes)
  ttest_stat  <- numeric(n_genes)
  wilcox_p    <- numeric(n_genes)
  wilcox_stat <- numeric(n_genes)
  mean_log2fc <- numeric(n_genes)

  for (i in seq_len(n_genes)) {
    x_g1 <- as.numeric(mat_g1[i, ])
    x_g2 <- as.numeric(mat_g2[i, ])
    mean_log2fc[i] <- mean(x_g1 - x_g2, na.rm = TRUE)
    t_res <- tryCatch(t.test(x_g1, x_g2, paired = TRUE),
                      error = function(e) list(statistic = NA_real_, p.value = NA_real_))
    ttest_p[i]    <- t_res$p.value
    ttest_stat[i] <- as.numeric(t_res$statistic)
    w_res <- tryCatch(wilcox.test(x_g1, x_g2, paired = TRUE, exact = FALSE),
                      error = function(e) list(statistic = NA_real_, p.value = NA_real_))
    wilcox_p[i]    <- w_res$p.value
    wilcox_stat[i] <- as.numeric(w_res$statistic)
  }

  data.frame(
    geneID            = rownames(norm_counts),
    paired_mean_log2FC = mean_log2fc,
    ttest_stat        = ttest_stat,
    ttest_p           = ttest_p,
    ttest_padj_BH     = p.adjust(ttest_p, method = "BH"),
    wilcox_stat       = wilcox_stat,
    wilcox_p          = wilcox_p,
    wilcox_padj_BH    = p.adjust(wilcox_p, method = "BH"),
    stringsAsFactors  = FALSE
  )
}

# Run all pairwise comparisons
for (pair in all_pairs) {
  g1 <- pair[1]
  g2 <- pair[2]
  comp_name <- paste0(g1, "_vs_", g2)
  message(sprintf("[DESeq2] %s", comp_name))

  res <- results(dds, contrast = c("class", g1, g2))
  res_df <- as_tibble(res, rownames = "geneID") %>% arrange(padj)
  fwrite(res_df, file.path(output_dir, paste0(comp_name, "_DE_Results.tsv")), sep = "\t")

  # Paired test (if applicable)
  norm_sub <- counts(dds, normalized = TRUE)
  sd_sub   <- sampleData[sampleData$class %in% c(g1, g2), , drop = FALSE]
  norm_sub <- norm_sub[, as.character(sd_sub$samples), drop = FALSE]

  paired_res <- my_paired_test(norm_sub, sd_sub, group1 = g1, group2 = g2)
  if (!is.null(paired_res)) {
    merged <- left_join(as.data.frame(res_df), paired_res, by = "geneID") %>%
      arrange(padj)
    fwrite(merged, file.path(output_dir, paste0(comp_name, "_DE_Results_with_paired.tsv")),
           sep = "\t")
  }

  # Volcano
  plot_volcano(res_df, gene_list,
               file.path(output_dir, comp_name),
               padj_thr = padj_thr, lfc_thr = lfc_thr)
}

cat("============================================\n")
cat(" Step 4: DESeq2 Analysis Complete\n")
cat("============================================\n")
