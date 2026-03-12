#!/usr/bin/env Rscript
# =============================================================================
# Step 6: Single-sample GSEA and selected-gene visualization
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config.sh")

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

  dat <- data.table::fread(
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

clean_key <- function(x) {
  toupper(trimws(x))
}

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "yes", "y")
}

parse_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    return(character())
  }
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals[nzchar(vals)]
}

scale_matrix <- function(mat, scale_mode) {
  if (scale_mode == "row") {
    scaled <- t(scale(t(mat)))
  } else if (scale_mode == "column") {
    scaled <- scale(mat)
  } else {
    scaled <- mat
  }
  scaled[!is.finite(scaled)] <- 0
  scaled
}

read_gmt_file <- function(path) {
  gmt_lines <- readLines(path, warn = FALSE)
  gmt_lines <- gmt_lines[nzchar(trimws(gmt_lines))]
  gmt_split <- strsplit(gmt_lines, "\t", fixed = TRUE)
  out <- list()
  for (entries in gmt_split) {
    if (length(entries) < 3) {
      next
    }
    out[[entries[1]]] <- unique(clean_key(entries[-c(1, 2)]))
  }
  out
}

resolve_gmt_paths <- function(gmt_files, project_dir, cfg = list()) {
  if (length(gmt_files) == 0) {
    return(character())
  }

  unique(vapply(gmt_files, function(item) {
    trimmed <- trimws(item)
    if (!nzchar(trimmed)) {
      return("")
    }
    if (grepl("^/", trimmed) || grepl("^~", trimmed) || grepl("\\$", trimmed) || grepl("/", trimmed)) {
      return(resolve_input_path(trimmed, project_dir, cfg))
    }
    file.path(project_dir, "Gene_set", trimmed)
  }, character(1)))
}

resolve_gene_sets <- function(target_gene_sets, msig_categories, species, gmt_paths = character()) {
  collections <- list()

  for (cat in msig_categories) {
    msig <- msigdbr::msigdbr(species = species, category = cat)
    if (nrow(msig) == 0) {
      next
    }
    msig$gs_name_clean <- clean_key(msig$gs_name)
    msig$gene_symbol_clean <- clean_key(msig$gene_symbol)
    split_genes <- split(msig$gene_symbol_clean, msig$gs_name_clean)
    split_genes <- lapply(split_genes, function(v) unique(v[nzchar(v)]))
    collections <- c(collections, split_genes[setdiff(names(split_genes), names(collections))])
  }

  for (gmt_path in unique(gmt_paths[nzchar(gmt_paths)])) {
    if (!file.exists(gmt_path)) {
      warning("ssGSEA gene set file not found: ", gmt_path)
      next
    }
    custom_sets <- read_gmt_file(gmt_path)
    custom_keys <- clean_key(names(custom_sets))
    names(custom_sets) <- custom_keys
    collections[names(custom_sets)] <- custom_sets
  }

  target_keys <- clean_key(target_gene_sets)
  matched <- collections[target_keys]
  names(matched) <- target_keys
  missing <- target_gene_sets[!(target_keys %in% names(collections))]

  list(
    gene_sets = matched[!vapply(matched, is.null, logical(1))],
    missing = missing
  )
}

build_annotation <- function(sample_info, sample_names) {
  ann_col <- data.frame(Group = factor(sample_info$group[match(sample_names, sample_info$sample_name)]))
  rownames(ann_col) <- sample_names
  ann_col
}

save_heatmap <- function(mat, ann_col, main, filepath, cfg, show_rownames = NULL, show_colnames = NULL) {
  cols <- colorRampPalette(c(
    cfg$colors$heatmap$low %||% "#313695",
    cfg$colors$heatmap$mid %||% "#FFFFBF",
    cfg$colors$heatmap$high %||% "#A50026"
  ))(100)

  width <- cfg$figure_size$heatmap$width %||% 10
  base_height <- cfg$figure_size$heatmap$height %||% 8
  height <- max(base_height, nrow(mat) * 0.22 + 2)

  grDevices::pdf(filepath, width = width, height = height)
  pheatmap::pheatmap(
    mat,
    color = cols,
    annotation_col = ann_col,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    cluster_rows = cfg$heatmap$cluster_rows %||% TRUE,
    cluster_cols = cfg$heatmap$cluster_cols %||% TRUE,
    show_rownames = show_rownames %||% (cfg$heatmap$show_rownames %||% FALSE),
    show_colnames = show_colnames %||% (cfg$heatmap$show_colnames %||% TRUE),
    fontsize_row = 8,
    fontsize_col = 8,
    main = main
  )
  grDevices::dev.off()
}

build_gene_plot <- function(data_gene, gene_name, group_colors, paired = TRUE, expression_label = "log1p(normalized counts)") {
  p <- ggplot2::ggplot(
    data_gene,
    ggplot2::aes(x = .data$group, y = .data$expr, fill = .data$group)
  ) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.35) +
    ggplot2::geom_point(shape = 21, size = 2.8, color = "black",
                        position = ggplot2::position_jitter(width = 0.08, height = 0)) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::labs(title = gene_name, x = NULL, y = expression_label) +
    ggplot2::theme(legend.position = "none")

  summary_tbl <- data.frame(
    Gene = gene_name,
    mode = if (paired) "paired" else "grouped",
    n = nrow(data_gene),
    p_value = NA_real_,
    statistic = NA_real_,
    stringsAsFactors = FALSE
  )

  if (paired) {
    p <- p + ggplot2::geom_line(ggplot2::aes(group = .data$patient), color = "gray60", alpha = 0.6)
    wide <- tidyr::pivot_wider(data_gene[, c("patient", "group", "expr")],
                               names_from = "group", values_from = "expr")
    group_levels <- unique(as.character(data_gene$group))
    if (length(group_levels) == 2 && all(group_levels %in% colnames(wide))) {
      t_res <- tryCatch(stats::t.test(wide[[group_levels[1]]], wide[[group_levels[2]]], paired = TRUE),
                        error = function(e) NULL)
      if (!is.null(t_res)) {
        summary_tbl$p_value <- t_res$p.value
        summary_tbl$statistic <- unname(t_res$statistic)
        p <- p + ggplot2::labs(subtitle = paste0("paired t-test p=", signif(t_res$p.value, 3)))
      }
    }
  }

  list(plot = p, summary = summary_tbl)
}

cfg <- parse_config(config_file)

project_dir <- getwd()
output_dir <- resolve_input_path(cfg$SSGSEA_OUTPUT_DIR %||% file.path(project_dir, "single_sample"), project_dir, cfg)
deseq2_dir <- file.path(project_dir, "DESeq2")
counts_file <- file.path(deseq2_dir, "Counts.normalized.txt")
custom_gmt <- resolve_input_path(cfg$CUSTOM_GMT %||% "", project_dir, cfg)

merge_samples_env <- Sys.getenv("PIPELINE_MERGE_SAMPLES", unset = "")
sample_files <- parse_path_list(merge_samples_env, project_dir, cfg)
if (length(sample_files) == 0) {
  sample_files <- file.path(project_dir, "samples.tsv")
}

ssgsea_categories <- parse_csv(cfg$SSGSEA_MSIG_CATEGORIES %||% "H,C2,C3,C4,C5,C6,C7")
ssgsea_gmt_files <- parse_csv(cfg$SSGSEA_GMT_FILES %||% "")
ssgsea_min_size <- as.integer(cfg$SSGSEA_MIN_SIZE %||% 10)
ssgsea_max_size <- as.integer(cfg$SSGSEA_MAX_SIZE %||% 5000)
ssgsea_min_expr <- as.numeric(cfg$SSGSEA_MIN_EXPR %||% 1)
ssgsea_min_fraction <- as.numeric(cfg$SSGSEA_MIN_SAMPLES_FRACTION %||% 0.5)
gene_plot_groups <- parse_csv(cfg$GENE_PLOT_GROUPS %||% "")
gene_plot_paired <- to_bool(cfg$GENE_PLOT_PAIRED %||% "true", default = TRUE)
genome <- cfg$GENOME %||% "hg38"
species <- ifelse(genome == "hg38", "Homo sapiens", "Mus musculus")
ssgsea_gmt_paths <- resolve_gmt_paths(ssgsea_gmt_files, project_dir, cfg)
if (nzchar(custom_gmt)) {
  ssgsea_gmt_paths <- unique(c(ssgsea_gmt_paths, custom_gmt))
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(GSVA)
  library(msigdbr)
  library(pheatmap)
})

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

cat("============================================\n")
cat(" Step 6: ssGSEA and selected-gene visualization\n")
cat(" Counts    : ", counts_file, "\n")
cat(" Samples   : ", paste(sample_files, collapse = ", "), "\n")
cat(" Output    : ", output_dir, "\n")
cat("============================================\n")

if (!file.exists(counts_file)) {
  stop("ERROR: Normalized counts not found. Run Step 4 first: ", counts_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "per_gene_plots"), showWarnings = FALSE, recursive = TRUE)

sample_info <- read_sample_tables(sample_files)
expr_df <- data.table::fread(counts_file) %>% as.data.frame()
colnames(expr_df)[1] <- "Gene"
rownames(expr_df) <- make.unique(expr_df$Gene)
expr_df$Gene <- NULL
expr_mat <- as.matrix(expr_df)
storage.mode(expr_mat) <- "double"

sample_match <- match(sample_info$sample_name, colnames(expr_mat))
if (all(is.na(sample_match))) {
  stop("ERROR: No sample names in samples.tsv match Counts.normalized.txt columns.")
}
sample_info <- sample_info[!is.na(sample_match), ]
expr_mat <- expr_mat[, sample_match[!is.na(sample_match)], drop = FALSE]
colnames(expr_mat) <- sample_info$sample_name

min_samples <- max(1L, ceiling(ncol(expr_mat) * ssgsea_min_fraction))
keep <- rowSums(expr_mat >= ssgsea_min_expr) >= min_samples
expr_keep <- expr_mat[keep, , drop = FALSE]
expr_log <- log1p(expr_keep)

expr_log_upper <- expr_log
rownames(expr_log_upper) <- clean_key(rownames(expr_log_upper))
if (any(duplicated(rownames(expr_log_upper)))) {
  expr_log_upper <- rowsum(expr_log_upper, group = rownames(expr_log_upper), reorder = FALSE)
}

selected_gene_sets <- pcfg$ssgsea$gene_sets %||% character()
selected_genes <- pcfg$selected_genes$genes %||% character()

if (length(selected_gene_sets) > 0) {
  gene_set_info <- resolve_gene_sets(selected_gene_sets, ssgsea_categories, species, ssgsea_gmt_paths)
  if (length(gene_set_info$missing) > 0) {
    warning("The following ssGSEA gene sets were not found: ", paste(gene_set_info$missing, collapse = ", "))
  }

  if (length(gene_set_info$gene_sets) > 0) {
    gsva_param <- gsvaParam(
      exprData = expr_log_upper,
      geneSets = gene_set_info$gene_sets,
      kcdf = "Gaussian",
      minSize = ssgsea_min_size,
      maxSize = ssgsea_max_size,
      tau = 1,
      maxDiff = TRUE,
      absRanking = FALSE
    )
    ssgsea_scores <- gsva(gsva_param)
    ssGSEA_z <- scale_matrix(ssgsea_scores, "row")
    ann_col <- build_annotation(sample_info, colnames(ssGSEA_z))

    data.table::fwrite(
      as.data.frame(ssgsea_scores) %>% tibble::rownames_to_column("Pathway"),
      file.path(output_dir, "ssgsea_scores.tsv"),
      sep = "\t"
    )
    data.table::fwrite(
      as.data.frame(ssGSEA_z) %>% tibble::rownames_to_column("Pathway"),
      file.path(output_dir, "ssgsea_scores_rowZ.tsv"),
      sep = "\t"
    )
    save_heatmap(
      ssGSEA_z,
      ann_col,
      main = "Single-sample GSEA",
      filepath = file.path(output_dir, "ssgsea_heatmap.pdf"),
      cfg = pcfg,
      show_rownames = TRUE,
      show_colnames = pcfg$selected_genes$show_colnames %||% TRUE
    )
  } else {
    warning("No valid gene sets resolved for ssGSEA. Skipping ssGSEA output.")
  }
}

if (length(selected_genes) > 0) {
  gene_map <- setNames(rownames(expr_mat), clean_key(rownames(expr_mat)))
  gene_keys <- clean_key(selected_genes)
  present_genes <- unique(unname(gene_map[gene_keys[gene_keys %in% names(gene_map)]]))
  missing_genes <- selected_genes[!(gene_keys %in% names(gene_map))]

  if (length(missing_genes) > 0) {
    warning("Selected genes not found: ", paste(missing_genes, collapse = ", "))
  }

  if (length(present_genes) > 0) {
    gene_mat <- expr_mat[present_genes, , drop = FALSE]
    gene_mat_log <- log1p(gene_mat)
    heatmap_scale <- pcfg$heatmap$scale %||% "row"
    gene_mat_scaled <- scale_matrix(gene_mat_log, heatmap_scale)
    ann_col <- build_annotation(sample_info, colnames(gene_mat_scaled))

    data.table::fwrite(
      as.data.frame(gene_mat) %>% tibble::rownames_to_column("Gene"),
      file.path(output_dir, "selected_genes_normalized_counts.tsv"),
      sep = "\t"
    )
    data.table::fwrite(
      as.data.frame(gene_mat_log) %>% tibble::rownames_to_column("Gene"),
      file.path(output_dir, "selected_genes_log1p.tsv"),
      sep = "\t"
    )
    data.table::fwrite(
      as.data.frame(gene_mat_scaled) %>% tibble::rownames_to_column("Gene"),
      file.path(output_dir, "selected_genes_scaled.tsv"),
      sep = "\t"
    )
    save_heatmap(
      gene_mat_scaled,
      ann_col,
      main = pcfg$selected_genes$heatmap_title %||% "Selected genes",
      filepath = file.path(output_dir, "selected_gene_heatmap.pdf"),
      cfg = pcfg,
      show_rownames = pcfg$selected_genes$show_rownames %||% TRUE,
      show_colnames = pcfg$selected_genes$show_colnames %||% TRUE
    )

    plot_sample_info <- sample_info
    if (length(gene_plot_groups) == 0) {
      unique_groups <- unique(as.character(plot_sample_info$group))
      if (length(unique_groups) == 2) {
        gene_plot_groups <- unique_groups
      }
    }

    plot_sample_info <- plot_sample_info[plot_sample_info$group %in% gene_plot_groups, , drop = FALSE]
    group_colors <- get_group_colors(gene_plot_groups, pcfg)
    summary_list <- list()

    if (length(gene_plot_groups) == 2 && nrow(plot_sample_info) > 0) {
      patient_ids <- gsub(
        paste0("_", gene_plot_groups[1], "$|_", gene_plot_groups[2], "$"),
        "",
        as.character(plot_sample_info$sample_name)
      )
      plot_sample_info$patient <- patient_ids
      plot_sample_info$group <- factor(plot_sample_info$group, levels = gene_plot_groups)

      if (gene_plot_paired) {
        common_patients <- intersect(
          plot_sample_info$patient[plot_sample_info$group == gene_plot_groups[1]],
          plot_sample_info$patient[plot_sample_info$group == gene_plot_groups[2]]
        )
        plot_sample_info <- plot_sample_info[plot_sample_info$patient %in% common_patients, , drop = FALSE]
      }

      if (nrow(plot_sample_info) > 0) {
        for (gene_name in present_genes) {
          data_gene <- data.frame(
            sample = plot_sample_info$sample_name,
            group = factor(plot_sample_info$group, levels = gene_plot_groups),
            patient = plot_sample_info$patient,
            expr = as.numeric(gene_mat_log[gene_name, plot_sample_info$sample_name]),
            stringsAsFactors = FALSE
          )
          plot_result <- build_gene_plot(
            data_gene,
            gene_name,
            group_colors,
            paired = gene_plot_paired,
            expression_label = pcfg$selected_genes$expression_label %||% "log1p(normalized counts)"
          )
          plot_result$plot <- plot_result$plot + theme_pipeline(pcfg)
          save_plot(
            plot_result$plot,
            paste0("per_gene_plots/", gsub("[^A-Za-z0-9._-]+", "_", gene_name), "_paired_lineplot"),
            type = "gene_plot",
            outdir = output_dir,
            cfg = pcfg
          )
          summary_list[[gene_name]] <- plot_result$summary
        }
      } else {
        warning("No samples available for selected per-gene plots after group/pair filtering.")
      }
    } else {
      warning("GENE_PLOT_GROUPS could not be resolved to exactly 2 groups. Skipping per-gene plots.")
    }

    if (length(summary_list) > 0) {
      data.table::fwrite(
        dplyr::bind_rows(summary_list),
        file.path(output_dir, "per_gene_plots", "gene_plot_summary.tsv"),
        sep = "\t"
      )
    }
  } else {
    warning("No selected genes matched the normalized count matrix.")
  }
}

cat("============================================\n")
cat(" Step 6: Complete\n")
cat("============================================\n")