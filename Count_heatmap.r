# 必要パッケージ
# install.packages(c("readr","dplyr","stringr","pheatmap"))
# BiocManager::install("ComplexHeatmap")
library(readr)
library(dplyr)
library(stringr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
setwd("C://Users/hikob/Dropbox/Research/PAI-I_inhibitor_project/250821_RNAseq/")
# ===== 入力 =====
# 正規化済みカウント表（1列目が遺伝子名、以降がサンプル）
count_path <- "Counts.normalized.txt"

# 解析したい遺伝子（例：適宜置き換えてください）
genes_of_interest <- c(
  "DDIT3","CFB","ARRDC3","FOS","ATF3","H1-3","H2BC3","SCIMP","BNIPL","LILRA4","CD69",
  "SDE2","IFI30","PAIP1P1","H3C13","SOCS3","H1-4","H3-3B")

genes_of_interest <- c(
  "DDIT3","CFB","ARRDC3","FOS","ATF3","H1-3","H2BC3","SCIMP","BNIPL","LILRA4","CD69",
  "SDE2","IFI30","PAIP1P1","H3C13","SOCS3","H1-4","H3-3B","CDKN1B","H4C6","EXO1","BRCA1",
  "VCAN","SMAD7","HMOX1","MMP9")



# オプション：行スケーリング(Zスコア)の有無
scale_by_gene <- TRUE

# ===== 読み込み & 整形 =====
dat <- read_tsv(count_path, show_col_types = FALSE)
# 1列目が遺伝子名のはず：列名を明示
colnames(dat)[1] <- "Gene"
# 重複遺伝子名があってもユニーク化
dat$Gene <- make.unique(dat$Gene)

# 行名に遺伝子を設定して数値行列へ
mat <- dat %>% as.data.frame()
rownames(mat) <- mat$Gene
mat$Gene <- NULL
mat <- as.matrix(mat)

# サンプル注釈（列名から群情報を抽出：例としてPBで始まる記号を「PB」に）
# 列名が「PB1, PB2, ...」のような場合の例。任意のパターンに合わせて調整ください。
sample_groups <- c(rep("Pre",3),rep("Post",3))
ann_col <- data.frame(Group = factor(sample_groups))
row.names(ann_col) <- colnames(mat)

# 指定遺伝子の存在確認（大文字小文字を無視して照合）
gene_map <- setNames(rownames(mat), toupper(rownames(mat)))
query_upper <- toupper(genes_of_interest)
present <- gene_map[query_upper]
present <- present[!is.na(present)]  # 存在するものだけ

missing <- genes_of_interest[!query_upper %in% names(present)]
if (length(present) == 0L) {
  stop("指定した遺伝子が見つかりませんでした。遺伝子名を確認してください。")
}
if (length(missing) > 0L) {
  message("見つからなかった遺伝子: ", paste(missing, collapse = ", "))
}

# 抽出
submat <- mat[present, , drop = FALSE]

# 変換：正規化済みとはいえ、レンジを整えるためにlog1pを推奨
submat_log <- log1p(submat)

# 行方向にZスコア（遺伝子ごとに中心化・スケーリング）
if (scale_by_gene) {
  submat_scaled <- t(scale(t(submat_log)))
} else {
  submat_scaled <- submat_log
}

# 欠損/NaN/Inf の安全対処
submat_scaled[!is.finite(submat_scaled)] <- 0

# ====== pheatmap 版 ======
# 距離: ユークリッド、クラスタ法: ward.D2（好みに応じて変更）
pdf("heatmap_pheatmap.pdf", width = 8, height = max(5, nrow(submat_scaled)*0.25 + 2))
pheatmap(
  submat_scaled,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  clustering_distance_cols = "euclidean",
  annotation_col = ann_col,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Selected genes (row-scaled)"
)
dev.off()

# ====== ComplexHeatmap 版 ======
# 注釈の色（任意）
group_levels <- levels(ann_col$Group)
group_cols <- setNames(colorRampPalette(brewer.pal(8,"Set2"))(length(group_levels)), group_levels)
ha <- HeatmapAnnotation(
  df = ann_col,
  col = list(Group = group_cols),
  which = "column"
)

# 行/列クラスタリング（階層的クラスタ）
row_hclust <- hclust(dist(submat_scaled, method = "euclidean"), method = "ward.D2")
col_hclust <- hclust(dist(t(submat_scaled), method = "euclidean"), method = "ward.D2")

pdf("heatmap_ComplexHeatmap.pdf", width = 8, height = max(5, nrow(submat_scaled)*0.3 + 2))
Heatmap(
  submat_scaled,
  name = "Z",
  top_annotation = ha,
  cluster_rows = row_hclust,
  cluster_columns = col_hclust,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Hierarchical clustering heatmap",
  row_title = paste0("Genes (n=", nrow(submat_scaled), ")"),
  col = colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027"))
)
dev.off()

# ===== 参考：遺伝子リストを外部ファイルから読む場合 =====
# genes_of_interest <- readLines("gene_list.txt")
# 空行や重複の掃除:
# genes_of_interest <- unique(genes_of_interest[nzchar(genes_of_interest)])

# 解析したい遺伝子（例：適宜置き換えてください）
genes_of_interest <- c(
  "DDIT3","CFB","ARRDC3","FOS","ATF3","H1-3","H2BC3","SCIMP","BNIPL","LILRA4","CD69",
  "SDE2","IFI30","PAIP1P1","H3C13","SOCS3","H1-4","H3-3B","CDKN1B","H4C6","EXO1","BRCA1",
   "VCAN","SMAD7","HMOX1","MMP9")



suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggpubr); library(readr)
})

# === 出力先フォルダ ===
out_dir <- "per_gene_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

submat_log <-  log1p(submat)

# === 列順が 'Pre→Post' の場合（基本形） ===
stopifnot(ncol(submat_log) %% 2 == 0)
n_pairs <- ncol(submat_log) / 2
sample_info <- data.frame(
  Sample = colnames(submat_log),
  Status = factor(c(rep("pre", n_pairs), rep("post", n_pairs)), levels = c("pre","post")),
  ID     = factor(rep(paste0("ID", seq_len(n_pairs)), times = 2))
)

# --- 列名から Pre/Post を判定したい場合（任意・必要なら置き換え） ---
# 例: 列名が "PB1_pre","PB1_post" などのとき
# sample_info <- tibble(Sample = colnames(submat_log)) |>
#   mutate(Status = case_when(
#     str_detect(tolower(Sample), "pre")  ~ "pre",
#     str_detect(tolower(Sample), "post") ~ "post",
#     TRUE ~ NA_character_
#   )) |>
#   tidyr::extract(Sample, into = "ID", regex = "(PB\\d+)", remove = FALSE) |>
#   mutate(Status = factor(Status, levels = c("pre","post")),
#          ID = factor(ID))
# stopifnot(!any(is.na(sample_info$Status)), !any(is.na(sample_info$ID)))

# === ロング化（1遺伝子=1ファセット想定だが、後で1遺伝子ずつ抽出） ===
expr_long <- submat_log |>
  as.data.frame() |>
  rownames_to_column("Gene") |>
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expr") |>
  left_join(sample_info, by = "Sample") |>
  filter(is.finite(Expr))

# === 実際に存在する遺伝子のみ ===
genes_ok <- unique(intersect(expr_long$Gene, present))

# === 保存関数 ===
save_one_gene_plot <- function(gene_name, data_gene) {
  # ペアド t 検定
  wide <- data_gene |>
    select(ID, Status, Expr) |>
    pivot_wider(names_from = Status, values_from = Expr) |>
    arrange(ID)
  t_res <- try(t.test(wide$pre, wide$post, paired = TRUE), silent = TRUE)
  pval  <- if (inherits(t_res, "try-error")) NA_real_ else t_res$p.value
  delta <- if (inherits(t_res, "try-error")) NA_real_ else mean(wide$post - wide$pre, na.rm = TRUE)
  n     <- nrow(wide)
  
  # プロット
  p <- ggplot(data_gene, aes(x = Status, y = Expr, group = ID, fill = Status)) +
    geom_line(color = "gray60", alpha = 0.6) +
    geom_point(size = 3, pch = 21) +
    scale_fill_manual(values = c("pre" = "#e41a1c", "post" = "#377eb8")) +
    theme_classic(base_size = 13) +
    labs(
      title = gene_name,
      subtitle = sprintf("paired t-test: p=%s",
                          signif(pval, 3)),
      x = NULL, y = "Log(Normalized counts)"
    ) +
    theme(legend.position = "none") +
    #stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
    #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, color = "black") +
    # 星表示（視覚用）：数値は subtitle に出してあるので任意
    stat_compare_means(paired = TRUE, method = "t.test", label = "p.signif")
  
  # ファイル名を安全化
  safe <- gsub("[^A-Za-z0-9._-]+", "_", gene_name)
  ggsave(file.path(out_dir, paste0(safe, "_paired_lineplot.png")), p, width = 4.5, height = 4.2, dpi = 300)
  ggsave(file.path(out_dir, paste0(safe, "_paired_lineplot.pdf")), p, width = 4.5, height = 4.2, device = cairo_pdf)
  
  # 数値結果も書き出し（任意）：個票 + 要約
  readr::write_tsv(
    wide |>
      mutate(delta = post - pre),
    file.path(out_dir, paste0(safe, "_paired_values.tsv"))
  )
  readr::write_tsv(
    tibble(Gene = gene_name, n = n, delta_mean = delta, p_value = pval),
    file.path(out_dir, paste0(safe, "_paired_summary.tsv"))
  )
}

# === ループ実行 ===
message(sprintf("Exporting %d genes into '%s' ...", length(genes_ok), out_dir))
for (g in genes_ok) {
  dg <- expr_long |> filter(Gene == g)
  if (n_distinct(dg$Status) == 2 && n_distinct(dg$ID) >= 2) {
    save_one_gene_plot(g, dg)
  } else {
    warning(sprintf("Skip '%s' (need both pre & post and >=2 pairs).", g))
  }
}
message("Done.")

