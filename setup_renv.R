#!/usr/bin/env Rscript
# =============================================================================
#  setup_renv.R — R パッケージ環境のセットアップ (Windows / Mac / Linux 共通)
#
#  使い方 (ターミナル or RStudio コンソールから):
#    Rscript setup_renv.R          # renv.lock に従ってパッケージを復元
#    Rscript setup_renv.R --update # 新パッケージ追加後にlockfileを更新
#
#  初回セットアップの流れ:
#    1. git clone https://github.com/takubo-lab/RNAseq_pipeline_takubo.git
#    2. cd RNAseq_pipeline_takubo
#    3. Rscript setup_renv.R
#    4. → renv.lock に記載された全パッケージが自動インストールされる
#
#  R 解析のみ (Windows/Mac) の場合:
#    Step 4 (DESeq2), Step 5 (fGSEA), Step 6 (DEXSeq) の R スクリプトを
#    RStudio や Rscript から直接実行可能
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# --- renv がインストールされていなければインストール ---
if (!requireNamespace("renv", quietly = TRUE)) {
  cat("[INFO] renv をインストール中...\n")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

cat("=== RNA-seq Pipeline — R 環境セットアップ ===\n\n")

if ("--update" %in% args) {
  # --- lockfile 更新モード ---
  cat("[INFO] renv.lock を更新します...\n")

  # Bioconductor リポジトリを設定
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    renv::install("BiocManager")
  }
  options(repos = BiocManager::repositories())

  renv::snapshot()
  cat("\n[INFO] renv.lock を更新しました。\n")
  cat("[INFO] git commit して共有してください。\n")

} else {
  # --- 復元モード (デフォルト) ---
  cat("[INFO] renv.lock に従ってパッケージを復元します...\n")
  cat("[INFO] 初回は 10〜30 分かかる場合があります。\n\n")

  renv::restore(prompt = FALSE)

  # 検証
  cat("\n=== パッケージ検証 ===\n")
  required <- c(
    "DESeq2", "DEXSeq", "fgsea", "msigdbr", "BiocParallel",
    "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
    "tidyverse", "data.table", "magrittr", "umap", "yaml"
  )

  all_ok <- TRUE
  for (pkg in required) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      ver <- as.character(packageVersion(pkg))
      cat(sprintf("  OK : %-20s %s\n", pkg, ver))
    } else {
      cat(sprintf("  NG : %-20s (not found)\n", pkg))
      all_ok <- FALSE
    }
  }

  if (all_ok) {
    cat("\n=== 全パッケージの復元が完了しました ===\n")
  } else {
    cat("\n[WARN] 一部パッケージが見つかりません。\n")
    cat("  renv::install('パッケージ名') で個別にインストールしてください。\n")
  }
}

cat("\n使い方:\n")
cat("  # RStudio: プロジェクトを開くだけで renv が自動アクティベート\n")
cat("  # ターミナル: cd RNAseq_pipeline_takubo && Rscript scripts/04_deseq2.R\n")
