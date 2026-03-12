# 必要パッケージ
# install.packages(c("data.table","dplyr","tibble","stringr","pheatmap"))
# install.packages("msigdbr")
# BiocManager::install(c("GSVA","GSEABase","ComplexHeatmap","circlize"))

library(data.table); library(dplyr); library(tibble); library(stringr)
library(GSVA); library(msigdbr)
library(pheatmap); library(ComplexHeatmap); library(circlize)

# ========= 0) 入力（正規化カウント） =========
path_input  <- "C://Users/hikob/Dropbox/Research/PAI-I_inhibitor_project/250821_RNAseq/"
path_output <- file.path(path_input, "DESeq2"); dir.create(path_output, showWarnings = FALSE)
expr <- fread(file.path(path_input, "Counts.normalized.txt")) %>% as.data.frame()
colnames(expr)[1] <- "Gene"
rownames(expr) <- make.unique(expr$Gene); expr$Gene <- NULL
expr <- as.matrix(expr); storage.mode(expr) <- "double"

# 軽いフィルタ & log 変換（GSVA は rank-based ではないので振幅整え推奨）
keep <- rowSums(expr >= 1) >= ceiling(ncol(expr)/2)
expr <- expr[keep, , drop = FALSE]
expr_log <- log1p(expr)

# 行名を大文字化＆重複集約（gene set 名と頑健に合わせるため）
expr_log_u <- expr_log; rownames(expr_log_u) <- toupper(rownames(expr_log_u))
if (any(duplicated(rownames(expr_log_u)))) {
  expr_log_u <- rowsum(expr_log_u, group = rownames(expr_log_u), reorder = FALSE)
}

# ========= 1) C2 から目的セットだけ抽出 =========
clean <- function(x) toupper(trimws(x))
target_sets_raw <- c(
  "WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP",
  "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE",
  "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
  "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
  "REACTOME_DNA_DOUBLE_STRAND_BREAK_RESPONSE",
  "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE",
  "JAATINEN_HEMATOPOIETIC_STEM_CELL_DN",
  "GERY_CEBP_TARGETS",
  "REACTOME_REGULATION_OF_PTEN_GENE_TRANSCRIPTION",
  "MCBRYAN_PUBERTAL_TGFB1_TARGETS_UP",
  "WONG_ADULT_TISSUE_STEM_MODULE ",
  "LEE_BMP2_TARGETS_UP",
  "LEE_NEURAL_CREST_STEM_CELL_UP",
  "BOQUEST_STEM_CELL_CULTURED_VS_FRESH_UP"
)
target_sets <- clean(target_sets_raw)

m_c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  mutate(gs_name_clean = clean(gs_name), gene_symbol_clean = clean(gene_symbol))

m_sel <- m_c2 %>% filter(gs_name_clean %in% target_sets)
genesets_C2_selected <- split(m_sel$gene_symbol_clean, m_sel$gs_name_clean)
genesets_C2_selected <- lapply(genesets_C2_selected, \(v) unique(na.omit(v)))

# ========= 2) GSVA（method = "gsva"） =========
param <- gsvaParam(
  exprData = expr_log_u,
  geneSets = genesets_C2_selected,
  kcdf     = "Gaussian",
  minSize  = 10,
  maxSize  = 5000,
  tau      = 1,
  maxDiff  = TRUE,   # = 旧mx.diff
  absRanking = FALSE # GSVA法では通常FALSEのままでOK
)

es <- gsva(param)  # 行=Pathway, 列=Sample のスコア行列
# es: 行=Pathway, 列=Sample

# ========= 3) 可視化（ヒートマップ） =========
# 行方向Zスケール
esZ <- t(scale(t(es))); esZ[!is.finite(esZ)] <- 0

# サンプル注釈（Pre/Post がある前提。無ければ適宜作成）
if (!exists("sampleData")) {
  smp <- colnames(esZ)
  class <- c(rep("Pre",3),rep("Post",3))
    
    
  sampleData <- data.frame(samples = smp, class = factor(class, levels = c("Pre","Post","Other")))
}
ann_col <- data.frame(Group = droplevels(sampleData$class[match(colnames(esZ), sampleData$samples)]))
rownames(ann_col) <- colnames(esZ)

# pheatmap（手軽）
pdf(file.path(path_output, "GSVA_C2_selected_heatmap_pheatmap.pdf"),
    width = 8, height = max(5, nrow(esZ)*0.25 + 2))
pheatmap(esZ,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",
         annotation_col = ann_col,
         main = "GSVA (C2 selected)",
         fontsize_row = 8, fontsize_col = 8)
dev.off()

# ========= 4) 出力 =========
fwrite(as.data.frame(es)  %>% rownames_to_column("Pathway"),
       file.path(path_output, "GSVA_scores_C2_selected.tsv"), sep = "\t")
fwrite(as.data.frame(esZ) %>% rownames_to_column("Pathway"),
       file.path(path_output, "GSVA_scores_C2_selected_rowZ.tsv"), sep = "\t")
