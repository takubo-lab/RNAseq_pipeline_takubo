# RNAseq Pipeline (Takubo Lab)

RNA-seq解析パイプライン。FASTQ → STAR mapping → featureCounts → DESeq2 → fGSEA、およびエクソンレベルのsplicing variant解析を一括実行。

## パイプライン概要

| Step | Script | 内容 |
|------|--------|------|
| 1 | `scripts/01_mapping.sh` | STARによるリファレンスゲノムへのマッピング |
| 2 | `scripts/02_featurecounts.sh` | featureCounts (遺伝子レベルカウント) |
| 3 | `scripts/03_featurecounts_exon.sh` | featureCounts (エクソンレベルカウント) |
| 4 | `scripts/04_deseq2.R` | DESeq2 発現変動解析 (全グループペア比較) |
| 5 | `scripts/05_fgsea.R` | fGSEA 遺伝子セットエンリッチメント解析 |
| 6 | `scripts/06_exon_de.R` | DEXSeq エクソンレベルDE (スプライシング変異解析) |

## クイックスタート

### 1. セットアップ

```bash
# リポジトリをクローン
git clone https://github.com/takubo-lab/RNAseq_pipeline_takubo.git
cd RNAseq_pipeline_takubo
```

### 2-A. プロジェクト初期化（推奨）

パイプラインコードとデータを分離して管理する方法:

```bash
# 新しいプロジェクトを初期化
bash init_project.sh /path/to/my_project hg38

# サンプル情報を編集
vim /path/to/my_project/samples.tsv

# FASTQファイルを配置
cp *.fq.gz /path/to/my_project/fastq/

# パイプライン実行
bash run_pipeline.sh --config /path/to/my_project/config.sh
```

### 2-B. 直接実行（シンプル）

パイプラインディレクトリ内で直接解析する方法:

### 2-B. 直接実行（シンプル）

パイプラインディレクトリ内で直接解析する方法:

**`config.sh`** — ゲノム、パス、パラメータを設定:
```bash
GENOME="hg38"              # "hg38" or "mm10"
FEATURE_COUNTS_PAIRED=true  # paired-end or single-end
```

**`samples.tsv`** — サンプル情報を記載:
```
sample_name	group	fq_prefix	lane_suffix
Sample1_pre	pre	Sample1	L5
Sample1_post	post	Sample1_post	L5
```

### 3. FASTQファイルを配置

```bash
mkdir fastq
# fastq/ ディレクトリにfq.gzファイルを配置
# ファイル名: {fq_prefix}_{lane_suffix}_{1,2}.fq.gz
```

### 4. パイプライン実行

```bash
# 全ステップ実行
bash run_pipeline.sh

# 特定のステップから開始
bash run_pipeline.sh --from 4

# 特定のステップのみ実行
bash run_pipeline.sh --only 5

# プロジェクト指定で実行
bash run_pipeline.sh --config /path/to/my_project/config.sh

# 複数 count matrix を simple merge して DESeq2 以降を実行
bash run_pipeline.sh --from 4 \
	--merge-samples /path/a/samples.tsv,/path/b/samples.tsv \
	--merge-counts /path/a/counts/count.txt,/path/b/counts/count.txt \
	--merge-exon-counts /path/a/counts/exon_count.txt,/path/b/counts/exon_count.txt
```

`--merge-samples`, `--merge-counts`, `--merge-exon-counts` は Step 4-6 で使用されます。
同じ annotation / gene definition で作成した count matrix を前提に、行名で単純結合し、サンプル列を横方向にマージします。

## バージョン管理・マルチユーザー運用

### バージョニング

パイプラインは `VERSION` ファイルでセマンティックバージョニングを管理:
- パッチ (1.0.x): バグ修正
- マイナー (1.x.0): 機能追加（後方互換）
- メジャー (x.0.0): 破壊的変更

リリース時は git tag を使用:
```bash
git tag -a v1.1.0 -m "Add comparisons.tsv support"
git push origin v1.1.0
```

### Provenance (実行記録)

パイプライン実行時に `provenance.yml` が自動生成され、以下を記録:
- パイプラインバージョン・コミット
- 実行日時・ユーザー・ホスト名
- ソフトウェアバージョン (STAR, R, featureCounts)

### ロック機構

同一プロジェクトでの同時実行を防止する `.pipeline.lock` ファイルによるロック機構を搭載。

## ディレクトリ構成

```
RNAseq_pipeline_takubo/
├── VERSION                # パイプラインバージョン
├── config.sh              # パラメータ設定
├── samples.tsv            # サンプル情報
├── init_project.sh        # プロジェクト初期化スクリプト
├── run_pipeline.sh        # メイン実行スクリプト
├── scripts/
│   ├── utils.sh           # ユーティリティ関数
│   ├── 01_mapping.sh
│   ├── 02_featurecounts.sh
│   ├── 03_featurecounts_exon.sh
│   ├── 04_deseq2.R
│   ├── 05_fgsea.R
│   └── 06_exon_de.R
├── Gene_set/              # カスタム遺伝子セット
│   └── HSC_set3.gmt
├── fastq/                 # 入力FASTQファイル (gitignore)
├── BAM/                   # マッピング結果 (gitignore)
├── counts/                # カウントデータ (gitignore)
├── DESeq2/                # DESeq2結果
├── fGSEA/                 # fGSEA結果
└── exon_DE/               # エクソンレベルDE結果
```

## 出力ファイル

### DESeq2
- `DESeq2/Counts.normalized.txt` — 正規化カウント
- `DESeq2/{g1}_vs_{g2}_DE_Results.tsv` — 発現変動解析結果
- `DESeq2/{g1}_vs_{g2}_DE_Results_with_paired.tsv` — 対応ありt検定結果を含む
- `DESeq2/PCA_class.png` — PCAプロット
- `DESeq2/UMAP_all.png` — UMAPプロット
- `DESeq2/*_volcano.png` — ボルケーノプロット

### fGSEA
- `fGSEA/fGSEA_{comparison}/fgsea_Results_*.tsv` — 各カテゴリのGSEA結果

### Exon-level DE
- `exon_DE/DEXSeq_all_results.tsv` — 全エクソン結果
- `exon_DE/DEXSeq_significant_exons.tsv` — 有意なエクソン
- `exon_DE/DEXSeq_gene_summary.tsv` — 遺伝子ごとのサマリー
- `exon_DE/plots/` — トップ遺伝子のエクソン使用率プロット

## 必要なソフトウェア

### コマンドラインツール
- [STAR](https://github.com/alexdobin/STAR) (v2.7.3+)
- [Subread/featureCounts](http://subread.sourceforge.net/)

### Rパッケージ
```r
# Bioconductor
BiocManager::install(c("DESeq2", "fgsea", "msigdbr", "DEXSeq", "BiocParallel"))

# CRAN
install.packages(c("tidyverse", "data.table", "magrittr", "umap", "ggrepel"))
```

## 対応ゲノム

| Genome | STAR Index | GTF |
|--------|-----------|-----|
| hg38 | `~/STAR2.7.3/hg38_sjdb149_RefSeq` | `~/Genome/hg38/hg38.ncbiRefSeq.gtf` |
| mm10 | `~/STAR2.7.3/mm10_sjdb149_RefSeq_ver241129` | `~/Genome/mm10/mm10.ncbiRefSeq.gtf` |

パスは `config.sh` で変更可能。

## License

Internal use — Takubo Lab
