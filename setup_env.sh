#!/usr/bin/env bash
# =============================================================================
#  setup_env.sh — RNA-seq パイプライン環境セットアップ (Linux / Mac)
#
#  使い方:
#    bash setup_env.sh          # 環境を新規作成
#    bash setup_env.sh --update # 既存環境を更新
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="rnaseq_takubo"
ENV_FILE="${SCRIPT_DIR}/environment.yml"

# --- 色付きメッセージ ---
info()  { echo -e "\033[1;34m[INFO]\033[0m  $*"; }
warn()  { echo -e "\033[1;33m[WARN]\033[0m  $*"; }
error() { echo -e "\033[1;31m[ERROR]\033[0m $*"; exit 1; }

cleanup() {
  if [[ -n "${TMP_ENV_FILE:-}" && -f "${TMP_ENV_FILE}" ]]; then
    rm -f "${TMP_ENV_FILE}"
  fi
  if [[ -n "${TMP_STAR_BUILD_DIR:-}" && -d "${TMP_STAR_BUILD_DIR}" ]]; then
    rm -rf "${TMP_STAR_BUILD_DIR}"
  fi
}

trap cleanup EXIT

# --- conda / mamba 検出 ---
detect_installer() {
  if command -v mamba &>/dev/null; then
    echo "mamba"
  elif command -v conda &>/dev/null; then
    echo "conda"
  else
    echo ""
  fi
}

install_miniforge() {
  info "Miniforge (conda + mamba) をインストールします..."
  local installer="Miniforge3-$(uname)-$(uname -m).sh"
  wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/${installer}" \
    -O /tmp/${installer}
  bash /tmp/${installer} -b -p "${HOME}/miniforge3"
  rm /tmp/${installer}

  # シェル初期化
  eval "$(${HOME}/miniforge3/bin/conda shell.bash hook)"
  conda init bash 2>/dev/null || true
  conda init zsh  2>/dev/null || true
  info "Miniforge をインストールしました。"
  info "シェルを再起動してから、再度このスクリプトを実行してください。"
  exit 0
}

make_env_file_without_star() {
  TMP_ENV_FILE="$(mktemp --suffix=.yml)"
  awk '
    /^[[:space:]]*-[[:space:]]*star([[:space:]]|=|$)/ { next }
    { print }
  ' "${ENV_FILE}" > "${TMP_ENV_FILE}"
  echo "${TMP_ENV_FILE}"
}

get_star_version() {
  local star_line
  star_line="$(awk '
    /^[[:space:]]*-[[:space:]]*star=/ {
      sub(/^[[:space:]]*-[[:space:]]*star=/, "", $0)
      print $0
      exit
    }
  ' "${ENV_FILE}")"

  if [[ -z "${star_line}" ]]; then
    error "environment.yml から STAR のバージョンを取得できません。"
  fi

  printf '%s\n' "${star_line}"
}

check_build_tools() {
  local missing=()
  for tool in make g++ tar; do
    if ! command -v "${tool}" &>/dev/null; then
      missing+=("${tool}")
    fi
  done

  if [[ ${#missing[@]} -gt 0 ]]; then
    error "STARソースビルドに必要なツールが不足しています: ${missing[*]}"
  fi
}

download_star_source() {
  local tarball="${1}"
  local star_source_url="${2}"

  if command -v wget &>/dev/null; then
    wget -q -O "${tarball}" "${star_source_url}"
  elif command -v curl &>/dev/null; then
    curl -L "${star_source_url}" -o "${tarball}"
  else
    error "STARソース取得には wget または curl が必要です。"
  fi
}

install_star_from_source() {
  local star_version="${1}"
  local star_source_url="https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz"

  info "STAR ${star_version} を公式ソースからインストールします..."

  eval "$(conda shell.bash hook)"
  conda activate "${ENV_NAME}"
  check_build_tools

  TMP_STAR_BUILD_DIR="$(mktemp -d)"
  local tarball="${TMP_STAR_BUILD_DIR}/${star_version}.tar.gz"
  local star_root="${TMP_STAR_BUILD_DIR}/STAR-${star_version}"
  local make_target="STAR"
  local -a make_args=()

  download_star_source "${tarball}" "${star_source_url}"
  tar -xzf "${tarball}" -C "${TMP_STAR_BUILD_DIR}"

  if [[ ! -d "${star_root}/source" ]]; then
    error "STARソース展開に失敗しました: ${star_root}/source"
  fi

  if conda list -n "${ENV_NAME}" 2>/dev/null | awk '{print $1}' | grep -qx 'star'; then
    info "既存の conda 版 STAR を削除します..."
    conda remove -n "${ENV_NAME}" -y star || true
  fi

  if [[ "$(uname)" == "Linux" ]] && ! grep -qm1 '\bavx\b' /proc/cpuinfo; then
    make_args+=("CXXFLAGS_SIMD=sse")
  fi

  if [[ "$(uname)" == "Darwin" ]]; then
    make_target="STARforMacStatic"
    local cxx="${STAR_CXX:-}"
    if [[ -z "${cxx}" ]]; then
      cxx="$(command -v g++-14 || command -v g++-13 || command -v g++-12 || command -v g++-11 || command -v g++-10 || command -v g++-9 || command -v g++)"
    fi
    [[ -n "${cxx}" ]] || error "macOSでは GNU g++ が必要です。STAR_CXX で指定してください。"
    make_args+=("CXX=${cxx}")
  fi

  info "STAR をビルド中..."
  (
    cd "${star_root}/source"
    make "${make_target}" "${make_args[@]}"
  )

  install -m 0755 "${star_root}/source/STAR" "${CONDA_PREFIX}/bin/STAR"

  if [[ -f "${star_root}/source/STARlong" ]]; then
    install -m 0755 "${star_root}/source/STARlong" "${CONDA_PREFIX}/bin/STARlong"
  fi

  if ! "${CONDA_PREFIX}/bin/STAR" --version >/dev/null 2>&1; then
    error "ソースビルドした STAR の動作確認に失敗しました。"
  fi

  info "STAR ${star_version} のソースインストール完了"
}

# --- ツール検証 ---
verify_tools() {
  info "インストールされたツールを検証中..."
  eval "$(conda shell.bash hook)"
  conda activate "${ENV_NAME}"

  local tools=("STAR" "samtools" "featureCounts" "Rscript")
  local all_ok=true
  for tool in "${tools[@]}"; do
    if command -v "${tool}" &>/dev/null; then
      local ver=""
      case "${tool}" in
        STAR)          ver=$(STAR --version 2>&1 || echo "OK") ;;
        samtools)      ver=$(samtools --version 2>&1 | head -1 | awk '{print $2}') ;;
        featureCounts) ver=$(featureCounts -v 2>&1 | head -1 || echo "OK") ;;
        Rscript)       ver=$(Rscript --version 2>&1 | grep -oP '[\d.]+' | head -1 || echo "OK") ;;
      esac
      info "  ✓ ${tool}: ${ver}"
    else
      warn "  ✗ ${tool}: 見つかりません"
      all_ok=false
    fi
  done

  # R パッケージ検証
  info "R パッケージを検証中..."
  Rscript -e '
    pkgs <- c("DESeq2", "DEXSeq", "fgsea", "GSVA", "msigdbr", "BiocParallel",
              "tidyverse", "data.table", "ggplot2", "ggrepel", "umap",
              "magrittr", "yaml", "pheatmap")
    for (p in pkgs) {
      ok <- requireNamespace(p, quietly = TRUE)
      if (ok) {
        ver <- as.character(packageVersion(p))
        cat(sprintf("  OK : %-20s %s\n", p, ver))
      } else {
        cat(sprintf("  NG : %-20s (not found)\n", p))
      }
    }
  ' 2>/dev/null || warn "R パッケージの検証中にエラーが発生しました。"

  if ${all_ok}; then
    info "全ツールの検証完了 ✓"
  else
    warn "一部ツールが見つかりません。手動で確認してください。"
  fi
}

# --- メイン ---
main() {
  info "=== RNA-seq Pipeline 環境セットアップ ==="
  info "環境名: ${ENV_NAME}"
  info "定義ファイル: ${ENV_FILE}"
  echo ""

  # 1. conda/mamba 確認
  INSTALLER=$(detect_installer)
  if [[ -z "${INSTALLER}" ]]; then
    warn "conda/mamba が見つかりません。"
    read -rp "Miniforge を自動インストールしますか？ [y/N]: " yn
    case "${yn}" in
      [Yy]*) install_miniforge ;;
      *)     error "conda または mamba をインストールしてから再実行してください。" ;;
    esac
  fi
  info "パッケージマネージャ: ${INSTALLER}"

  local resolved_env_file
  local star_version
  resolved_env_file="$(make_env_file_without_star)"
  star_version="$(get_star_version)"
  info "environment.yml から conda 版 STAR を除外して環境を構築します。"

  # 2. 既存環境の確認
  if conda env list 2>/dev/null | grep -qw "${ENV_NAME}"; then
    warn "環境 '${ENV_NAME}' は既に存在します。"
    if [[ "${1:-}" == "--update" ]]; then
      info "環境を更新します..."
      ${INSTALLER} env update -n "${ENV_NAME}" -f "${resolved_env_file}" --prune
      install_star_from_source "${star_version}"
      verify_tools
      info "=== 環境の更新が完了しました ==="
      info "使用方法: conda activate ${ENV_NAME}"
      exit 0
    fi
    read -rp "再作成しますか（既存環境を削除）？ [y/N]: " yn
    case "${yn}" in
      [Yy]*)
        info "既存環境を削除中..."
        conda env remove -n "${ENV_NAME}" -y
        ;;
      *)
        info "既存環境を更新します..."
        ${INSTALLER} env update -n "${ENV_NAME}" -f "${resolved_env_file}" --prune
        install_star_from_source "${star_version}"
        verify_tools
        info "=== 環境の更新が完了しました ==="
        info "使用方法: conda activate ${ENV_NAME}"
        exit 0
        ;;
    esac
  fi

  # 3. 環境作成
  info "環境 '${ENV_NAME}' を作成中（10〜20分かかる場合があります）..."
  ${INSTALLER} env create -f "${resolved_env_file}"

  # 4. STAR を公式ソースから導入
  install_star_from_source "${star_version}"

  # 5. 検証
  verify_tools

  echo ""
  info "=== セットアップ完了 ==="
  info ""
  info "使用方法:"
  info "  conda activate ${ENV_NAME}"
  info "  bash run_pipeline.sh"
}

main "$@"
