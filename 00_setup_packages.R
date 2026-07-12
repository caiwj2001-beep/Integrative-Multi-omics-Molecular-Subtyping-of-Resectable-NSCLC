# =============================================================
# NSCLC多组学分子分型项目 - R包安装脚本
# Project: Resectable NSCLC Multi-omics Molecular Subtyping
# =============================================================

# 用户级库目录(避免Program Files权限问题)
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

cat("Library paths:\n"); print(.libPaths())

# ---- CRAN 包 ----
cran_pkgs <- c(
  "data.table",      # 快速IO
  "survival",        # 生存分析
  "survminer",       # KM曲线
  "rms",             # 回归建模
  "ggplot2",         # 绘图
  "dplyr", "tidyr",  # 数据处理
  "readr",
  "pheatmap",        # 热图
  "RColorBrewer",
  "ComplexHeatmap",  # 复杂热图(也在Bioc)
  "circlize",
  "cluster",         # 聚类
  "factoextra",
  "umap",
  "Rtsne",
  "matrixStats",
  "corrplot",
  "ggpubr",
  "cowplot",
  "gridExtra",
  "openxlsx",        # Excel输出
  "officer",         # DOCX输出
  "flextable",       # DOCX表格
  "BiocManager"
)

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      cat(sprintf("Installing CRAN: %s\n", p))
      tryCatch(install.packages(p, lib = user_lib, dependencies = TRUE),
               error = function(e) cat(sprintf("  FAILED %s: %s\n", p, conditionMessage(e))))
    } else {
      cat(sprintf("  OK: %s\n", p))
    }
  }
}

cat("\n=== Installing CRAN packages ===\n")
install_if_missing(cran_pkgs)

# ---- Bioconductor 包 ----
cat("\n=== Installing Bioconductor packages ===\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = user_lib)

bioc_pkgs <- c(
  "maftools",              # 突变景观
  "ConsensusClusterPlus",  # 共识聚类
  "GSVA",                  # 基因集变异分析
  "clusterProfiler",       # 富集分析
  "org.Hs.eg.db",          # 基因注释
  "limma",                 # 差异表达
  "edgeR",                 # RNA-seq DE
  "DESeq2",
  "estimate",              # 免疫基质评分 (可能需其他源)
  "ComplexHeatmap",
  "MOFA2",                 # 多组学整合
  "sva"                    # 批次校正
)

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("Installing Bioc: %s\n", p))
    tryCatch(BiocManager::install(p, lib = user_lib, update = FALSE, ask = FALSE),
             error = function(e) cat(sprintf("  FAILED %s: %s\n", p, conditionMessage(e))))
  } else {
    cat(sprintf("  OK: %s\n", p))
  }
}

cat("\n=== Package installation complete ===\n")

# 验证核心包
core <- c("data.table","survival","survminer","maftools","GSVA",
          "ConsensusClusterPlus","ComplexHeatmap","MOFA2","limma","officer","flextable")
cat("\n=== Core package availability ===\n")
for (p in core) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %s: %s\n", p, ifelse(ok, "AVAILABLE", "MISSING")))
}
