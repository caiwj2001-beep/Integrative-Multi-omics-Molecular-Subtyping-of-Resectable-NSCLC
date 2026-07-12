# =============================================================
# 04_immune_pathway.R
# 免疫微环境(ESTIMATE) + 通路活性(ssGSEA/GSVA) + 免疫细胞浸润
# =============================================================
source("config.R")
suppressMessages({
  library(data.table); library(GSVA); library(matrixStats)
})

# ---- 读取处理后的表达 ----
ep <- readRDS(file.path(RESULTS, "expression_processed.rds"))
logmat <- ep$logmat  # log2(TPM+1), gene_id rows
id2name <- ep$gene_id2name

# 转成gene_symbol矩阵(取每symbol方差最大的id)
sym <- id2name[rownames(logmat)]
dt <- data.table(sym=sym, idx=seq_len(nrow(logmat)), v=rowVars(logmat))
dt <- dt[!is.na(sym) & sym != ""]
setorder(dt, sym, -v)
keep <- dt[, .SD[1], by=sym]
symmat <- logmat[keep$idx, ]
rownames(symmat) <- keep$sym
cat("Symbol-level expression matrix:", nrow(symmat), "x", ncol(symmat), "\n")

# ================= ESTIMATE =================
cat("\n=== Running ESTIMATE ===\n")
tryCatch({
  library(estimate)
  # estimate需要写入文件
  est_in <- file.path(RESULTS, "estimate_input.gct")
  tmp_tsv <- file.path(RESULTS, "estimate_expr.txt")
  df_out <- data.frame(GeneSymbol=rownames(symmat), symmat, check.names=FALSE)
  write.table(df_out, tmp_tsv, sep="\t", quote=FALSE, row.names=FALSE)
  filterCommonGenes(input.f=tmp_tsv, output.f=est_in, id="GeneSymbol")
  est_out <- file.path(RESULTS, "estimate_score.gct")
  estimateScore(est_in, est_out, platform="illumina")
  est <- read.table(est_out, skip=2, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  est_scores <- as.data.frame(t(est[,-1]))
  est_scores$patient_id <- rownames(est_scores)
  fwrite(est_scores, file.path(RESULTS, "estimate_scores.csv"))
  cat("ESTIMATE done. Columns:", paste(colnames(est_scores), collapse=", "), "\n")
}, error=function(e) cat("ESTIMATE error:", conditionMessage(e), "\n"))

# ================= ssGSEA: 免疫细胞 + Hallmark通路 =================
cat("\n=== ssGSEA immune cell signatures ===\n")
# Bindea et al. 28免疫细胞签名 (常用, 内置一个精简版)
# 使用msigdbr获取Hallmark基因集
tryCatch({
  library(msigdbr)
  hallmark <- msigdbr(species="Homo sapiens", category="H")
  hm_list <- split(hallmark$gene_symbol, hallmark$gs_name)
  cat("Hallmark gene sets:", length(hm_list), "\n")

  # ssGSEA (GSVA新API)
  gsva_par <- ssgseaParam(exprData=symmat, geneSets=hm_list)
  hm_scores <- gsva(gsva_par, verbose=FALSE)
  hm_dt <- as.data.table(t(hm_scores), keep.rownames="patient_id")
  fwrite(hm_dt, file.path(RESULTS, "hallmark_ssgsea.csv"))
  cat("Hallmark ssGSEA done:", nrow(hm_scores), "pathways x", ncol(hm_scores), "samples\n")
}, error=function(e) cat("Hallmark ssGSEA error:", conditionMessage(e), "\n"))

# ================= 免疫检查点 & 关键标志基因表达 =================
immune_genes <- c("CD274","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","CD8A","CD8B",
                  "GZMB","PRF1","IFNG","CXCL9","CXCL10","CXCL13","FOXP3","CD4",
                  "TGFB1","IL10","CD68","ITGAX","HLA-A","HLA-B","B2M","TAP1")
present <- intersect(immune_genes, rownames(symmat))
ig_dt <- as.data.table(t(symmat[present, ]), keep.rownames="patient_id")
fwrite(ig_dt, file.path(RESULTS, "immune_gene_expression.csv"))
cat("\nImmune marker genes extracted:", length(present), "\n")

cat("\n=== Immune/pathway analysis complete ===\n")
