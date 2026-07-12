# =============================================================
# 02_expression_clustering.R
# 表达谱共识聚类 (ConsensusClusterPlus) - 分型主轴
# =============================================================
source("config.R")
suppressMessages({
  library(data.table); library(ConsensusClusterPlus); library(matrixStats)
})
set.seed(20260707)

# ---- 读取TPM矩阵 ----
cat("Reading TPM matrix...\n")
tpm <- fread(file.path(MATDIR, "expr_tpm.tsv.gz"))
gene_names <- tpm$gene_name
gene_ids <- tpm$gene_id
mat <- as.matrix(tpm[, -c("gene_id","gene_name")])
rownames(mat) <- gene_ids
cat("TPM matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")

# ---- log2转换 ----
logmat <- log2(mat + 1)

# ---- 选高变基因 (top 2000 MAD) ----
mads <- rowMads(logmat)
names(mads) <- rownames(logmat)
top_genes <- names(sort(mads, decreasing = TRUE))[1:2000]
sub <- logmat[top_genes, ]
# 中位中心化
sub <- sweep(sub, 1, rowMedians(sub), "-")
cat("High-variance genes selected:", length(top_genes), "\n")

# ---- 共识聚类 ----
cat("Running ConsensusClusterPlus...\n")
cc_dir <- file.path(RESULTS, "cc_expression")
dir.create(cc_dir, showWarnings = FALSE)
# pam + pearson 是TCGA表达分型标准, 比km快很多且更稳健
cc <- ConsensusClusterPlus(sub, maxK = 8, reps = 500, pItem = 0.8, pFeature = 1,
                           clusterAlg = "pam", distance = "pearson",
                           seed = 20260707, plot = "png", title = cc_dir,
                           writeTable = FALSE)

# ---- 评估最优K: 计算共识CDF的PAC (Proportion of Ambiguous Clustering) ----
calc_PAC <- function(cc_res, k, x1 = 0.1, x2 = 0.9) {
  M <- cc_res[[k]]$consensusMatrix
  vals <- M[lower.tri(M)]
  Fn <- ecdf(vals)
  Fn(x2) - Fn(x1)
}
pac <- sapply(2:8, function(k) calc_PAC(cc, k))
names(pac) <- 2:8
cat("\n=== PAC by K (lower=better) ===\n"); print(round(pac, 4))
best_k <- as.integer(names(which.min(pac)))
cat("Best K by PAC:", best_k, "\n")

# 也看cluster稳定性: 通常选3-5之间平衡可解释性
# 输出多个K的分配供后续选择
for (k in 2:6) {
  cl <- cc[[k]]$consensusClass
  tbl <- table(cl)
  cat(sprintf("K=%d cluster sizes: %s\n", k, paste(tbl, collapse=", ")))
}

# ---- 保存聚类分配 (K=2..6) ----
assign_dt <- data.table(patient_id = colnames(sub))
for (k in 2:6) {
  assign_dt[[paste0("expr_k", k)]] <- cc[[k]]$consensusClass
}
fwrite(assign_dt, file.path(RESULTS, "expression_cluster_assignments.csv"))
saveRDS(cc, file.path(RESULTS, "cc_expression_object.rds"))
saveRDS(list(logmat=logmat, top_genes=top_genes, gene_id2name=setNames(gene_names, gene_ids)),
        file.path(RESULTS, "expression_processed.rds"))

cat("\nPAC values saved. Best K =", best_k, "\n")
fwrite(data.table(K=2:8, PAC=pac), file.path(RESULTS, "expression_PAC.csv"))
cat("=== Expression clustering complete ===\n")
