# =============================================================
# 07_mofa_integration.R
# 多组学整合分型 (MOFA2) + 综合分子亚型定义
# =============================================================
source("config.R")
suppressMessages({
  library(data.table); library(reticulate)
  # Use a Python interpreter that has mofapy2 installed (avoids basilisk download).
  # Set env NSCLC_PYTHON to your python.exe/python path; otherwise reticulate autodetects.
  .py <- Sys.getenv("NSCLC_PYTHON", unset = "")
  if (nzchar(.py)) use_python(.py, required = TRUE)
  library(MOFA2); library(matrixStats)
})
set.seed(20260707)

clin <- fread(file.path(DATA, "alchemist_master_clinical.csv"))

# ---- 载入各组学 ----
# 表达 (处理后)
ep <- readRDS(file.path(RESULTS, "expression_processed.rds"))
logexpr <- ep$logmat
id2name <- ep$gene_id2name

# CNA
cp <- readRDS(file.path(RESULTS, "cna_processed.rds"))
cna_call <- cp$call_mat  # gene_id rows

# 突变 (binary)
mut <- fread(file.path(MATDIR, "mutation_binary.tsv.gz"))
mut_genes <- mut$gene
mutmat <- as.matrix(mut[, -"gene"]); rownames(mutmat) <- mut_genes

# miRNA
mir <- fread(file.path(MATDIR, "mirna_rpm.tsv.gz"))
mir_ids <- mir$miRNA_ID
mirmat <- as.matrix(mir[, -"miRNA_ID"]); rownames(mirmat) <- mir_ids
mirmat <- log2(mirmat + 1)

# ---- 全组学交集样本 (core4: expr+mut+cna+mirna) ----
common <- Reduce(intersect, list(colnames(logexpr), colnames(cna_call),
                                 colnames(mutmat), colnames(mirmat)))
cat("Common samples (all 4 omics):", length(common), "\n")

# ---- 特征选择 (降维, MOFA适宜规模) ----
# 表达: top 3000 高变
ev <- rowVars(logexpr[, common]); names(ev) <- rownames(logexpr)
expr_feat <- names(sort(ev, decreasing=TRUE))[1:3000]
expr_v <- logexpr[expr_feat, common]
# 中心化
expr_v <- t(scale(t(expr_v), center=TRUE, scale=FALSE))

# CNA: top 3000 高变(去NA), 用绝对拷贝数中心化
cna_num <- cp$cn_mat[, common]
cv <- rowVars(cna_num, na.rm=TRUE); names(cv) <- rownames(cna_num)
cv <- cv[!is.na(cv)]
cna_feat <- names(sort(cv, decreasing=TRUE))[1:3000]
cna_v <- cna_num[cna_feat, common]
cna_v <- t(scale(t(cna_v), center=TRUE, scale=FALSE))

# 突变: 频率>=3%的基因
mut_v <- mutmat[, common]
mut_freq <- rowMeans(mut_v)
mut_v <- mut_v[mut_freq >= 0.03, ]
cat("Mutation features (>=3%):", nrow(mut_v), "\n")

# miRNA: top 300 高变
mv <- rowVars(mirmat[, common]); names(mv) <- rownames(mirmat)
mir_feat <- names(sort(mv, decreasing=TRUE))[1:300]
mir_v <- mirmat[mir_feat, common]
mir_v <- t(scale(t(mir_v), center=TRUE, scale=FALSE))

# ---- 构建MOFA对象 ----
data_list <- list(
  Expression = as.matrix(expr_v),
  CNA        = as.matrix(cna_v),
  Mutation   = as.matrix(mut_v),
  miRNA      = as.matrix(mir_v)
)
cat("View dimensions:\n")
for (v in names(data_list)) cat(sprintf("  %s: %d x %d\n", v, nrow(data_list[[v]]), ncol(data_list[[v]])))

MOFAobject <- create_mofa(data_list)

# 数据选项: 突变为二元(bernoulli)
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
model_opts$likelihoods["Mutation"] <- "bernoulli"
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 20260707
train_opts$maxiter <- 1000

MOFAobject <- prepare_mofa(MOFAobject, data_options=data_opts,
                           model_options=model_opts, training_options=train_opts)

outfile <- file.path(RESULTS, "mofa_model.hdf5")
MOFAobject <- run_mofa(MOFAobject, outfile=outfile, use_basilisk=FALSE)
saveRDS(MOFAobject, file.path(RESULTS, "mofa_object.rds"))

# ---- 方差分解 ----
vexp <- get_variance_explained(MOFAobject)
saveRDS(vexp, file.path(RESULTS, "mofa_variance.rds"))
cat("\n=== MOFA variance explained (per factor/view) ===\n")
print(round(vexp$r2_per_factor[[1]], 2))

# ---- 因子聚类 → 综合亚型 ----
factors <- get_factors(MOFAobject, factors="all")[[1]]
fwrite(data.table(patient_id=rownames(factors), factors),
       file.path(RESULTS, "mofa_factors.csv"))

# 用MOFA因子做共识聚类
suppressMessages(library(ConsensusClusterPlus))
fm <- t(factors)  # factors x samples
fm <- fm[complete.cases(fm), , drop=FALSE]
cc_dir <- file.path(RESULTS, "cc_mofa"); dir.create(cc_dir, showWarnings=FALSE)
ccf <- ConsensusClusterPlus(fm, maxK=8, reps=1000, pItem=0.9, pFeature=1,
                            clusterAlg="km", distance="euclidean",
                            seed=20260707, plot="png", title=cc_dir, writeTable=FALSE)
calc_PAC <- function(cc,k,x1=0.1,x2=0.9){M<-cc[[k]]$consensusMatrix;v<-M[lower.tri(M)];Fn<-ecdf(v);Fn(x2)-Fn(x1)}
pac <- sapply(2:8, calc_PAC, cc=ccf); names(pac)<-2:8
cat("\nMOFA-factor cluster PAC:\n"); print(round(pac,4))

assign_dt <- data.table(patient_id=colnames(fm))
for (k in 2:6) assign_dt[[paste0("mofa_k",k)]] <- ccf[[k]]$consensusClass
fwrite(assign_dt, file.path(RESULTS, "mofa_cluster_assignments.csv"))
fwrite(data.table(K=2:8, PAC=pac), file.path(RESULTS, "mofa_PAC.csv"))
saveRDS(ccf, file.path(RESULTS, "cc_mofa_object.rds"))

cat("\n=== MOFA integration complete ===\n")
