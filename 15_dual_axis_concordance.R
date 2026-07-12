# =============================================================
# V3_01_dual_axis_concordance.R
# 核心: 在同一批TCGA患者上, 分别用"基因组轴"和"转录组轴"分型,
#       计算Cohen's kappa + 混淆矩阵, 实证支撑"不同轴捕捉不同维度"
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(matrixStats) })
set.seed(20260709)
subtypes <- c("IMS1","IMS2","IMS3","IMS4")

# ============ 方法A: 基因组轴分型 (与MSK迁移同方法) ============
centroid_g <- fread(file.path(RESULTS, "ims_genomic_centroids.csv"))
shared_drivers <- gsub("^mut_","", grep("^mut_", centroid_g$feature, value=TRUE))
feat_mut <- paste0("mut_", shared_drivers)
cent_features <- c(feat_mut, "TMB", "FGA")
cent_g <- as.matrix(centroid_g[, ..subtypes]); rownames(cent_g) <- centroid_g$feature
cent_g <- cent_g[cent_features, ]

tcga <- fread(file.path(TCGADIR, "tcga_analysis_dataset.csv"))
Xg <- as.matrix(tcga[, ..feat_mut])
Xg <- cbind(Xg, TMB=as.numeric(tcga$TMB_NONSYNONYMOUS), FGA=as.numeric(tcga$FRACTION_GENOME_ALTERED))
rownames(Xg) <- tcga$patient
Xg[, feat_mut][is.na(Xg[, feat_mut])] <- 0
for (cc in c("TMB","FGA")) Xg[is.na(Xg[,cc]),cc] <- median(Xg[,cc],na.rm=TRUE)
Xgs <- scale(Xg)
cent_gs <- scale(t(cent_g)); cent_gs <- t(cent_gs)
assign_genomic <- apply(Xgs, 1, function(x){
  cors <- sapply(subtypes, function(s) cor(x, cent_gs[,s], use="complete.obs"))
  subtypes[which.max(cors)]
})
cat("基因组轴分型:\n"); print(table(assign_genomic))

# ============ 方法B: 转录组轴分型 (与V2 TCGA迁移同方法) ============
sig <- readRDS(file.path(RESULTS, "expression_signature.rds"))
centroid_e <- sig$centroid_expr; all_markers <- sig$all_markers
te <- fread(file.path(TCGADIR, "tcga_signature_expr_matrix.csv"))
tcga_e <- as.matrix(te[,-1]); rownames(tcga_e) <- te$gene
tcga_e[tcga_e<0] <- 0; tcga_log <- log2(tcga_e+1)
shared_markers <- intersect(all_markers, rownames(tcga_log))
z_alch <- t(scale(t(centroid_e[shared_markers,])))
z_tcga <- t(scale(t(tcga_log[shared_markers,])))
assign_expr <- apply(z_tcga, 2, function(x){
  if (all(is.na(x))) return(NA)
  cors <- sapply(subtypes, function(s) cor(x, z_alch[,s], method="spearman", use="complete.obs"))
  subtypes[which.max(cors)]
})
cat("\n转录组轴分型:\n"); print(table(assign_expr))

# ============ 一致性: 共同患者 ============
dg <- data.table(patient=names(assign_genomic), genomic=assign_genomic)
de <- data.table(patient=names(assign_expr), expr=assign_expr)
D <- merge(dg, de, by="patient")
D <- D[!is.na(genomic) & !is.na(expr)]
cat(sprintf("\n共同患者(双轴均分型): %d\n", nrow(D)))

# 混淆矩阵
cm <- table(Genomic=factor(D$genomic,levels=subtypes), Transcriptomic=factor(D$expr,levels=subtypes))
cat("\n=== 混淆矩阵 (行=基因组轴, 列=转录组轴) ===\n"); print(cm)
cat("\n行归一化(%):\n"); print(round(prop.table(cm,1)*100,1))

# Cohen's kappa
if (!requireNamespace("irr", quietly=TRUE)) install.packages("irr", repos="https://cloud.r-project.org")
suppressMessages(library(irr))
kap <- kappa2(D[,.(genomic, expr)], weight="unweighted")
cat(sprintf("\n=== Cohen's kappa = %.3f (z=%.2f, p=%.2e) ===\n", kap$value, kap$statistic, kap$p.value))
# 简单一致率
agree_pct <- 100*mean(D$genomic == D$expr)
cat(sprintf("原始一致率 (identical label): %.1f%%\n", agree_pct))

# kappa解读
interp <- if(kap$value<0.20)"slight" else if(kap$value<0.40)"fair" else if(kap$value<0.60)"moderate" else if(kap$value<0.80)"substantial" else "almost perfect"
cat(sprintf("Kappa强度: %s agreement\n", interp))

# ============ 关键: IMS4在两轴的成员是否重叠 ============
cat("\n=== IMS4 跨轴成员分析 (回答'为何IMS4在两队列表现相反') ===\n")
g4 <- D[genomic=="IMS4", patient]; e4 <- D[expr=="IMS4", patient]
overlap4 <- length(intersect(g4,e4))
cat(sprintf("基因组轴IMS4: n=%d; 转录组轴IMS4: n=%d; 重叠: n=%d (Jaccard=%.2f)\n",
            length(g4), length(e4), overlap4, overlap4/length(union(g4,e4))))

# 保存
fwrite(D, file.path(RESULTS, "tcga_dual_axis_assignments.csv"))
saveRDS(list(cm=cm, kappa=kap, agree_pct=agree_pct, D=D), file.path(RESULTS, "dual_axis_concordance.rds"))
cat("\n=== 双轴一致性分析完成 ===\n")
