# =============================================================
# V2_01_tcga_validation.R
# 短板1补救: TCGA(含分期+生存)验证IMS分型
#   (a) 迁移IMS分型至TCGA (基因组质心)
#   (b) 分期在亚型间的分布 (证明分型≠分期替身)
#   (c) 校正分期的多变量Cox (证明分型独立于分期的预后价值)
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(survival); library(survminer) })
set.seed(20260709)

# ---- 载入ALCHEMIST基因组质心(V1已建) ----
centroid <- fread(file.path(RESULTS, "ims_genomic_centroids.csv"))
subtypes <- c("IMS1","IMS2","IMS3","IMS4")
shared_drivers <- gsub("^mut_","", grep("^mut_", centroid$feature, value=TRUE))
cat("Shared drivers:", length(shared_drivers), "\n")

# ---- 载入TCGA分析数据集 ----
tcga <- fread(file.path(TCGADIR, "tcga_analysis_dataset.csv"))
cat("TCGA patients:", nrow(tcga), "\n")

# 构建TCGA同构特征 (mut_* + TMB + FGA)
# TMB: TCGA用TMB_NONSYNONYMOUS (mut/Mb, 与ALCHEMIST tmb_per_mb同尺度)
# FGA: TCGA FRACTION_GENOME_ALTERED (与ALCHEMIST FGA同定义)
feat_mut <- paste0("mut_", shared_drivers)
tcga_feat <- tcga[, c("patient", feat_mut), with=FALSE]
tcga_feat[, TMB := tcga$TMB_NONSYNONYMOUS]
tcga_feat[, FGA := tcga$FRACTION_GENOME_ALTERED]

# 质心矩阵
cent_m <- as.matrix(centroid[, ..subtypes]); rownames(cent_m) <- centroid$feature
# 特征顺序对齐质心
cent_features <- c(feat_mut, "TMB", "FGA")
cent_aligned <- cent_m[cent_features, ]

# TCGA特征矩阵
X <- as.matrix(tcga_feat[, c(feat_mut,"TMB","FGA"), with=FALSE])
rownames(X) <- tcga_feat$patient
# 缺失填充: 突变NA->0, TMB/FGA->中位
X[, feat_mut][is.na(X[, feat_mut])] <- 0
for (cc in c("TMB","FGA")) X[is.na(X[,cc]),cc] <- median(X[,cc],na.rm=TRUE)

# z标准化(按特征)
Xs <- scale(X)
cent_s <- scale(t(cent_aligned)); cent_s <- t(cent_s)  # features x subtypes

# 最近质心分配(相关性最高)
assign_tcga <- apply(Xs, 1, function(x){
  cors <- sapply(subtypes, function(s) cor(x, cent_s[,s], use="complete.obs"))
  subtypes[which.max(cors)]
})
tcga[, IMS := assign_tcga[patient]]
cat("\n=== TCGA IMS assignment ===\n"); print(table(tcga$IMS))

# ---- (b) 分期 × 亚型 分布 ----
cat("\n=== 分期 × 亚型 (关键: 证明分型不是分期替身) ===\n")
stage_tab <- table(tcga$IMS, tcga$stage_simple)
print(stage_tab)
print(round(prop.table(stage_tab, 1)*100, 1))
chisq_stage <- chisq.test(stage_tab)
cat(sprintf("Chi-square p=%.3f\n", chisq_stage$p.value))
# 早期(I-II)比例 by 亚型
tcga[, early := stage_simple %in% c("I","II")]
cat("\n各亚型早期(I-II)比例:\n")
print(tcga[, .(early_pct=round(100*mean(early,na.rm=TRUE),1), n=.N), by=IMS][order(IMS)])

# ---- (c) 生存分析: 全部 + 仅可切除(I-III) ----
tcga[, os_event := as.integer(grepl("DECEASED", OS_STATUS))]
tcga[, OS_MONTHS := as.numeric(OS_MONTHS)]
surv <- tcga[!is.na(OS_MONTHS) & OS_MONTHS>=0 & !is.na(IMS)]
cat(sprintf("\n=== TCGA OS cohort: n=%d, events=%d ===\n", nrow(surv), sum(surv$os_event)))

fit <- survfit(Surv(OS_MONTHS, os_event) ~ IMS, data=surv)
sd <- survdiff(Surv(OS_MONTHS, os_event) ~ IMS, data=surv)
cat(sprintf("OS by IMS log-rank p=%.2e\n", 1-pchisq(sd$chisq, length(sd$n)-1)))
print(summary(fit)$table[, c("records","events","median")])

# 单变量Cox
surv[, IMS := relevel(factor(IMS), ref="IMS1")]
cox_uni <- coxph(Surv(OS_MONTHS, os_event) ~ IMS, data=surv)
cat("\n=== 单变量Cox (IMS only) ===\n"); print(summary(cox_uni)$conf.int)

# 多变量Cox: 校正分期+组织学+年龄
surv[, stage_f := factor(stage_simple, levels=c("I","II","III","IV"))]
cox_multi <- coxph(Surv(OS_MONTHS, os_event) ~ IMS + stage_f + histology, data=surv)
cat("\n=== 多变量Cox (校正分期+组织学) ===\n")
print(summary(cox_multi))
cat(sprintf("\nC-index: %.3f\n", summary(cox_multi)$concordance[1]))

# ---- DFS分析(可切除人群更相关) ----
tcga[, DFS_MONTHS := as.numeric(DFS_MONTHS)]
tcga[, dfs_event := as.integer(grepl("Recur|Progres|1:", DFS_STATUS))]
dfs <- tcga[!is.na(DFS_MONTHS) & DFS_MONTHS>=0 & !is.na(IMS)]
if (nrow(dfs) > 50) {
  dfs[, IMS := relevel(factor(IMS), ref="IMS1")]
  dfs[, stage_f := factor(stage_simple, levels=c("I","II","III","IV"))]
  cat(sprintf("\n=== TCGA DFS cohort: n=%d, events=%d ===\n", nrow(dfs), sum(dfs$dfs_event)))
  sd2 <- survdiff(Surv(DFS_MONTHS, dfs_event) ~ IMS, data=dfs)
  cat(sprintf("DFS by IMS log-rank p=%.2e\n", 1-pchisq(sd2$chisq, length(sd2$n)-1)))
  cox_dfs <- coxph(Surv(DFS_MONTHS, dfs_event) ~ IMS + stage_f + histology, data=dfs)
  cat("=== DFS多变量Cox(校正分期) ===\n"); print(summary(cox_dfs)$conf.int)
}

# ---- 保存 ----
fwrite(tcga, file.path(RESULTS, "tcga_with_ims.csv"))
saveRDS(list(fit=fit, cox_uni=cox_uni, cox_multi=cox_multi, surv=surv,
             stage_tab=stage_tab, chisq_stage=chisq_stage),
        file.path(RESULTS, "tcga_validation_objects.rds"))
cat("\n=== TCGA validation complete ===\n")
