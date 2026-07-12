# =============================================================
# 08_msk_validation.R
# MSK-IMPACT外部验证: 基因组特征迁移分型 + 真实OS生存分析
# 关键: MSK无RNA-seq, 验证分型的"基因组轴"可重复性与预后价值
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(survival); library(survminer) })
set.seed(20260707)

MSKDIR <- file.path(DATA, "msk")

# ===== 1. ALCHEMIST基因组特征定义IMS的"基因组签名" =====
# IMS由: TMB, FGA(基因组不稳定), 关键驱动突变, 焦点CNA 定义
# 提取ALCHEMIST core3各亚型的基因组特征profile
alch <- fread(file.path(RESULTS, "core3_subtypes_final.csv"))
cat("ALCHEMIST core3 subtyped:", nrow(alch), "\n")

# ALCHEMIST突变(binary)
amut <- fread(file.path(MATDIR, "mutation_binary.tsv.gz"))
amut_genes <- amut$gene; amut_m <- as.matrix(amut[,-"gene"]); rownames(amut_m)<-amut_genes

shared_drivers <- c("TP53","EGFR","KRAS","STK11","KEAP1","SMARCA4","NF1","BRAF",
                    "PIK3CA","RB1","CDKN2A","ATM","ARID1A","PTEN","NFE2L2","MET","RBM10")
shared_drivers <- intersect(shared_drivers, rownames(amut_m))

# 各亚型驱动突变频率 + 基因组特征均值 (作为质心)
subtypes <- c("IMS1","IMS2","IMS3","IMS4")
centroid <- data.table(feature=c(paste0("mut_",shared_drivers),"TMB","FGA"))
for (s in subtypes) {
  pts <- intersect(alch[subtype==s, patient_id], colnames(amut_m))
  mut_freq <- rowMeans(amut_m[shared_drivers, pts])
  tmb_m <- mean(alch[subtype==s, tmb_per_mb], na.rm=TRUE)
  fga_m <- mean(alch[subtype==s, FGA], na.rm=TRUE)
  centroid[[s]] <- c(mut_freq, tmb_m, fga_m)
}
fwrite(centroid, file.path(RESULTS, "ims_genomic_centroids.csv"))
cat("\n=== IMS genomic centroids ===\n"); print(centroid)

# ===== 2. MSK-IMPACT: 构建同样的基因组特征 =====
msk_clin <- fread(file.path(MSKDIR, "msk_nsclc_clinical.csv"))
msk_mut <- fread(file.path(MSKDIR, "msk_mutation_binary.csv"))
mg <- msk_mut[[1]]; msk_mut_m <- as.matrix(msk_mut[,-1]); rownames(msk_mut_m)<-mg

# MSK FGA: 从gene CNA计算(非0比例)
msk_cna <- fread(file.path(MSKDIR, "msk_gene_cna.csv"))
cna_genes <- msk_cna$Hugo_Symbol; cna_m <- as.matrix(msk_cna[,-"Hugo_Symbol"])
rownames(cna_m) <- cna_genes
msk_fga <- colMeans(abs(cna_m)>=1, na.rm=TRUE)  # panel-based FGA proxy

# 组装MSK特征表
msk_samples <- intersect(msk_clin$SAMPLE_ID, colnames(msk_mut_m))
msk_feat <- data.table(SAMPLE_ID=msk_samples)
for (g in shared_drivers) {
  msk_feat[[paste0("mut_",g)]] <- if(g %in% rownames(msk_mut_m)) msk_mut_m[g, msk_samples] else 0L
}
msk_feat[, TMB := msk_clin[match(msk_samples, SAMPLE_ID), TMB_SCORE]]
msk_feat[, FGA := msk_fga[msk_samples]]
cat("\nMSK feature matrix:", nrow(msk_feat), "x", ncol(msk_feat), "\n")

# ===== 3. 最近质心分类 (基于基因组特征) =====
# 标准化特征后计算到各质心的相关性
cent_m <- as.matrix(centroid[,..subtypes]); rownames(cent_m)<-centroid$feature
feat_cols <- c(paste0("mut_",shared_drivers),"TMB","FGA")
msk_X <- as.matrix(msk_feat[, ..feat_cols]); rownames(msk_X)<-msk_feat$SAMPLE_ID
# z-score(用ALCHEMIST尺度近似: 这里各特征标准化)
msk_Xs <- scale(msk_X)
cent_s <- scale(t(cent_m))  # subtypes x features
cent_s <- t(cent_s)  # features x subtypes
# 分配到相关性最高的亚型
assign_msk <- apply(msk_Xs, 1, function(x){
  cors <- sapply(subtypes, function(s) cor(x, cent_s[,s], use="complete.obs"))
  subtypes[which.max(cors)]
})
msk_feat[, IMS := assign_msk]
cat("\n=== MSK subtype assignment ===\n"); print(table(msk_feat$IMS))

# ===== 4. 真实OS生存分析 =====
surv_dt <- merge(msk_feat[,.(SAMPLE_ID, IMS)],
                 msk_clin[,.(SAMPLE_ID, OS_months, OS_event, SAMPLE_TYPE, histology_group)],
                 by="SAMPLE_ID")
surv_dt <- surv_dt[!is.na(OS_months) & !is.na(OS_event) & OS_months>=0]
cat("\nSurvival cohort:", nrow(surv_dt), " events:", sum(surv_dt$OS_event), "\n")
cat("By subtype:\n"); print(table(surv_dt$IMS))

# 全队列KM
fit <- survfit(Surv(OS_months, OS_event) ~ IMS, data=surv_dt)
sd <- survdiff(Surv(OS_months, OS_event) ~ IMS, data=surv_dt)
p_all <- 1 - pchisq(sd$chisq, length(sd$n)-1)
cat(sprintf("\n=== OS by IMS (all NSCLC, n=%d): log-rank p=%.2e ===\n", nrow(surv_dt), p_all))
print(summary(fit)$table[, c("records","events","median")])

# 仅原发灶
surv_prim <- surv_dt[SAMPLE_TYPE=="Primary"]
if (nrow(surv_prim)>50) {
  sdp <- survdiff(Surv(OS_months, OS_event) ~ IMS, data=surv_prim)
  pp <- 1-pchisq(sdp$chisq, length(sdp$n)-1)
  cat(sprintf("\n=== OS by IMS (Primary only, n=%d): log-rank p=%.2e ===\n", nrow(surv_prim), pp))
  print(summary(survfit(Surv(OS_months,OS_event)~IMS,data=surv_prim))$table[,c("records","events","median")])
}

# 多变量Cox (校正组织学, 样本类型)
surv_dt[, IMS := relevel(factor(IMS), ref="IMS1")]
cox <- coxph(Surv(OS_months, OS_event) ~ IMS + histology_group + SAMPLE_TYPE, data=surv_dt)
cat("\n=== Multivariable Cox (ref=IMS1) ===\n"); print(summary(cox))

fwrite(surv_dt, file.path(RESULTS, "msk_survival_data.csv"))
fwrite(msk_feat, file.path(RESULTS, "msk_subtype_assignments.csv"))
saveRDS(list(fit=fit, cox=cox, surv_dt=surv_dt), file.path(RESULTS, "msk_survival_objects.rds"))
cat("\n=== MSK validation complete ===\n")
