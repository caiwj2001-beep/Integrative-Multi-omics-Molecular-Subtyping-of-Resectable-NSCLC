# =============================================================
# V2_03_tcga_transcriptomic_transfer.R
# 用表达signature将IMS分型迁移到TCGA (转录组, 比基因组轴更忠实)
# 然后: 分期分布 + 校正分期的Cox
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(survival); library(matrixStats) })
set.seed(20260709)

# ---- 载入ALCHEMIST表达signature质心 ----
sig <- readRDS(file.path(RESULTS, "expression_signature.rds"))
centroid_expr <- sig$centroid_expr  # genes x subtypes (ALCHEMIST log2 TPM均值)
all_markers <- sig$all_markers

# ---- 载入TCGA表达(Python已导出CSV) ----
te <- fread(file.path(TCGADIR, "tcga_signature_expr_matrix.csv"))
tcga_e <- as.matrix(te[,-1]); rownames(tcga_e) <- te$gene
cat("TCGA expression matrix:", nrow(tcga_e), "genes x", ncol(tcga_e), "samples\n")

# RSEM -> log2
tcga_e[tcga_e < 0] <- 0
tcga_log <- log2(tcga_e + 1)

# ---- 共同marker基因 ----
shared_markers <- intersect(all_markers, rownames(tcga_log))
cat("Shared markers:", length(shared_markers), "\n")

# 两个队列各自基因内z-score(消除平台差异)
z_alch <- t(scale(t(centroid_expr[shared_markers,])))  # ALCHEMIST质心z
z_tcga <- t(scale(t(tcga_log[shared_markers,])))        # TCGA样本z

# ---- 最近质心(Spearman相关)分配 ----
subtypes <- colnames(centroid_expr)
assign_tcga <- apply(z_tcga, 2, function(x){
  if (all(is.na(x))) return(NA)
  cors <- sapply(subtypes, function(s) cor(x, z_alch[,s], method="spearman", use="complete.obs"))
  subtypes[which.max(cors)]
})
tcga_ims <- data.table(patient=colnames(tcga_log), IMS_expr=assign_tcga)
cat("\n=== TCGA transcriptomic IMS assignment ===\n"); print(table(tcga_ims$IMS_expr))

# ---- 合并临床(分期+生存) ----
tcga <- fread(file.path(TCGADIR, "tcga_analysis_dataset.csv"))
tcga <- merge(tcga, tcga_ims, by="patient", how="left")
tcga <- tcga[!is.na(IMS_expr)]
cat("Merged with clinical:", nrow(tcga), "\n")

# ---- 分期 × 亚型 ----
cat("\n=== 分期 × 转录组IMS ===\n")
st <- table(tcga$IMS_expr, tcga$stage_simple); print(st)
cat("Chi-sq p =", round(chisq.test(st)$p.value,3), "\n")
tcga[, early := stage_simple %in% c("I","II")]
cat("各亚型早期比例:\n"); print(tcga[,.(early_pct=round(100*mean(early,na.rm=T),1),n=.N),by=IMS_expr][order(IMS_expr)])

# ---- 基因组特征验证(亚型是否重现ALCHEMIST生物学) ----
cat("\n=== 亚型基因组特征(验证生物学一致性) ===\n")
tcga[, TMB_NONSYNONYMOUS:=as.numeric(TMB_NONSYNONYMOUS)]
tcga[, FRACTION_GENOME_ALTERED:=as.numeric(FRACTION_GENOME_ALTERED)]
print(tcga[,.(N=.N, TMB=round(median(TMB_NONSYNONYMOUS,na.rm=T),1),
              FGA=round(median(FRACTION_GENOME_ALTERED,na.rm=T),2),
              TP53=round(100*mean(mut_TP53),0), EGFR=round(100*mean(mut_EGFR),0),
              KRAS=round(100*mean(mut_KRAS),0)), by=IMS_expr][order(IMS_expr)])

# ---- OS生存 ----
tcga[, os_event := as.integer(grepl("DECEASED", OS_STATUS))]
tcga[, OS_MONTHS := as.numeric(OS_MONTHS)]
surv <- tcga[!is.na(OS_MONTHS) & OS_MONTHS>=0]
surv[, IMS := relevel(factor(IMS_expr), ref="IMS1")]
surv[, stage_f := factor(stage_simple, levels=c("I","II","III","IV"))]
cat(sprintf("\n=== OS cohort n=%d events=%d ===\n", nrow(surv), sum(surv$os_event)))
sd <- survdiff(Surv(OS_MONTHS,os_event)~IMS, data=surv)
cat(sprintf("Log-rank p=%.3f\n", 1-pchisq(sd$chisq,length(sd$n)-1)))
print(summary(survfit(Surv(OS_MONTHS,os_event)~IMS,data=surv))$table[,c("records","events","median")])

cat("\n=== 多变量Cox (校正分期+组织学) ===\n")
cox <- coxph(Surv(OS_MONTHS,os_event)~IMS+stage_f+histology, data=surv)
print(summary(cox)$conf.int); cat("C-index:", round(summary(cox)$concordance[1],3),"\n")
print(summary(cox)$coefficients[,c("coef","Pr(>|z|)")])

# ---- DFS ----
tcga[, DFS_MONTHS:=as.numeric(DFS_MONTHS)]
tcga[, dfs_event:=as.integer(grepl("Recur|Progres|1:", DFS_STATUS))]
dfs <- tcga[!is.na(DFS_MONTHS)&DFS_MONTHS>=0]
dfs[, IMS:=relevel(factor(IMS_expr),ref="IMS1")]; dfs[, stage_f:=factor(stage_simple,levels=c("I","II","III","IV"))]
if(nrow(dfs)>50){
  cat(sprintf("\n=== DFS cohort n=%d events=%d ===\n", nrow(dfs), sum(dfs$dfs_event)))
  sd2<-survdiff(Surv(DFS_MONTHS,dfs_event)~IMS,data=dfs)
  cat(sprintf("DFS log-rank p=%.3f\n",1-pchisq(sd2$chisq,length(sd2$n)-1)))
  coxd<-coxph(Surv(DFS_MONTHS,dfs_event)~IMS+stage_f+histology,data=dfs)
  cat("DFS多变量Cox(校正分期):\n"); print(summary(coxd)$conf.int)
}

fwrite(tcga, file.path(RESULTS,"tcga_transcriptomic_ims.csv"))
saveRDS(list(surv=surv,cox=cox,stage_tab=st), file.path(RESULTS,"tcga_transcriptomic_objects.rds"))
cat("\n=== TCGA transcriptomic transfer complete ===\n")
