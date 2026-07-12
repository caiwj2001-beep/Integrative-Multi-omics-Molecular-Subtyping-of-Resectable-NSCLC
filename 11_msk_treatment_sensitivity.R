# =============================================================
# V2_05_msk_treatment_sensitivity.R
# 短板2补救: MSK无治疗信息 → 治疗混杂的敏感性分析
# 逻辑: 若IMS预后差异在(a)测序年代/panel (b)样本类型 各层一致,
#       则预后信号不太可能纯由治疗差异驱动
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(survival) })

msk <- fread(file.path(RESULTS, "msk_subtype_assignments.csv"))
surv <- fread(file.path(RESULTS, "msk_survival_data.csv"))
clin <- fread(file.path(DATA, "msk", "msk_nsclc_clinical.csv"))

# 合并panel(年代代理)
D <- merge(surv, clin[,.(SAMPLE_ID, GENE_PANEL)], by="SAMPLE_ID", all.x=TRUE)
D <- D[!is.na(OS_months) & OS_months>=0 & !is.na(IMS)]
D[, IMS := relevel(factor(IMS), ref="IMS1")]

# 年代分组: IMPACT341/410 = 早期(pre/early-IO), IMPACT468/505 = 近期(IO时代)
D[, era := ifelse(GENE_PANEL %in% c("IMPACT341","IMPACT410"), "Early (2014-2016)", "Recent (2017+)")]
cat("=== Era (panel proxy) distribution ===\n"); print(table(D$era, D$IMS))

sink(file.path(RESULTS, "msk_sensitivity_analysis.txt"))
cat("=== MSK-IMPACT 治疗混杂敏感性分析 ===\n\n")
cat("原理: MSK队列无治疗注释。若IMS的预后分层在不同测序年代(治疗时代代理)\n")
cat("和样本类型中一致, 则该预后信号不太可能是治疗差异的假象。\n\n")

# 1. 全队列基线
cox_all <- coxph(Surv(OS_months, OS_event) ~ IMS + histology_group + SAMPLE_TYPE, data=D)
cat("--- (1) 全队列多变量Cox (基线) ---\n")
print(summary(cox_all)$conf.int); cat("C-index:", round(summary(cox_all)$concordance[1],3),"\n\n")

# 2. 按年代分层
for (e in c("Early (2014-2016)", "Recent (2017+)")) {
  sub <- D[era==e]
  if (nrow(sub) > 100 && length(unique(sub$IMS))>=3) {
    sub[, IMS := relevel(factor(IMS), ref="IMS1")]
    cox_e <- coxph(Surv(OS_months, OS_event) ~ IMS + histology_group, data=sub)
    cat(sprintf("--- (2) %s 层 (n=%d, events=%d) ---\n", e, nrow(sub), sum(sub$OS_event)))
    ci <- summary(cox_e)$conf.int
    print(round(ci[grep("IMS", rownames(ci)), c("exp(coef)","lower .95","upper .95")],2))
    cat("\n")
  }
}

# 3. 按样本类型分层(原发vs转移)
for (st in c("Primary","Metastasis")) {
  sub <- D[SAMPLE_TYPE==st]
  if (nrow(sub) > 100) {
    sub[, IMS := relevel(factor(IMS), ref="IMS1")]
    cox_s <- coxph(Surv(OS_months, OS_event) ~ IMS + histology_group, data=sub)
    cat(sprintf("--- (3) %s 样本层 (n=%d, events=%d) ---\n", st, nrow(sub), sum(sub$OS_event)))
    ci <- summary(cox_s)$conf.int
    print(round(ci[grep("IMS", rownames(ci)), c("exp(coef)","lower .95","upper .95")],2))
    cat("\n")
  }
}

# 4. IMS1(最佳预后)在各层的HR一致性总结
cat("--- (4) 结论 ---\n")
cat("若IMS1的生存优势(其他亚型HR>1)在早期与近期测序年代、原发与转移样本中\n")
cat("方向一致, 则支持该分层反映肿瘤内在生物学, 而非治疗时代或治疗选择的混杂。\n")
sink()

# 交互检验: IMS × era
D[, IMS_bin := ifelse(IMS=="IMS1","IMS1","other")]
cox_int <- coxph(Surv(OS_months,OS_event) ~ IMS_bin*era + histology_group, data=D)
cat("\n=== IMS×era交互检验 ===\n")
print(summary(cox_int)$coefficients)
int_p <- summary(cox_int)$coefficients[grep(":", rownames(summary(cox_int)$coefficients)), "Pr(>|z|)"]
cat(sprintf("\n交互项p值: %.3f (>0.05表示各年代效应一致,无显著治疗时代混杂)\n", int_p[1]))

fwrite(D, file.path(RESULTS, "msk_sensitivity_data.csv"))
cat("\n=== MSK sensitivity analysis complete ===\n")
