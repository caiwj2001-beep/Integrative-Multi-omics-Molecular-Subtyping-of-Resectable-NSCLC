# =============================================================
# 06_subtype_definition.R
# 整合确定最终分子亚型 + 亚型生物学特征刻画 + 统计检验
# 运行前提: 02(表达聚类) + 03(CNA) + 04(免疫) + 05(MOFA) 已完成
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(matrixStats) })

clin <- fread(file.path(DATA, "alchemist_master_clinical.csv"))

# ---- 载入MOFA综合聚类 (核心分型) ----
mofa_assign <- fread(file.path(RESULTS, "mofa_cluster_assignments.csv"))
mofa_pac <- fread(file.path(RESULTS, "mofa_PAC.csv"))
cat("=== MOFA PAC ===\n"); print(mofa_pac)

# 选择K: 优先PAC低且亚型数在3-5(可解释) ; 默认取K=4, 可根据PAC调整
best_k <- 4
kcol <- paste0("mofa_k", best_k)
cat("\nUsing integrated K =", best_k, "\n")

subtypes <- mofa_assign[, .(patient_id, cluster = get(kcol))]
subtypes[, subtype := paste0("IMS", cluster)]  # Integrated Molecular Subtype
cat("\n亚型分布:\n"); print(table(subtypes$subtype))

# ---- 合并所有分子特征 ----
D <- merge(clin, subtypes[, .(patient_id, subtype)], by="patient_id")

# TMB
tmb <- fread(file.path(MATDIR, "tmb.csv"))
D <- merge(D, tmb[, .(patient_id, tmb_per_mb, nonsyn_count)], by="patient_id", all.x=TRUE)
# FGA
fga <- fread(file.path(RESULTS, "cna_genomic_instability.csv"))
D <- merge(D, fga, by="patient_id", all.x=TRUE)
# ESTIMATE
est <- fread(file.path(RESULTS, "estimate_scores.csv"))
D <- merge(D, est, by="patient_id", all.x=TRUE)

fwrite(D, file.path(RESULTS, "master_with_subtypes.csv"))
cat("\nMaster table with subtypes:", nrow(D), "patients,", sum(!is.na(D$subtype)), "with subtype\n")

# ---- 亚型间差异检验 ----
sink(file.path(RESULTS, "subtype_characterization.txt"))
cat("=== Integrated Molecular Subtypes (IMS) Characterization ===\n\n")
Dsub <- D[!is.na(subtype)]
cat("N per subtype:\n"); print(table(Dsub$subtype)); cat("\n")

# 组织学分布 (卡方)
cat("--- Histology x Subtype ---\n")
ht <- table(Dsub$subtype, Dsub$histology_group); print(ht)
print(chisq.test(table(Dsub$subtype, Dsub$histology_group=="LUSC"))); cat("\n")

# 连续变量: TMB, FGA, ImmuneScore (Kruskal-Wallis)
for (v in c("tmb_per_mb","FGA","ImmuneScore","StromalScore","ESTIMATEScore","age_years")) {
  if (v %in% colnames(Dsub)) {
    kw <- kruskal.test(Dsub[[v]] ~ Dsub$subtype)
    med <- Dsub[, .(median=round(median(get(v),na.rm=TRUE),3)), by=subtype][order(subtype)]
    cat(sprintf("--- %s (Kruskal-Wallis p=%.2e) ---\n", v, kw$p.value))
    print(med); cat("\n")
  }
}

# 治疗臂分布
cat("--- Treatment arm x Subtype ---\n")
print(table(Dsub$subtype, Dsub$treatment_arm))
print(chisq.test(table(Dsub$subtype, Dsub$treatment_arm))); cat("\n")

# 性别
cat("--- Sex x Subtype ---\n")
print(table(Dsub$subtype, Dsub$sex)); cat("\n")
sink()

cat("Characterization written to subtype_characterization.txt\n")

# ---- 驱动基因在亚型间的频率 ----
mut <- fread(file.path(MATDIR, "mutation_binary.tsv.gz"))
mut_genes <- mut$gene
mutmat <- as.matrix(mut[, -"gene"]); rownames(mutmat) <- mut_genes
drivers <- c("TP53","EGFR","KRAS","STK11","KEAP1","SMARCA4","NF1","BRAF","MET",
             "RBM10","PIK3CA","RB1","CDKN2A","ATM","ARID1A","SETD2","PTEN","NFE2L2")
drivers <- intersect(drivers, rownames(mutmat))
res <- list()
for (g in drivers) {
  row <- data.table(subtype="", freq=numeric())
  freqs <- sapply(sort(unique(Dsub$subtype)), function(s){
    pts <- Dsub[subtype==s, patient_id]
    pts <- intersect(pts, colnames(mutmat))
    if (length(pts)==0) return(NA)
    round(100*mean(mutmat[g, pts]),1)
  })
  # Fisher检验: 该基因突变是否在亚型间不同
  pts_all <- intersect(Dsub$patient_id, colnames(mutmat))
  mstat <- mutmat[g, pts_all]
  styp <- Dsub[match(pts_all, patient_id), subtype]
  ptest <- tryCatch(fisher.test(table(styp, mstat), simulate.p.value=TRUE, B=10000)$p.value,
                    error=function(e) NA)
  res[[g]] <- c(freqs, p=ptest)
}
driver_bysubtype <- as.data.table(do.call(rbind, res), keep.rownames="gene")
fwrite(driver_bysubtype, file.path(RESULTS, "driver_freq_by_subtype.csv"))
cat("\n=== Driver frequency by subtype ===\n"); print(driver_bysubtype)

cat("\n=== Subtype definition & characterization complete ===\n")
