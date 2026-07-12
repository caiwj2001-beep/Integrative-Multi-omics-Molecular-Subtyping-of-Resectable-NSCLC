# =============================================================
# 01_mutation_landscape.R
# ALCHEMIST突变景观分析 (maftools)
# =============================================================
source("config.R")
suppressMessages({
  library(maftools); library(data.table)
})

# ---- 读取合并MAF ----
maf_file <- file.path(MATDIR, "alchemist_merged.maf.gz")
maf_dt <- fread(maf_file)
cat("Merged MAF rows:", nrow(maf_dt), " patients:", length(unique(maf_dt$patient_id)), "\n")

# 读取临床(作为maftools注释)
clin <- fread(file.path(DATA, "alchemist_master_clinical.csv"))
clin_maf <- clin[has_maf == TRUE, .(
  Tumor_Sample_Barcode = patient_id,
  treatment_arm, histology_group, sex, age_years
)]

# maftools需要的最小列: Hugo_Symbol, Variant_Classification, Variant_Type, Tumor_Sample_Barcode
# 已有Tumor_Sample_Barcode=patient_id
maf_dt[, Tumor_Sample_Barcode := patient_id]

# 过滤: 仅保留有临床注释的样本
maf_dt <- maf_dt[Tumor_Sample_Barcode %in% clin_maf$Tumor_Sample_Barcode]

# 创建maf对象
laml <- read.maf(maf = maf_dt, clinicalData = clin_maf,
                 vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
                               "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site",
                               "Translation_Start_Site","Nonstop_Mutation"),
                 verbose = FALSE)

cat("\n=== MAF Summary ===\n")
print(laml)

# ---- 保存突变摘要 ----
gene_summary <- getGeneSummary(laml)
fwrite(gene_summary, file.path(RESULTS, "mutation_gene_summary.csv"))
sample_summary <- getSampleSummary(laml)
fwrite(sample_summary, file.path(RESULTS, "mutation_sample_summary.csv"))

cat("\n=== Top 20 mutated genes ===\n")
print(gene_summary[1:20, .(Hugo_Symbol, MutatedSamples, total)])

# ---- 已知NSCLC驱动基因 ----
nsclc_drivers <- c("TP53","EGFR","KRAS","STK11","KEAP1","NF1","BRAF","MET","RBM10",
                   "SMARCA4","ARID1A","PIK3CA","RB1","CDKN2A","ATM","PTEN","U2AF1",
                   "SETD2","NFE2L2","PTPRD","MGA","APC","CTNNB1","ERBB2","ALK","ROS1","RET")
drivers_present <- intersect(nsclc_drivers, gene_summary$Hugo_Symbol)
cat("\nNSCLC driver genes present:", length(drivers_present), "\n")

# 驱动基因频率表
driver_freq <- gene_summary[Hugo_Symbol %in% nsclc_drivers,
                            .(Hugo_Symbol, MutatedSamples,
                              Freq_pct = round(100*MutatedSamples/as.numeric(laml@summary$summary[3]), 1))]
setorder(driver_freq, -MutatedSamples)
fwrite(driver_freq, file.path(RESULTS, "driver_gene_frequency.csv"))
cat("\n=== NSCLC Driver Frequencies ===\n"); print(driver_freq)

# ---- 突变互斥/共现分析 ----
tryCatch({
  som_int <- somaticInteractions(maf = laml, top = 30, pvalue = c(0.05, 0.01),
                                 returnAll = TRUE)
  if (!is.null(som_int)) {
    si_dt <- as.data.table(som_int)
    fwrite(si_dt, file.path(RESULTS, "somatic_interactions.csv"))
    cat("\nSomatic interactions saved:", nrow(si_dt), "pairs\n")
  }
}, error = function(e) cat("somaticInteractions error:", conditionMessage(e), "\n"))

# ---- 按组织学分层的驱动基因频率 ----
for (hg in c("LUAD","LUSC")) {
  samples_hg <- clin_maf[histology_group == hg, Tumor_Sample_Barcode]
  maf_hg <- subsetMaf(laml, tsb = samples_hg, mafObj = TRUE)
  gs_hg <- getGeneSummary(maf_hg)
  n_hg <- length(samples_hg)
  dfreq <- gs_hg[Hugo_Symbol %in% nsclc_drivers,
                 .(Hugo_Symbol, MutatedSamples,
                   Freq_pct = round(100*MutatedSamples/n_hg,1))]
  setorder(dfreq, -MutatedSamples)
  fwrite(dfreq, file.path(RESULTS, paste0("driver_freq_", hg, ".csv")))
  cat(sprintf("\n%s (n=%d) top drivers:\n", hg, n_hg))
  print(head(dfreq, 12))
}

# ---- 保存maf对象供后续 ----
saveRDS(laml, file.path(RESULTS, "maf_object.rds"))
cat("\n=== Mutation landscape analysis complete ===\n")
