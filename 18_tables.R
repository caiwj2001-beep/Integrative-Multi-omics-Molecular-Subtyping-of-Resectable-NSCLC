# =============================================================
# 11_tables_docx.R
# 投稿表格 (DOCX格式, officer + flextable)
# Table1: 队列临床特征 by subtype
# Table2: MSK多变量Cox
# STable1: 驱动基因频率 by subtype
# STable2: 亚型分子特征汇总
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(officer); library(flextable); library(survival) })

final <- fread(file.path(RESULTS, "core3_subtypes_final.csv"))[!is.na(subtype)]

# ============ TABLE 1: 临床特征 by subtype ============
make_row_cat <- function(dt, var, label, levels_use=NULL) {
  if (is.null(levels_use)) levels_use <- sort(unique(dt[[var]]))
  rows <- list()
  for (lv in levels_use) {
    r <- c(paste0("  ", lv))
    for (s in c("IMS1","IMS2","IMS3","IMS4")) {
      sub <- dt[subtype==s]
      n <- sum(sub[[var]]==lv, na.rm=TRUE); pct <- 100*n/nrow(sub)
      r <- c(r, sprintf("%d (%.1f)", n, pct))
    }
    # 总体
    n <- sum(dt[[var]]==lv,na.rm=TRUE); r <- c(r, sprintf("%d (%.1f)", n, 100*n/nrow(dt)))
    rows[[lv]] <- r
  }
  # p值(卡方)
  tb <- table(dt$subtype, dt[[var]])
  p <- tryCatch(chisq.test(tb)$p.value, error=function(e) NA)
  list(header=c(label,"","","","",""), rows=rows, p=p)
}

make_row_cont <- function(dt, var, label) {
  r <- c(label)
  for (s in c("IMS1","IMS2","IMS3","IMS4")) {
    v <- dt[subtype==s][[var]]
    r <- c(r, sprintf("%.1f (%.1f-%.1f)", median(v,na.rm=TRUE),
                      quantile(v,0.25,na.rm=TRUE), quantile(v,0.75,na.rm=TRUE)))
  }
  v <- dt[[var]]
  r <- c(r, sprintf("%.1f (%.1f-%.1f)", median(v,na.rm=TRUE),
                    quantile(v,0.25,na.rm=TRUE), quantile(v,0.75,na.rm=TRUE)))
  p <- tryCatch(kruskal.test(dt[[var]]~dt$subtype)$p.value, error=function(e) NA)
  list(row=r, p=p)
}

# 构建表格数据
ns <- final[, .N, by=subtype][order(subtype)]
n_hdr <- sapply(c("IMS1","IMS2","IMS3","IMS4"), function(s) sum(final$subtype==s))
tbl1 <- data.table(
  Characteristic = character(), IMS1=character(), IMS2=character(),
  IMS3=character(), IMS4=character(), Overall=character(), P=character())

add_line <- function(tbl, cells, p="") {
  rbind(tbl, as.data.table(as.list(setNames(c(cells, p),
    c("Characteristic","IMS1","IMS2","IMS3","IMS4","Overall","P")))), fill=TRUE)
}
# N行
tbl1 <- add_line(tbl1, c("No. of patients", as.character(n_hdr), as.character(nrow(final))))
# 年龄
ar <- make_row_cont(final, "age_years", "Age, median (IQR)")
tbl1 <- add_line(tbl1, ar$row, sprintf("%.2g", ar$p))
# 性别
tbl1 <- add_line(tbl1, c("Sex, n (%)","","","","",""))
for (lv in c("female","male")) {
  r <- c(paste0("  ",lv))
  for (s in c("IMS1","IMS2","IMS3","IMS4")) { sub<-final[subtype==s]; n<-sum(sub$sex==lv); r<-c(r,sprintf("%d (%.1f)",n,100*n/nrow(sub))) }
  n<-sum(final$sex==lv); r<-c(r,sprintf("%d (%.1f)",n,100*n/nrow(final)))
  tbl1 <- add_line(tbl1, r)
}
psex <- chisq.test(table(final$subtype,final$sex))$p.value
tbl1[Characteristic=="Sex, n (%)", P:=sprintf("%.2g",psex)]
# 组织学
tbl1 <- add_line(tbl1, c("Histology, n (%)","","","","",""))
for (lv in c("LUAD","LUSC","Adenosquamous","Other")) {
  r <- c(paste0("  ",lv))
  for (s in c("IMS1","IMS2","IMS3","IMS4")) { sub<-final[subtype==s]; n<-sum(sub$histology_group==lv); r<-c(r,sprintf("%d (%.1f)",n,100*n/nrow(sub))) }
  n<-sum(final$histology_group==lv); r<-c(r,sprintf("%d (%.1f)",n,100*n/nrow(final)))
  tbl1 <- add_line(tbl1, r)
}
phist <- chisq.test(table(final$subtype,final$histology_group))$p.value
tbl1[Characteristic=="Histology, n (%)", P:=sprintf("%.1e",phist)]
# 治疗臂
tbl1 <- add_line(tbl1, c("Treatment arm, n (%)","","","","",""))
for (lv in c("EA5142","A081105","E4512")) {
  r <- c(paste0("  ",lv))
  for (s in c("IMS1","IMS2","IMS3","IMS4")) { sub<-final[subtype==s]; n<-sum(sub$treatment_arm==lv); r<-c(r,sprintf("%d (%.1f)",n,100*n/nrow(sub))) }
  n<-sum(final$treatment_arm==lv); r<-c(r,sprintf("%d (%.1f)",n,100*n/nrow(final)))
  tbl1 <- add_line(tbl1, r)
}
parm <- chisq.test(table(final$subtype,final$treatment_arm))$p.value
tbl1[Characteristic=="Treatment arm, n (%)", P:=sprintf("%.1e",parm)]
# TMB, FGA, Immune
for (vinfo in list(c("tmb_per_mb","TMB, median (IQR)"), c("FGA","FGA, median (IQR)"),
                   c("ImmuneScore","ImmuneScore, median (IQR)"))) {
  rr <- make_row_cont(final, vinfo[1], vinfo[2])
  tbl1 <- add_line(tbl1, rr$row, sprintf("%.1e", rr$p))
}

ft1 <- flextable(tbl1)
ft1 <- set_header_labels(ft1, Characteristic="Characteristic",
                         IMS1=sprintf("IMS1\n(n=%d)",n_hdr[1]), IMS2=sprintf("IMS2\n(n=%d)",n_hdr[2]),
                         IMS3=sprintf("IMS3\n(n=%d)",n_hdr[3]), IMS4=sprintf("IMS4\n(n=%d)",n_hdr[4]),
                         Overall=sprintf("Overall\n(n=%d)",nrow(final)), P="P value")
ft1 <- autofit(ft1); ft1 <- fontsize(ft1, size=8, part="all")
ft1 <- bold(ft1, part="header")

doc1 <- read_docx()
doc1 <- body_add_par(doc1, "Table 1. Clinical and molecular characteristics of the ALCHEMIST discovery cohort by integrated molecular subtype (n=935).", style="Normal")
doc1 <- body_add_flextable(doc1, ft1)
doc1 <- body_add_par(doc1, "IMS, integrated molecular subtype; TMB, tumor mutational burden (mutations/Mb); FGA, fraction of genome altered; IQR, interquartile range. P values by Kruskal-Wallis test (continuous) or chi-square test (categorical).", style="Normal")
print(doc1, target=file.path(TABDIR,"Table1_cohort_characteristics.docx"))
cat("Saved Table1_cohort_characteristics.docx\n")

# ============ TABLE 2: MSK多变量Cox ============
msk <- readRDS(file.path(RESULTS, "msk_survival_objects.rds"))
cox <- msk$cox
sc <- summary(cox)
cox_dt <- data.table(
  Variable=rownames(sc$coefficients),
  HR=sprintf("%.2f", sc$coefficients[,"exp(coef)"]),
  CI=sprintf("%.2f-%.2f", sc$conf.int[,"lower .95"], sc$conf.int[,"upper .95"]),
  P=sapply(sc$coefficients[,"Pr(>|z|)"], function(p) if(p<0.001) "<0.001" else sprintf("%.3f",p)))
# 美化变量名
cox_dt[, Variable := gsub("IMS","IMS ", Variable)]
cox_dt[, Variable := gsub("histology_group","Histology: ", Variable)]
cox_dt[, Variable := gsub("SAMPLE_TYPE","Sample: ", Variable)]
ft2 <- flextable(cox_dt)
ft2 <- set_header_labels(ft2, Variable="Variable", HR="Hazard ratio", CI="95% CI", P="P value")
ft2 <- autofit(ft2); ft2 <- fontsize(ft2, size=9, part="all"); ft2 <- bold(ft2, part="header")
doc2 <- read_docx()
doc2 <- body_add_par(doc2, "Table 2. Multivariable Cox proportional hazards analysis of overall survival by integrated molecular subtype in the MSK-IMPACT validation cohort (n=5,958; reference IMS1).", style="Normal")
doc2 <- body_add_flextable(doc2, ft2)
doc2 <- body_add_par(doc2, sprintf("Concordance index = %.3f. HR, hazard ratio; CI, confidence interval. Subtypes assigned by transfer of the genomic-axis signature.", sc$concordance[1]), style="Normal")
print(doc2, target=file.path(TABDIR,"Table2_msk_cox.docx"))
cat("Saved Table2_msk_cox.docx\n")

# ============ STable1: 驱动基因频率 ============
dbs <- fread(file.path(RESULTS, "driver_freq_by_subtype.csv"))
setnames(dbs, "p", "P_value")
dbs[, P_value := sapply(P_value, function(p) if(is.na(p)) "NA" else if(p<0.001) "<0.001" else sprintf("%.3f",p))]
fts1 <- flextable(dbs); fts1 <- autofit(fts1); fts1 <- fontsize(fts1, size=9, part="all"); fts1 <- bold(fts1, part="header")
docs1 <- read_docx()
docs1 <- body_add_par(docs1, "Supplementary Table S1. Driver gene mutation frequency (%) by integrated molecular subtype in ALCHEMIST (n=935). P values by Fisher exact test with Monte Carlo simulation.", style="Normal")
docs1 <- body_add_flextable(docs1, fts1)
print(docs1, target=file.path(TABDIR,"TableS1_driver_frequency.docx"))
cat("Saved TableS1_driver_frequency.docx\n")

cat("\n=== Tables complete ===\n")


# ---- (merged from 12_supp_tables.R) ----
# =============================================================
# 12_supp_tables.R  —  STable2 (Hallmark pathways by subtype)
# =============================================================
suppressMessages({ library(data.table); library(officer); library(flextable) })

hm <- fread(file.path(RESULTS, "hallmark_by_subtype.csv"))
setorder(hm, fdr)
hm_top <- hm[1:25]
hm_top[, `:=`(IMS1=round(IMS1,3), IMS2=round(IMS2,3), IMS3=round(IMS3,3), IMS4=round(IMS4,3))]
hm_top[, FDR := sapply(fdr, function(p) if(p<0.001) format(p,scientific=TRUE,digits=2) else sprintf("%.3f",p))]
hm_top[, p := NULL][, fdr := NULL]
setnames(hm_top, "pathway", "Hallmark_pathway")

ft <- flextable(hm_top)
ft <- autofit(ft); ft <- fontsize(ft, size=8, part="all"); ft <- bold(ft, part="header")
doc <- read_docx()
doc <- body_add_par(doc, "Supplementary Table S2. Top 25 differentially active Hallmark pathways across integrated molecular subtypes (ALCHEMIST, ssGSEA median scores; Kruskal-Wallis with Benjamini-Hochberg FDR).", style="Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target=file.path(TABDIR,"TableS2_pathways.docx"))
cat("Saved TableS2_pathways.docx\n")
