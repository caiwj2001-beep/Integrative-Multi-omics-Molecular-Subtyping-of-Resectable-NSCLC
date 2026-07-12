# =============================================================
# 09_figures.R
# 投稿主图生成 (PNG + TIFF 300dpi)
# Fig1: 队列与多组学概览 + MOFA方差分解
# Fig2: 突变景观(oncoplot)
# Fig3: 亚型定义热图(多组学)
# Fig4: 亚型生物学特征(TMB/FGA/免疫/通路)
# Fig5: 驱动基因与CNA的亚型特异性
# Fig6: MSK外部验证 + 生存
# =============================================================
source("config.R")
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(survival); library(survminer)
  library(maftools); library(RColorBrewer); library(grid)
})

IMS_COL <- c(IMS1="#1B5E9B", IMS2="#E8A73C", IMS3="#B23A48", IMS4="#3E8E7E")
final <- fread(file.path(RESULTS, "core3_subtypes_final.csv"))
final <- final[!is.na(subtype)]
setkey(final, patient_id)

# ============ FIG 1: 队列概览 + MOFA方差 ============
# 1a: 组学可用性堆叠
clin <- fread(file.path(DATA, "alchemist_master_clinical.csv"))
omics_counts <- data.table(
  Omics=c("WXS\nMutation","RNA-Seq\nExpression","WGS\nCNA","miRNA","isoform"),
  N=c(sum(clin$has_maf), sum(clin$has_expr), sum(clin$has_cna_gene),
      sum(clin$has_mirna), sum(clin$has_isoform)))
omics_counts[, Omics:=factor(Omics, levels=Omics)]
p1a <- ggplot(omics_counts, aes(x=Omics, y=N, fill=Omics)) +
  geom_col(width=0.7, show.legend=FALSE) +
  geom_text(aes(label=N), vjust=-0.4, size=3.2) +
  scale_fill_manual(values=colorRampPalette(c("#003366","#6BAED6"))(5)) +
  labs(title="Multi-omics data availability", x=NULL, y="Patients") +
  ylim(0,1250) + theme_pub()

# 1b: MOFA方差分解
vexp <- readRDS(file.path(RESULTS, "mofa_variance.rds"))
r2 <- vexp$r2_per_factor[[1]]
r2dt <- as.data.table(reshape2::melt(r2))
setnames(r2dt, c("Factor","View","R2"))
r2dt <- r2dt[Factor %in% paste0("Factor",1:8)]
p1b <- ggplot(r2dt, aes(x=View, y=Factor, fill=R2)) +
  geom_tile(color="white") +
  scale_fill_gradient(low="white", high="#003366", name="Var%") +
  labs(title="MOFA variance decomposition", x=NULL, y=NULL) +
  theme_pub() + theme(axis.text.x=element_text(angle=45,hjust=1))

# 1c: 亚型分布饼/柱
sub_n <- final[, .N, by=subtype][order(subtype)]
p1c <- ggplot(sub_n, aes(x=subtype, y=N, fill=subtype)) +
  geom_col(width=0.7, show.legend=FALSE) +
  geom_text(aes(label=N), vjust=-0.4, size=3.5) +
  scale_fill_manual(values=IMS_COL) +
  labs(title="Integrated molecular subtypes", x=NULL, y="Patients") +
  ylim(0,720) + theme_pub()

# 1d: 亚型×组织学
ht <- final[, .N, by=.(subtype, histology_group)]
p1d <- ggplot(ht, aes(x=subtype, y=N, fill=histology_group)) +
  geom_col(position="fill", width=0.7) +
  scale_fill_brewer(palette="Set2", name="Histology") +
  labs(title="Histology composition", x=NULL, y="Proportion") +
  theme_pub()

fig1 <- (p1a | p1b) / (p1c | p1d) + plot_annotation(tag_levels='A')
save_fig(fig1, "Figure1_cohort_overview", width=11, height=9)

# ============ FIG 4: 亚型生物学特征 ============
plot_box <- function(var, ylab, title, log=FALSE) {
  d <- final[!is.na(get(var))]
  p <- ggplot(d, aes(x=subtype, y=get(var), fill=subtype)) +
    geom_boxplot(outlier.size=0.4, width=0.6, show.legend=FALSE) +
    scale_fill_manual(values=IMS_COL) +
    labs(title=title, x=NULL, y=ylab) + theme_pub()
  if (log) p <- p + scale_y_log10()
  # 加Kruskal p
  kw <- kruskal.test(d[[var]] ~ d$subtype)
  p + annotate("text", x=1, y=Inf, hjust=0, vjust=1.3,
               label=sprintf("KW p=%.1e", kw$p.value), size=2.8)
}
p4a <- plot_box("tmb_per_mb","TMB (mut/Mb)","Tumor mutational burden", log=TRUE)
p4b <- plot_box("FGA","Fraction genome altered","Genomic instability")
p4c <- plot_box("ImmuneScore","ESTIMATE ImmuneScore","Immune infiltration")
p4d <- plot_box("StromalScore","ESTIMATE StromalScore","Stromal content")
fig4 <- (p4a|p4b)/(p4c|p4d) + plot_annotation(tag_levels='A')
save_fig(fig4, "Figure4_subtype_biology", width=10, height=8)

# ============ FIG 6: MSK验证 生存 ============
msk_surv <- readRDS(file.path(RESULTS, "msk_survival_objects.rds"))
surv_dt <- msk_surv$surv_dt
fit <- survfit(Surv(OS_months, OS_event) ~ IMS, data=surv_dt)
km <- ggsurvplot(fit, data=surv_dt, palette=unname(IMS_COL),
                 risk.table=TRUE, pval=TRUE, conf.int=FALSE,
                 xlab="Months", ylab="Overall survival",
                 legend.title="Subtype", legend.labs=names(IMS_COL),
                 risk.table.height=0.28, ggtheme=theme_pub(),
                 title="MSK-IMPACT external validation (n=5,958)")
png(file.path(FIGDIR,"Figure6_msk_survival.png"), width=7, height=7, units="in", res=300)
print(km); dev.off()
tiff(file.path(FIGDIR,"Figure6_msk_survival.tiff"), width=7, height=7, units="in", res=300, compression="lzw")
print(km); dev.off()
cat("Saved: Figure6_msk_survival (.png/.tiff)\n")

cat("\n=== Figures 1,4,6 complete ===\n")


# ---- (merged from 10_figures_part2.R) ----
# =============================================================
# 10_figures_part2.R
# Fig2: 突变景观oncoplot
# Fig3: 多组学整合热图 (亚型定义)
# Fig5: 驱动基因&CNA亚型特异性 + 通路活性
# =============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(maftools); library(grid)
})
IMS_COL <- c(IMS1="#1B5E9B", IMS2="#E8A73C", IMS3="#B23A48", IMS4="#3E8E7E")
final <- fread(file.path(RESULTS, "core3_subtypes_final.csv"))[!is.na(subtype)]

# ============ FIG 2: Oncoplot ============
laml <- readRDS(file.path(RESULTS, "maf_object.rds"))
# 注释: 亚型
anno <- final[, .(Tumor_Sample_Barcode=patient_id, IMS=subtype,
                  Histology=histology_group)]
laml2 <- tryCatch(subsetMaf(laml, tsb=anno$Tumor_Sample_Barcode, mafObj=TRUE),
                  error=function(e) laml)
# 更新临床注释
laml2@clinical.data <- merge(laml2@clinical.data, anno, by="Tumor_Sample_Barcode", all.x=TRUE)

top_genes <- c("TP53","EGFR","KRAS","KEAP1","STK11","NF1","RBM10","SMARCA4","BRAF",
               "PIK3CA","RB1","CDKN2A","ATM","ARID1A","SETD2","PTEN","NFE2L2","MET","MGA","MET")
top_genes <- unique(top_genes)
ann_cols <- list(IMS=IMS_COL,
                 Histology=c(LUAD="#66C2A5",LUSC="#FC8D62",Adenosquamous="#8DA0CB",Other="#E78AC3"))

png(file.path(FIGDIR,"Figure2_oncoplot.png"), width=10, height=8, units="in", res=300)
oncoplot(maf=laml2, genes=top_genes, clinicalFeatures=c("IMS","Histology"),
         sortByAnnotation=TRUE, annotationColor=ann_cols,
         gene_mar=7, fontSize=0.6, titleText="ALCHEMIST mutation landscape by subtype")
dev.off()
tiff(file.path(FIGDIR,"Figure2_oncoplot.tiff"), width=10, height=8, units="in", res=300, compression="lzw")
oncoplot(maf=laml2, genes=top_genes, clinicalFeatures=c("IMS","Histology"),
         sortByAnnotation=TRUE, annotationColor=ann_cols,
         gene_mar=7, fontSize=0.6, titleText="ALCHEMIST mutation landscape by subtype")
dev.off()
cat("Saved: Figure2_oncoplot\n")

# ============ FIG 3: 多组学整合热图 ============
# 用730全组学样本, 表达top差异基因 + 注释(亚型/组织学/TMB/FGA/免疫)
ep <- readRDS(file.path(RESULTS, "expression_processed.rds"))
logmat <- ep$logmat; id2name <- ep$gene_id2name
library(matrixStats)
sym <- id2name[rownames(logmat)]
dt <- data.table(sym=sym, idx=seq_len(nrow(logmat)), v=rowVars(logmat))
dt <- dt[!is.na(sym)&sym!=""]; setorder(dt,sym,-v); keep<-dt[,.SD[1],by=sym]
symmat <- logmat[keep$idx,]; rownames(symmat)<-keep$sym

master <- fread(file.path(RESULTS, "master_with_subtypes.csv"))[!is.na(subtype)]
ids <- intersect(master$patient_id, colnames(symmat))
master <- master[patient_id %in% ids]
setorder(master, subtype)
ids <- master$patient_id

# top50高变基因 for viz
vg <- names(sort(rowVars(symmat[,ids]),decreasing=TRUE))[1:50]
hm_mat <- t(scale(t(symmat[vg, ids])))
hm_mat[hm_mat > 3] <- 3; hm_mat[hm_mat < -3] <- -3

top_anno <- HeatmapAnnotation(
  Subtype=master$subtype,
  Histology=master$histology_group,
  TMB=master$tmb_per_mb,
  FGA=master$FGA,
  ImmuneScore=master$ImmuneScore,
  col=list(Subtype=IMS_COL,
           Histology=c(LUAD="#66C2A5",LUSC="#FC8D62",Adenosquamous="#8DA0CB",Other="#E78AC3"),
           TMB=colorRamp2(c(0,5,15),c("white","#FDBB84","#B30000")),
           FGA=colorRamp2(c(0,0.4,0.8),c("white","#9ECAE1","#08519C")),
           ImmuneScore=colorRamp2(c(0,2000,4000),c("#2166AC","white","#B2182B"))),
  annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(3,"mm"))

ht <- Heatmap(hm_mat, name="Expr\n(z)", top_annotation=top_anno,
              show_column_names=FALSE, cluster_columns=FALSE,
              column_split=master$subtype,
              col=colorRamp2(c(-3,0,3),c("#2166AC","white","#B2182B")),
              row_names_gp=gpar(fontsize=6), row_km=4,
              column_title_gp=gpar(fontsize=10, fontface="bold"))
png(file.path(FIGDIR,"Figure3_multiomics_heatmap.png"), width=11, height=8, units="in", res=300)
draw(ht, merge_legend=TRUE); dev.off()
tiff(file.path(FIGDIR,"Figure3_multiomics_heatmap.tiff"), width=11, height=8, units="in", res=300, compression="lzw")
draw(ht, merge_legend=TRUE); dev.off()
cat("Saved: Figure3_multiomics_heatmap\n")

# ============ FIG 5: 驱动基因/CNA/通路 亚型特异性 ============
# 5a: 驱动基因频率热图
dbs <- fread(file.path(RESULTS, "driver_freq_by_subtype.csv"))
dbs_m <- as.matrix(dbs[, .(IMS1,IMS2,IMS3,IMS4)]); rownames(dbs_m)<-dbs$gene
p5a_ht <- Heatmap(dbs_m, name="Mut %", col=colorRamp2(c(0,25,80),c("white","#FDBB84","#B30000")),
                  cluster_columns=FALSE, cluster_rows=TRUE,
                  cell_fun=function(j,i,x,y,w,h,f) grid.text(sprintf("%.0f",dbs_m[i,j]),x,y,gp=gpar(fontsize=7)),
                  column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=8),
                  column_title="Driver mutation frequency by subtype")
png(file.path(FIGDIR,"Figure5_driver_by_subtype.png"), width=6, height=8, units="in", res=300)
draw(p5a_ht); dev.off()
tiff(file.path(FIGDIR,"Figure5_driver_by_subtype.tiff"), width=6, height=8, units="in", res=300, compression="lzw")
draw(p5a_ht); dev.off()
cat("Saved: Figure5_driver_by_subtype\n")

cat("\n=== Figures 2,3,5 complete ===\n")


# ---- (merged from V3_04_redraw_figure6.R: redraw Fig6 with censor.size=1) ----
# =============================================================
# V3_04_redraw_figure6.R  —  重画Figure6 (方案B)
# 保留censoring记号但改小(censor.size=1) + 加粗生存曲线(size=1.4)
# 输出覆盖 V3/figures/Figure6_msk_survival.{png,tiff}
# =============================================================
suppressMessages({ library(data.table); library(survival); library(survminer) })

# 数据源: V1结果里的MSK生存对象(与原Figure6一致)
msk_surv <- readRDS(file.path(RESULTS, "msk_survival_objects.rds"))
surv_dt <- msk_surv$surv_dt
fit <- survfit(Surv(OS_months, OS_event) ~ IMS, data=surv_dt)

km <- ggsurvplot(
  fit, data = surv_dt,
  palette      = unname(IMS_COL),
  size         = 1.4,        # 加粗生存曲线 (默认1.0)
  censor.size  = 1,          # 缩小删失记号 (默认~4.5)
  censor.shape = "|",        # 标准竖线删失记号
  risk.table   = TRUE, pval = TRUE, conf.int = FALSE,
  xlab = "Months", ylab = "Overall survival",
  legend.title = "Subtype", legend.labs = names(IMS_COL),
  risk.table.height = 0.28, ggtheme = theme_pub(),
  title = "MSK-IMPACT external validation (n=5,958)")

png(file.path(FIGDIR, "Figure6_msk_survival.png"), width=7, height=7, units="in", res=300)
print(km); dev.off()
tiff(file.path(FIGDIR, "Figure6_msk_survival.tiff"), width=7, height=7, units="in", res=300, compression="lzw")
print(km); dev.off()
cat("Redrawn: Figure6_msk_survival (censor.size=1, line size=1.4)\n")
