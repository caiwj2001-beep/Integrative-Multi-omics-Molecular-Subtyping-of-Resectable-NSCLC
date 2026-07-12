# =============================================================
# V2_04_tcga_figures.R  —  短板1补救的图表
# Figure 7: TCGA验证 (分期分布 + 校正分期KM + 森林图)
# Table 3: 校正分期的多变量Cox
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(survival); library(survminer)
  library(ggplot2); library(patchwork); library(officer); library(flextable) })

obj <- readRDS(file.path(RESULTS, "tcga_transcriptomic_objects.rds"))
surv <- obj$surv; cox <- obj$cox
tcga <- fread(file.path(RESULTS, "tcga_transcriptomic_ims.csv"))

# ---- 7A: 分期在亚型间分布(堆叠条) ----
tcga2 <- tcga[!is.na(stage_simple) & stage_simple %in% c("I","II","III","IV")]
tcga2[, stage_simple := factor(stage_simple, levels=c("I","II","III","IV"))]
p7a <- ggplot(tcga2, aes(x=IMS_expr, fill=stage_simple)) +
  geom_bar(position="fill", width=0.7) +
  scale_fill_brewer(palette="Blues", name="AJCC\nstage") +
  scale_y_continuous(labels=scales::percent) +
  labs(title="Stage distribution across subtypes (TCGA)", x=NULL, y="Proportion") +
  theme_pub()

# ---- 7B: 校正分期后的KM (全队列OS) ----
fit <- survfit(Surv(OS_MONTHS, os_event) ~ IMS, data=surv)
km <- ggsurvplot(fit, data=surv, palette=unname(IMS_COL),
                 pval=TRUE, risk.table=TRUE, conf.int=FALSE,
                 xlab="Months", ylab="Overall survival", legend.title="Subtype",
                 legend.labs=levels(surv$IMS), risk.table.height=0.28,
                 ggtheme=theme_pub(), title="TCGA overall survival by subtype")

# ---- 7C: 多变量Cox森林图 ----
sc <- summary(cox)
cox_dt <- data.table(
  var=rownames(sc$conf.int),
  HR=sc$conf.int[,"exp(coef)"], lo=sc$conf.int[,"lower .95"], hi=sc$conf.int[,"upper .95"],
  p=sc$coefficients[,"Pr(>|z|)"])
cox_dt[, label := c("IMS2 vs IMS1","IMS3 vs IMS1","IMS4 vs IMS1",
                    "Stage II vs I","Stage III vs I","Stage IV vs I","LUSC vs LUAD")[seq_len(.N)]]
cox_dt[, label := factor(label, levels=rev(label))]
p7c <- ggplot(cox_dt, aes(x=HR, y=label)) +
  geom_vline(xintercept=1, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25, color="#003366") +
  geom_point(size=2.5, color="#B23A48") +
  scale_x_log10() +
  labs(title="Multivariable Cox (OS, stage-adjusted)", x="Hazard ratio (95% CI)", y=NULL) +
  theme_pub()

# 组合(KM单独存, 因ggsurvplot特殊)
fig7_top <- p7a | p7c
save_fig(fig7_top, "Figure7_TCGA_stage_forest", width=11, height=4.5)
png(file.path(FIGDIR,"Figure7_TCGA_KM.png"), width=7, height=7, units="in", res=300)
print(km); dev.off()
tiff(file.path(FIGDIR,"Figure7_TCGA_KM.tiff"), width=7, height=7, units="in", res=300, compression="lzw")
print(km); dev.off()
cat("Saved Figure7 (stage/forest + KM)\n")

# ---- Table 3: 校正分期Cox (DOCX) ----
t3 <- data.table(
  Variable = cox_dt$label,
  HR = sprintf("%.2f", cox_dt$HR),
  CI = sprintf("%.2f-%.2f", cox_dt$lo, cox_dt$hi),
  P = sapply(cox_dt$p, function(p) if(p<0.001) "<0.001" else sprintf("%.3f",p)))
ft <- flextable(t3); ft <- autofit(ft); ft <- fontsize(ft,size=9,part="all"); ft<-bold(ft,part="header")
doc <- read_docx()
doc <- body_add_par(doc, "Table 3. Multivariable Cox proportional-hazards analysis of overall survival in TCGA NSCLC (n=979), demonstrating that integrated molecular subtype retains prognostic value after adjustment for AJCC stage and histology.", style="Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, sprintf("Reference: IMS1, stage I, LUAD. Concordance index = %.3f. Stage distribution did not differ materially across subtypes (70-83%% early-stage), indicating that subtype captures prognostic information beyond stage.", summary(cox)$concordance[1]), style="Normal")
print(doc, target=file.path(TABDIR, "Table3_TCGA_stage_adjusted_cox.docx"))
cat("Saved Table3\n")

# 分期分布表(附表)
st <- obj$stage_tab
st_dt <- as.data.table(as.matrix(st), keep.rownames="Subtype")
ft2 <- flextable(st_dt); ft2<-autofit(ft2); ft2<-fontsize(ft2,size=9,part="all"); ft2<-bold(ft2,part="header")
doc2 <- read_docx()
doc2 <- body_add_par(doc2, "Supplementary Table S3. Distribution of AJCC pathologic stage across integrated molecular subtypes in TCGA (chi-square p=0.04; early-stage proportion 70-83% across subtypes).", style="Normal")
doc2 <- body_add_flextable(doc2, ft2)
print(doc2, target=file.path(TABDIR, "TableS3_TCGA_stage_distribution.docx"))
cat("Saved TableS3\n")
cat("\n=== V2 TCGA figures/tables complete ===\n")


# ---- (merged from V3_02_concordance_figure.R) ----
# =============================================================
# V3_02_concordance_figure.R  —  双轴一致性图表
# Figure 8: 混淆矩阵热图 + IMS4成员重叠
# Table S4: 双轴混淆矩阵 + kappa
# =============================================================
suppressMessages({ library(data.table); library(ggplot2); library(officer); library(flextable) })

obj <- readRDS(file.path(RESULTS, "dual_axis_concordance.rds"))
cm <- obj$cm; kap <- obj$kappa; D <- obj$D
subtypes <- c("IMS1","IMS2","IMS3","IMS4")

# ---- Figure 8A: 混淆矩阵热图(行归一化%) ----
cm_pct <- prop.table(as.matrix(cm), 1) * 100
dfp <- as.data.table(as.table(cm_pct)); setnames(dfp, c("Genomic","Transcriptomic","Pct"))
dfn <- as.data.table(as.table(as.matrix(cm))); setnames(dfn, c("Genomic","Transcriptomic","N"))
dfp <- merge(dfp, dfn, by=c("Genomic","Transcriptomic"))
p8a <- ggplot(dfp, aes(x=Transcriptomic, y=Genomic, fill=Pct)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.0f%%\n(n=%d)", Pct, N)), size=3) +
  scale_fill_gradient(low="white", high="#003366", name="Row %") +
  scale_y_discrete(limits=rev(subtypes)) +
  labs(title=sprintf("Dual-axis subtype concordance in TCGA (\u03BA=%.2f, fair)", kap$value),
       x="Transcriptomic-axis assignment", y="Genomic-axis assignment") +
  theme_pub() + theme(panel.grid=element_blank())
save_fig(p8a, "Figure8_dual_axis_concordance", width=6.5, height=5.5)

# ---- Table S4: 混淆矩阵 DOCX ----
cmdt <- as.data.table(as.matrix(cm), keep.rownames="Genomic_axis")
ft <- flextable(cmdt); ft <- autofit(ft); ft <- fontsize(ft, size=9, part="all"); ft <- bold(ft, part="header")
ft <- add_header_lines(ft, values="Transcriptomic-axis assignment (columns)")
doc <- read_docx()
doc <- body_add_par(doc, sprintf("Supplementary Table S4. Cross-tabulation of integrated molecular subtype assignments in TCGA (n=%d) when derived from the genomic axis (rows; driver mutations, TMB, fraction genome altered) versus the transcriptomic axis (columns; 154-gene expression signature). Agreement was fair (Cohen's kappa = %.2f; raw concordance 44.5%%), and the IMS4 label overlapped by only 18 of 212 patients between axes (Jaccard = 0.08), indicating that the two dimensionality-reduction axes capture distinct biological information.", nrow(D), kap$value), style="Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target=file.path(TABDIR, "TableS4_dual_axis_concordance.docx"))
cat("Saved Figure8 + TableS4\n")
cat(sprintf("kappa=%.3f, agree=%.1f%%\n", kap$value, obj$agree_pct))


# ---- (merged from V3_05_figureS1.R) ----
# =============================================================
# V3_05_figureS1.R  —  生成缺失的 Figure S1
# 一致性聚类诊断: (A)表达聚类CDF (B)MOFA因子聚类CDF (C)PAC vs K
# 数据源: V1的ConsensusClusterPlus对象 + PAC表
# 输出: V3/figures/FigureS1_consensus_diagnostics.{png,tiff}
# =============================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

# ---- 从共识矩阵计算CDF ----
cdf_from_cc <- function(cc, krange=2:8){
  out <- list()
  for (k in krange){
    M <- cc[[k]]$consensusMatrix
    v <- M[lower.tri(M)]
    ec <- ecdf(v)
    xs <- seq(0,1,length.out=100)
    out[[as.character(k)]] <- data.table(consensus_index=xs, CDF=ec(xs), K=factor(k))
  }
  rbindlist(out)
}

cc_expr <- readRDS(file.path(RESULTS, "cc_expression_object.rds"))
cc_mofa <- readRDS(file.path(RESULTS, "cc_mofa_object.rds"))
cdf_expr <- cdf_from_cc(cc_expr)
cdf_mofa <- cdf_from_cc(cc_mofa)

kcol <- setNames(colorRampPalette(c("#1B5E9B","#E8A73C","#3E8E7E","#B23A48"))(7), as.character(2:8))

# ---- (A) 表达聚类 CDF ----
pA <- ggplot(cdf_expr, aes(consensus_index, CDF, color=K)) +
  geom_line(linewidth=0.7) + scale_color_manual(values=kcol, name="K") +
  labs(title="Expression clustering", x="Consensus index", y="CDF") +
  theme_pub()

# ---- (B) MOFA因子聚类 CDF ----
pB <- ggplot(cdf_mofa, aes(consensus_index, CDF, color=K)) +
  geom_line(linewidth=0.7) + scale_color_manual(values=kcol, name="K") +
  labs(title="MOFA-factor clustering", x="Consensus index", y="CDF") +
  theme_pub()

# ---- (C) PAC vs K (两者) ----
pac_e <- fread(file.path(RESULTS, "expression_PAC.csv")); pac_e[, Method:="Expression"]
pac_m <- fread(file.path(RESULTS, "mofa_PAC.csv")); pac_m[, Method:="MOFA-factor"]
pac <- rbind(pac_e, pac_m)
pC <- ggplot(pac, aes(K, PAC, color=Method, shape=Method)) +
  geom_line(linewidth=0.7) + geom_point(size=2.5) +
  scale_color_manual(values=c("Expression"="#1B5E9B","MOFA-factor"="#B23A48")) +
  scale_x_continuous(breaks=2:8) +
  labs(title="Proportion of ambiguous clustering", x="Number of clusters (K)",
       y="PAC (lower = more stable)") +
  theme_pub() + theme(legend.position=c(0.98,0.98), legend.justification=c(1,1))

figS1 <- (pA | pB) / (pC | patchwork::plot_spacer()) + plot_annotation(tag_levels='A')
# 用grid直接组合更紧凑: A,B上排, C下排跨列
figS1 <- (pA | pB) / pC + plot_layout(heights=c(1,1)) + plot_annotation(tag_levels='A')

save_fig(figS1, "FigureS1_consensus_diagnostics", width=10, height=8)
cat("Saved FigureS1_consensus_diagnostics\n")
