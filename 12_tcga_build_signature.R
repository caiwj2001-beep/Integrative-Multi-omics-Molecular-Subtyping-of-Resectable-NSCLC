# =============================================================
# V2_02_build_expression_signature.R
# 从ALCHEMIST表达分型提取各亚型marker基因(用于TCGA转录组迁移)
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(matrixStats); library(limma) })

# ALCHEMIST表达(symbol级) + 分型
ep <- readRDS(file.path(RESULTS, "expression_processed.rds"))
logmat <- ep$logmat; id2name <- ep$gene_id2name
sym <- id2name[rownames(logmat)]
dt <- data.table(sym=sym, idx=seq_len(nrow(logmat)), v=rowVars(logmat))
dt <- dt[!is.na(sym)&sym!=""]; setorder(dt,sym,-v); keep<-dt[,.SD[1],by=sym]
symmat <- logmat[keep$idx,]; rownames(symmat)<-keep$sym

master <- fread(file.path(RESULTS, "master_with_subtypes.csv"))[!is.na(subtype)]
ids <- intersect(master$patient_id, colnames(symmat))
master <- master[match(ids, patient_id)]
X <- symmat[, ids]
grp <- factor(master$subtype)
cat("Signature discovery cohort:", length(ids), "samples,", nlevels(grp), "subtypes\n")
print(table(grp))

# limma: 每个亚型 vs 其余, 提取上调marker
design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
fit <- lmFit(X, design)
markers <- list()
for (s in levels(grp)) {
  others <- setdiff(levels(grp), s)
  contr <- makeContrasts(contrasts=paste0(s, " - (", paste(others, collapse="+"), ")/", length(others)),
                         levels=design)
  f2 <- contrasts.fit(fit, contr); f2 <- eBayes(f2)
  tt <- topTable(f2, number=Inf, sort.by="t")
  up <- rownames(tt)[tt$logFC > 0.5 & tt$adj.P.Val < 0.01][1:50]  # top50上调
  markers[[s]] <- up[!is.na(up)]
  cat(sprintf("%s: %d marker genes\n", s, length(markers[[s]])))
}

# 保存marker基因 + 各亚型质心(表达谱)
all_markers <- unique(unlist(markers))
cat("Total unique markers:", length(all_markers), "\n")
# 各亚型在marker基因上的质心(均值表达)
centroid_expr <- sapply(levels(grp), function(s){
  rowMeans(X[all_markers, master$subtype==s])
})
saveRDS(list(markers=markers, all_markers=all_markers, centroid_expr=centroid_expr),
        file.path(RESULTS, "expression_signature.rds"))
fwrite(data.table(gene=all_markers), file.path(RESULTS, "signature_genes.csv"))
cat("Saved expression signature (", length(all_markers), "genes)\n")
