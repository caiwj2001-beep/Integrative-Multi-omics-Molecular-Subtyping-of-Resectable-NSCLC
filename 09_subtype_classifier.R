# =============================================================
# 07_subtype_classifier.R (v2)
# 多组学特征分类器: 表达 + FGA + TMB + 焦点CNA + 免疫评分
# IMS3(高FGA)/IMS4(免疫热)由CNA/免疫定义, 必须纳入这些特征
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(matrixStats); library(randomForest) })
set.seed(20260707)

clin <- fread(file.path(DATA, "alchemist_master_clinical.csv"))
master <- fread(file.path(RESULTS, "master_with_subtypes.csv"))
core3 <- clin[core3 == TRUE, patient_id]
labeled <- master[!is.na(subtype), .(patient_id, subtype)]
unlabeled <- setdiff(core3, labeled$patient_id)
cat("core3:", length(core3), " labeled:", nrow(labeled), " to classify:", length(unlabeled), "\n")

# ---- 特征1: 表达 top1500高变 ----
ep <- readRDS(file.path(RESULTS, "expression_processed.rds"))
logmat <- ep$logmat; id2name <- ep$gene_id2name
sym <- id2name[rownames(logmat)]
dt <- data.table(sym=sym, idx=seq_len(nrow(logmat)), v=rowVars(logmat))
dt <- dt[!is.na(sym) & sym!=""]; setorder(dt, sym, -v)
keep <- dt[, .SD[1], by=sym]; symmat <- logmat[keep$idx, ]; rownames(symmat) <- keep$sym
vg <- names(sort(rowVars(symmat), decreasing=TRUE))[1:1500]
Xexpr <- t(symmat[vg, ])

# ---- 特征2: 基因组不稳定性 + 免疫 ----
tmb <- fread(file.path(MATDIR, "tmb.csv"))
fga <- fread(file.path(RESULTS, "cna_genomic_instability.csv"))
est <- fread(file.path(RESULTS, "estimate_scores.csv"))
genomic <- Reduce(function(a,b) merge(a,b,by="patient_id",all=TRUE),
                  list(tmb[,.(patient_id, tmb_per_mb, nonsyn_count)],
                       fga[,.(patient_id, FGA, gain_fraction, loss_fraction)],
                       est[,.(patient_id, ImmuneScore, StromalScore, ESTIMATEScore)]))

# ---- 特征3: 焦点CNA (定义IMS3的关键) ----
cp <- readRDS(file.path(RESULTS, "cna_processed.rds"))
call_mat <- cp$call_mat; cna_n2i <- setNames(names(cp$gene_id2name), cp$gene_id2name)
focal <- c("EGFR","MYC","MDM2","CDK4","MET","TERT","CCND1","CCNE1","SOX2","FGFR1",
           "KRAS","NKX2-1","CDKN2A","PTEN","RB1","SMAD4","MTAP","ERBB2")
focal_ids <- cna_n2i[focal]; focal_ids <- focal_ids[!is.na(focal_ids)]
Xcna <- t(call_mat[focal_ids, , drop=FALSE])
colnames(Xcna) <- paste0("cna_", names(focal_ids))

# ---- 组装特征矩阵 ----
common_all <- Reduce(intersect, list(rownames(Xexpr), genomic$patient_id, rownames(Xcna)))
common_all <- intersect(common_all, core3)
gm <- as.data.frame(genomic[match(common_all, patient_id)]); rownames(gm) <- common_all
gm$patient_id <- NULL
X <- cbind(scale(Xexpr[common_all,]), scale(as.matrix(gm)), Xcna[common_all,])
X[is.na(X)] <- 0  # CNA NA填0(中性)
cat("Feature matrix:", nrow(X), "x", ncol(X), "\n")

train_ids <- intersect(labeled$patient_id, rownames(X))
pred_ids <- intersect(unlabeled, rownames(X))
ytr <- labeled[match(train_ids, patient_id), subtype]

# ---- RF (平衡类权重, 因IMS4小) ----
cw <- 1/table(ytr); cw <- cw/sum(cw)
rf <- randomForest(x=X[train_ids,], y=factor(ytr), ntree=2000,
                   classwt=as.numeric(cw), importance=TRUE)
cat("\n=== RF OOB confusion ===\n"); print(rf$confusion)
cat("OOB error:", round(rf$err.rate[nrow(rf$err.rate),"OOB"],3), "\n")

# ---- 预测 ----
if (length(pred_ids)>0) {
  pred <- predict(rf, X[pred_ids,,drop=FALSE])
  pc <- predict(rf, X[pred_ids,,drop=FALSE], type="prob")
  new_labels <- data.table(patient_id=pred_ids, subtype=as.character(pred),
                           confidence=round(apply(pc,1,max),3), source="RF_predicted")
} else new_labels <- data.table()
labeled[, `:=`(confidence=1.0, source="MOFA")]
all_subtypes <- rbind(labeled, new_labels)
cat("\n=== Final core3 subtype distribution ===\n"); print(table(all_subtypes$subtype))
cat("By source:\n"); print(table(all_subtypes$subtype, all_subtypes$source))

# ---- 合并 & 保存 ----
final <- merge(clin, all_subtypes, by="patient_id")
final <- merge(final, tmb[,.(patient_id,tmb_per_mb,nonsyn_count)], by="patient_id", all.x=TRUE)
final <- merge(final, fga, by="patient_id", all.x=TRUE)
final <- merge(final, est, by="patient_id", all.x=TRUE)
fwrite(final, file.path(RESULTS, "core3_subtypes_final.csv"))
saveRDS(rf, file.path(RESULTS, "rf_classifier.rds"))

cat("\n=== 亚型×治疗臂 (core3) ===\n"); tt<-table(final$subtype,final$treatment_arm); print(tt)
print(suppressWarnings(chisq.test(tt)))
cat("\n=== 亚型×组织学 ===\n"); print(table(final$subtype, final$histology_group))

# 变量重要性top20
imp <- importance(rf, type=1)
imp_dt <- data.table(feature=rownames(imp), MDA=imp[,1])[order(-MDA)][1:20]
fwrite(imp_dt, file.path(RESULTS, "classifier_importance.csv"))
cat("\n=== Top classifier features ===\n"); print(imp_dt)
cat("\n=== Classifier v2 complete ===\n")
