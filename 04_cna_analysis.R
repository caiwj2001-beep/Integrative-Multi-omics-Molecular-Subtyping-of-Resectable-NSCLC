# =============================================================
# 03_cna_analysis.R
# 拷贝数分析: 基因级增益/缺失, 基因组不稳定性(FGA), 臂级事件
# =============================================================
source("config.R")
suppressMessages({ library(data.table); library(matrixStats) })

# ---- 读取CNA基因级矩阵 (绝对拷贝数, ASCAT) ----
cat("Reading CNA gene matrix...\n")
cna <- fread(file.path(MATDIR, "cna_gene.tsv.gz"))
gene_ids <- cna$gene_id
gene_names <- cna$gene_name
mat <- as.matrix(cna[, -c("gene_id","gene_name")])
rownames(mat) <- gene_ids
cat("CNA matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")
cat("NA fraction:", round(mean(is.na(mat)),3), "\n")

# ---- 相对拷贝数 (相对倍性2) ----
# ASCAT copy_number已是整数绝对拷贝数; 阈值定义增益/缺失
# 缺失: CN<=1; 增益: CN>=3; 高扩增: CN>=5; 深缺失: CN=0
call_mat <- matrix(0L, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
call_mat[mat == 0] <- -2L           # deep deletion
call_mat[mat == 1] <- -1L           # loss
call_mat[mat >= 3 & mat < 5] <- 1L  # gain
call_mat[mat >= 5] <- 2L            # amplification
call_mat[is.na(mat)] <- NA

# ---- 基因组改变分数 (FGA proxy): 有CN变化基因比例 ----
fga <- colMeans(abs(call_mat) >= 1, na.rm = TRUE)
gain_frac <- colMeans(call_mat >= 1, na.rm = TRUE)
loss_frac <- colMeans(call_mat <= -1, na.rm = TRUE)

fga_dt <- data.table(patient_id = colnames(mat),
                     FGA = round(fga,4),
                     gain_fraction = round(gain_frac,4),
                     loss_fraction = round(loss_frac,4))
fwrite(fga_dt, file.path(RESULTS, "cna_genomic_instability.csv"))
cat("\n=== Genomic instability (FGA) ===\n")
cat(sprintf("FGA: median=%.3f, mean=%.3f, range=%.3f-%.3f\n",
            median(fga), mean(fga), min(fga), max(fga)))

# ---- 已知NSCLC焦点CNA靶基因 ----
focal_targets <- c("EGFR","MYC","MDM2","CDK4","MET","TERT","CCND1","CCNE1","SOX2",
                   "FGFR1","KRAS","NKX2-1","MECOM","CDKN2A","CDKN2B","PTEN","RB1",
                   "STK11","SMAD4","MTAP","ERBB2","MDM4","YEATS4","PIK3CA")
# gene_name索引
name2id <- setNames(gene_ids, gene_names)
focal_freq <- rbindlist(lapply(focal_targets, function(g){
  gid <- name2id[g]
  if (is.na(gid) || !(gid %in% rownames(call_mat))) return(NULL)
  row <- call_mat[gid, ]
  data.table(gene=g,
             amp_pct = round(100*mean(row>=2, na.rm=TRUE),1),
             gain_pct = round(100*mean(row>=1, na.rm=TRUE),1),
             loss_pct = round(100*mean(row<=-1, na.rm=TRUE),1),
             deepdel_pct = round(100*mean(row<=-2, na.rm=TRUE),1))
}))
setorder(focal_freq, -amp_pct)
fwrite(focal_freq, file.path(RESULTS, "cna_focal_frequency.csv"))
cat("\n=== Focal CNA target frequencies ===\n"); print(focal_freq)

# ---- 保存calls供整合分型 ----
saveRDS(list(call_mat=call_mat, cn_mat=mat, gene_id2name=setNames(gene_names, gene_ids)),
        file.path(RESULTS, "cna_processed.rds"))
cat("\n=== CNA analysis complete ===\n")
