suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))

# 使用annotated矩阵（含表头）
csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE57338_expression_matrix_annotated.csv"
if(!file.exists(csv)) {
  csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE57338_expression_matrix_raw.csv"
}
message("▶ 使用: ", csv)
df <- read.csv(csv, stringsAsFactors=FALSE, check.names=FALSE)
message("▶ 维度: ", nrow(df), " x ", ncol(df))
message("▶ 列名前5: ", paste(head(colnames(df), 5), collapse=", "))
message("▶ 行名样例: ", paste(head(df[,1], 3), collapse=", "))

# 假设第一列是基因名，其余是样本
genes <- df[,1]
expr <- as.matrix(df[, -1])
rownames(expr) <- genes
message("▶ 表达矩阵: ", nrow(expr), " 基因 x ", ncol(expr), " 样本")

# 转置为基因×样本（WGCNA标准格式）
datExpr <- t(expr)
message("▶ WGCNA输入: ", ncol(datExpr), " 样本 x ", nrow(datExpr), " 基因")

# 过滤低表达/低变异
gsg <- goodSamplesGenes(datExpr, verbose=3)
datExpr <- datExpr[, gsg$goodGenes]
message("▶ 过滤后: ", ncol(datExpr), " 基因")

# 软阈值
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5)
message("▶ 建议power: ", sft$powerEstimate)

# 构建网络
net <- blockwiseModules(datExpr, power=ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate),
  TOMType="unsigned", minModuleSize=30, reassignThreshold=0, mergeCutHeight=0.25,
  numericLabels=TRUE, pamRespectsDendro=FALSE, verbose=3)

message("▶ 模块数: ", length(unique(net$colors)))
message("▶ 模块大小分布:")
print(table(net$colors))

# NDUFB7模块归属
if("NDUFB7" %in% colnames(datExpr)) {
  ndufb7_mod <- net$colors["NDUFB7"]
  message("▶ NDUFB7所属模块: ", ndufb7_mod)
  mod_genes <- names(net$colors)[net$colors==ndufb7_mod]
  message("▶ 同模块基因数: ", length(mod_genes))
  message("▶ 同模块基因样例: ", paste(head(mod_genes, 20), collapse=", "))
  
  # 检查同模块是否有铁死亡基因
  ferro <- c("FTL","FTH1","SLC7A11","GPX4","SAT1","ACSL4","NFE2L2","KEAP1","LPCAT3","PTGS2","ALOX15","GLS2")
  ferro_in_mod <- ferro[ferro %in% mod_genes]
  message("▶ 同模块铁死亡基因: ", paste(ferro_in_mod, collapse=", "))
  
  # 模块特征基因（ME）与NDUFB7表达的相关
  MEs <- net$MEs
  ndufb7_expr <- datExpr[, "NDUFB7"]
  me_cor <- apply(MEs, 2, function(me) cor(me, ndufb7_expr, use="complete.obs"))
  top_me <- names(me_cor)[order(-abs(me_cor))][1:3]
  message("▶ 与NDUFB7最相关的模块特征基因: ", paste(top_me, "r=", round(me_cor[top_me], 3), collapse="; "))
}

# 保存
out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_WGCNA_Modules_Fix.RDS"
saveRDS(net, out)
message("✅ 保存: ", out)

# 模块基因列表
mod_df <- data.frame(Gene=names(net$colors), Module=net$colors, stringsAsFactors=FALSE)
out2 <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_WGCNA_Module_Genes_Fix.csv"
write.csv(mod_df, out2, row.names=FALSE)
message("✅ 保存模块基因: ", out2)
