suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))

# 加载GSE57338
csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE57338_NDUFB7_expression.csv"
if(!file.exists(csv)) {
  message("❌ GSE57338不存在"); quit(status=1)
}
df <- read.csv(csv, stringsAsFactors=FALSE)
message("▶ GSE57338: ", nrow(df), " 样本 x ", ncol(df), " 基因")

# 提取表达矩阵（除样本名列外的所有数值列）
expr <- as.matrix(df[, -1])
rownames(expr) <- df[,1]
message("▶ 表达矩阵: ", nrow(expr), " x ", ncol(expr))

# 转置为基因×样本
datExpr <- t(expr)
message("▶ WGCNA输入: ", nrow(datExpr), " 基因 x ", ncol(datExpr), " 样本")

# 过滤低变异基因
var_genes <- apply(datExpr, 1, var, na.rm=TRUE)
datExpr <- datExpr[var_genes > quantile(var_genes, 0.25, na.rm=TRUE), ]
message("▶ 过滤后: ", nrow(datExpr), " 基因")

# 选择软阈值
powers <- c(1:20)
sft <- pickSoftThreshold(t(datExpr), powerVector=powers, verbose=5)
message("▶ 建议软阈值: ", sft$powerEstimate)

# 构建网络
net <- blockwiseModules(t(datExpr), power=ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate),
  TOMType="unsigned", minModuleSize=30, reassignThreshold=0, mergeCutHeight=0.25,
  numericLabels=TRUE, pamRespectsDendro=FALSE, verbose=3)
message("▶ 模块数: ", length(unique(net$colors)))
message("▶ 模块大小:")
print(table(net$colors))

# NDUFB7模块归属
ndufb7_mod <- net$colors["NDUFB7"]
message("▶ NDUFB7所属模块: ", ndufb7_mod)
mod_genes <- names(net$colors)[net$colors==ndufb7_mod]
message("▶ 同模块基因数: ", length(mod_genes))
message("▶ 同模块Top基因 (按模块内连接度):")
kme <- net$MEs
if("NDUFB7" %in% rownames(kme)) {
  message("  NDUFB7模块特征基因表达: ", paste(head(kme["NDUFB7", ], 5), collapse=", "))
}

# 保存
out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_WGCNA_Modules.RDS"
saveRDS(net, out)
message("✅ 保存: ", out)

# 模块-性状关联（如果临床数据可用）
# 简化：输出模块基因列表
mod_df <- data.frame(Gene=names(net$colors), Module=net$colors, stringsAsFactors=FALSE)
out2 <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_WGCNA_Module_Genes.csv"
write.csv(mod_df, out2, row.names=FALSE)
message("✅ 保存模块基因: ", out2)
