suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

# 加载评分后的对象（如果B1已完成）或原始对象
rds_scored <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Scored.RDS"
rds_raw <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
obj <- NULL
if(file.exists(rds_scored)) {
  obj <- readRDS(rds_scored); message("▶ 从评分RDS加载")
} else if(file.exists(rds_raw)) {
  obj <- readRDS(rds_raw); message("▶ 从原始RDS加载")
} else {
  env <- new.env()
  load("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/GSE183852_DCM_Integrated.Robj", envir=env)
  for(v in ls(env)) if(inherits(env[[v]], "Seurat")) { obj <- env[[v]]; break }
  message("▶ 从.Robj加载")
}
if(is.null(obj)) { message("❌ 无对象"); quit(status=1) }

# 获取NDUFB7表达
ndufb7_expr <- tryCatch(LayerData(obj, assay="RNA", layer="data")["NDUFB7", ], error=function(e) {
  tryCatch(GetAssayData(obj, assay="RNA", slot="data")["NDUFB7", ], error=function(e2) NULL)
})
if(is.null(ndufb7_expr)) { message("❌ 无法提取NDUFB7"); quit(status=1) }

message("▶ 计算NDUFB7共表达 (全基因组，可能需要5-10分钟)...")
# 采样加速（如果细胞数>50000）
set.seed(42)
cells <- colnames(obj)
if(length(cells) > 50000) {
  cells <- sample(cells, 50000)
  message("▶ 采样50,000细胞加速计算")
}
ndufb7_sub <- ndufb7_expr[cells]

# 计算相关
mat <- tryCatch(LayerData(obj, assay="RNA", layer="data")[, cells], error=function(e) {
  tryCatch(GetAssayData(obj, assay="RNA", slot="data")[, cells], error=function(e2) NULL)
})
if(is.null(mat)) { message("❌ 无法提取矩阵"); quit(status=1) }

cor_genes <- apply(mat, 1, function(g) cor(g, ndufb7_sub, use="complete.obs"))
cor_df <- data.frame(Gene=names(cor_genes), Correlation=cor_genes, stringsAsFactors=FALSE)
cor_df <- cor_df[!is.na(cor_df$Correlation), ]
cor_df <- cor_df[order(-abs(cor_df$Correlation)), ]

# 标记预期KO效应
cor_df$KO_Expected <- ifelse(cor_df$Correlation > 0, "Down_regulated", "Up_regulated")
cor_df$KO_Expected[abs(cor_df$Correlation) < 0.1] <- "Unchanged"

message("▶ Top 10 正相关 (KO后预期下调):")
print(head(cor_df[cor_df$KO_Expected=="Down_regulated", c("Gene","Correlation")], 10))
message("▶ Top 10 负相关 (KO后预期上调):")
print(head(cor_df[cor_df$KO_Expected=="Up_regulated", c("Gene","Correlation")], 10))

# 铁死亡基因子集
ferro_genes <- c("FTL","FTH1","SLC7A11","GPX4","SAT1","ACSL4","NFE2L2","KEAP1","LPCAT3","PTGS2","ALOX15","GLS2")
ferro_df <- cor_df[cor_df$Gene %in% ferro_genes, ]
message("▶ 铁死亡基因共表达:")
print(ferro_df[, c("Gene","Correlation","KO_Expected")])

# 保存
out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_NDUFB7_Coexpression_KO.csv"
write.csv(cor_df, out, row.names=FALSE)
message("✅ 保存: ", out)

# 与铁死亡评分的相关性验证
if("Ferroptosis1" %in% colnames(obj@meta.data)) {
  message("▶ NDUFB7 vs 铁死亡评分相关: r=", round(cor(ndufb7_expr, obj$Ferroptosis1, use="complete.obs"), 3))
}
