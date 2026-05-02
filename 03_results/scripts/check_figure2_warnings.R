library(data.table)

cat("诊断Figure 2 warning来源...\n")

# 检查各数据集NA
g1 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
g2 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_stratified_fixed.txt")
g3 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_stratified.txt")
g4 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_stratified.txt")

cat("GSE57338 NA:", sum(is.na(as.numeric(unlist(g1[,-1])))), "\n")
cat("GDS4772 NA:", sum(is.na(g2$Expression)), "\n")
cat("GSE116250 NA:", sum(is.na(g3$Expression)), "\n")
cat("GSE55296 NA:", sum(is.na(g4$Expression)), "\n")

cat("✅ 诊断完成\n")
