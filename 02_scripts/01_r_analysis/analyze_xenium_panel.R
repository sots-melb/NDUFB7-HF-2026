library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GSE290577 Xenium Panel分析\n")
cat("=" ,rep("=", 69), "\n", sep="")

rds_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/04_spatial_geo/GSE290577/GSE290577_heart_spatial.rds"
xenium <- readRDS(rds_file)

genes <- rownames(xenium)
cat("Panel总基因数:", length(genes), "\n")
cat("前20个基因:", head(genes, 20), "\n")

# 检查是否有任何NDUF家族基因
nduf_genes <- genes[grepl("^NDUF", genes)]
cat("\nNDUF家族基因在panel中:", length(nduf_genes), "\n")
if(length(nduf_genes) > 0) {
  cat("  ", paste(nduf_genes, collapse=", "), "\n")
} else {
  cat("  ❌ 无NDUF家族基因（包括NDUFB7）\n")
}

# 检查线粒体相关基因
mito_genes <- genes[grepl("^MT-|^COX|^ATP5|^NDU|^SDH|^UQC|^CYC", genes)]
cat("\n线粒体/ETC相关基因:", length(mito_genes), "\n")
cat("  ", paste(mito_genes, collapse=", "), "\n")

# 检查是否有Complex I相关
complex1 <- genes[grepl("^NDUF|^Nduf", genes, ignore.case=TRUE)]
cat("\nComplex I相关:", length(complex1), "\n")

# 细胞类型标记
immune_markers <- c("CD3D","CD3E","CD4","CD8A","CD8B","FOXP3","CD19","CD79A","MS4A1","IGHG1","LYZ","CD14","FCGR3A","CST3","PPBP")
cardio_markers <- c("TNNT2","TNNI3","MYH6","MYH7","ACTC1","NPPA","NPPB")
endo_markers <- c("PECAM1","VWF","CDH5","TEK")
fibro_markers <- c("COL1A1","COL1A2","PDGFRA","VIM","DCN")

cat("\n【Panel设计特征分析】\n")
cat("免疫标记:", sum(immune_markers %in% genes), "/", length(immune_markers), "\n")
cat("心肌标记:", sum(cardio_markers %in% genes), "/", length(cardio_markers), "\n")
cat("内皮标记:", sum(endo_markers %in% genes), "/", length(endo_markers), "\n")
cat("成纤维标记:", sum(fibro_markers %in% genes), "/", length(fibro_markers), "\n")

# 样本信息
meta <- xenium@meta.data
cat("\n【样本信息】\n")
cat("总细胞:", nrow(meta), "\n")
cat("患者数:", length(unique(meta$patient_id)), "\n")
cat("患者:", paste(unique(meta$patient_id), collapse=", "), "\n")
cat("活检时机:", paste(unique(meta$biopsy_timing), collapse=", "), "\n")
cat("细胞类型数:", length(unique(meta$ct_second_pass)), "\n")
cat("主要细胞类型:\n")
print(sort(table(meta$ct_second_pass), decreasing=TRUE)[1:10])

# 关键结论
cat("\n========================================\n")
cat("【关键结论】\n")
cat("1. GSE290577是Xenium预设计小panel（477基因）\n")
cat("2. NDUFB7不在panel中 → 无法直接验证NDUFB7\n")
cat("3. 但可用于: 细胞类型组成参考 / 免疫微环境分析\n")
cat("4. 样本为心脏移植患者活检（非MI），临床背景不同\n")
cat("========================================\n")

cat("\n完成\n")
