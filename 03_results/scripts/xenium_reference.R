library(data.table)

cat("========================================\n")
cat("GSE290577 Xenium: 心脏细胞类型参考\n")
cat("========================================\n")

rds_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/04_spatial_geo/GSE290577/GSE290577_heart_spatial.rds"

if(!file.exists(rds_file)) {
  cat("❌ RDS文件不存在\n")
  quit(status=1)
}

xenium <- readRDS(rds_file)
meta <- xenium@meta.data

# 细胞类型统计
ct_table <- sort(table(meta$ct_second_pass), decreasing=TRUE)
cat("\n【细胞类型组成】\n")
print(ct_table)

# 百分比
total_cells <- sum(ct_table)
cat("\n【百分比】\n")
for(i in 1:length(ct_table)) {
  pct <- round(ct_table[i] / total_cells * 100, 1)
  cat(names(ct_table)[i], ":", pct, "% (n=", ct_table[i], ")\n")
}

# 患者信息
cat("\n【患者信息】\n")
cat("患者数:", length(unique(meta$patient_id)), "\n")
cat("患者ID:", paste(unique(meta$patient_id), collapse=", "), "\n")
cat("活检时机:", paste(unique(meta$biopsy_timing), collapse=", "), "\n")
cat("总细胞:", format(nrow(meta), big.mark=","), "\n")

# 关键结论
cardio_pct <- round(ct_table["Cardiomyocytes"] / total_cells * 100, 1)
fibro_pct <- round(ct_table["Fibroblasts"] / total_cells * 100, 1)
immune_pct <- round(sum(ct_table[grepl("T cells|Macrophages|NK cells|DC|B cells|Plasma", names(ct_table))]) / total_cells * 100, 1)

cat("\n========================================\n")
cat("【关键结论】\n")
cat("1. 心肌细胞占比:", cardio_pct, "% — NDUFB7主要表达细胞类型\n")
cat("2. 成纤维细胞占比:", fibro_pct, "% — 纤维化区主要细胞\n")
cat("3. 免疫细胞总占比:", immune_pct, "% — 炎症微环境\n")
cat("4. Panel限制: 477预设计基因，无NDUFB7/线粒体基因\n")
cat("5. 临床背景: 心脏移植患者（非MI），但与Kuppe互补\n")
cat("========================================\n")

# 保存参考数据
ref <- data.table(
  CellType = names(ct_table),
  Count = as.integer(ct_table),
  Percentage = round(as.numeric(ct_table) / total_cells * 100, 2)
)
fwrite(ref, "~/Projects/NDUFB7_HF_2026_04_20/03_results/xenium_cell_type_reference.txt", sep="\t")
cat("\n✅ 已保存: xenium_cell_type_reference.txt\n")

cat("\n完成\n")
