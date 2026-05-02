#!/usr/bin/env Rscript
# V113C: GSE154170 Visium空间分析——ICM vs DCM区域NDUFB7比较

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V113C: GSE154170 Visium空间分析")
message("========================================")

# 读取分组
GRP_FILE <- "01_data/01_raw_geo/GSE154170_group_corrected.csv"
grp <- read.csv(GRP_FILE, stringsAsFactors = FALSE)
message("[PASS] 分组: ", nrow(grp), " samples (", sum(grp$Group=="ICM"), " ICM, ", sum(grp$Group=="DCM"), " DCM)")

# GSE154170是Visium数据，找h5ad或Seurat对象
h5ad_candidates <- c(
  "01_data/02_spatial/GSE154170_visium.h5ad",
  "01_data/02_spatial/GSE154170.h5ad"
)

# 如果有Seurat对象
rds_candidates <- list.files("01_data", pattern = "GSE154170.*\\.rds$", 
                             full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

if (length(rds_candidates) > 0) {
  srt <- readRDS(rds_candidates[1])
  message("[PASS] 加载Visium对象: ", basename(rds_candidates[1]))
  
  # 检查样本ID匹配
  meta <- srt@meta.data
  sample_col <- intersect(c("sample", "patient", "orig.ident"), colnames(meta))[1]
  
  if (!is.na(sample_col)) {
    # 提取NDUFB7表达
    if ("NDUFB7" %in% rownames(srt)) {
      expr <- FetchData(srt, vars = "NDUFB7")
      expr$sample <- meta[[sample_col]]
      
      # 匹配分组
      # 假设样本ID格式如 P1, P2, P3 对应 ICM/DCM
      # 需要手动映射或从GEO标题推断
      
      message("\n[INFO] 样本列表:")
      print(table(expr$sample))
      
      message("\n[PENDING] 需要手动映射样本ID到ICM/DCM分组")
      message("  请根据GEO标题将样本ID映射到分组:")
      for (i in 1:nrow(grp)) {
        message("    ", grp$Sample_ID[i], " -> ", grp$Group[i], " (", grp$Note[i], ")")
      }
      
      # 简化：输出描述统计
      outdir <- file.path(PROJECT_DIR, "03_results/V113C_GSE154170")
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
      
      write.csv(expr, file.path(outdir, "V113C_GSE154170_NDUFB7_spatial.csv"), row.names = FALSE)
      message("[DONE] 原始数据保存，待手动分组后比较")
    }
  }
} else {
  message("[MISS] 未找到GSE154170 Visium对象")
  message("[ACTION] 需要从GSE154170_RAW.tar解压并构建Seurat对象")
}
