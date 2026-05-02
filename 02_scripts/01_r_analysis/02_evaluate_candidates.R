#!/usr/bin/env Rscript
setwd("~/Projects/NDUFB7_HF_{2026_04_20}")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     空间转录组候选数据集评估报告                         ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

library(dplyr)

# 读取候选清单
candidates <- read.csv("03_results/15_spatial_candidates/52_priority_candidates.csv")

# 评估维度
eval_report <- data.frame(
  GSE = candidates$GSE编号,
  疾病 = candidates$疾病模型,
  样本量 = candidates$样本量,
  优先级 = candidates$备注,
  
  # 用户需要手动填写的评估项（浏览器访问GEO后填写）
  数据可用性 = c("待检查", "待检查", "待检查", "待检查", "待检查"),
  文件格式 = c("待确认", "待确认", "待确认", "待确认", "待确认"),
  病理注释 = c("待确认", "待确认", "待确认", "待确认", "待确认"),
  是否含NDUFB7 = c("待确认", "待确认", "待确认", "待确认", "待确认"),
  推荐度 = c("★★★★★", "★★★★☆", "★★★★☆", "★★★☆☆", "★★★☆☆"),
  
  GEO页面 = candidates$GEO链接,
  
  检查要点 = c(
    "1. 是否有.tar.gz/Seurat对象？ 2. 是否含infarct/border/remote注释？ 3. 样本是否配对？",
    "1. 时间点？ 2. 是否有sham对照？ 3. Visium分辨率？",
    "1. 病因类型（缺血性/非缺血性）？ 2. NYHA分级？ 3. 是否含纤维化区域注释？",
    "1. 急性/慢性MI？ 2. 是否含远端区域对照？",
    "1. 样本量较小，是否可合并？ 2. 数据质量？"
  )
)

cat("========== 优先候选评估表 ==========\n")
print(eval_report, row.names = FALSE)

# 生成浏览器检查清单
cat("\n========== 浏览器检查步骤 ==========\n")
cat("【Step 1】打开首选GSE225295:\n")
cat("  ", eval_report$GEO页面[1], "\n")
cat("  → 查看Summary确认疾病类型\n")
cat("  → 点击'Samples'查看样本列表和注释\n")
cat("  → 点击'SUPPLEMENTARY'查看文件列表\n")
cat("  → 确认是否有: spatial coordinates, tissue image, count matrix\n\n")

cat("【Step 2】检查文件大小和格式:\n")
cat("  → 理想的supplementary文件: *_filtered_feature_bc_matrix.h5, spatial/tissue_positions_list.csv\n")
cat("  → 或者: Seurat对象(.rds/.Robj), SpaceRanger输出目录\n\n")

cat("【Step 3】确认NDUFB7检测:\n")
cat("  → 查看Platform是否为10x Genomics Visium (GPLxxx)\n")
cat("  → 确认基因注释包含Symbol (非仅Ensembl)\n\n")

# 保存评估模板
write.csv(eval_report, "03_results/15_spatial_candidates/53_evaluation_template.csv", row.names = FALSE)

cat("✅ 评估模板已保存: 03_results/15_spatial_candidates/53_evaluation_template.csv\n")
cat("\n【建议】现在用浏览器打开前3个GEO链接，10分钟内完成评估，\n")
cat("      将结果填入CSV的'数据可用性'、'文件格式'、'病理注释'列\n")
