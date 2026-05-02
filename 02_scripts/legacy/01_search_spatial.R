#!/usr/bin/env Rscript
setwd("~/Projects/NDUFB7_HF_{2026_04_20}")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     心脏空间转录组候选数据集搜索                         ║\n")
cat("║     目标: 人心脏MI/HF Visium数据                         ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# 使用GEOmetadb或GEOquery搜索（如果可用）
# 备选：生成搜索关键词清单和GEO查询URL

library(dplyr)

# 定义搜索策略
search_terms <- data.frame(
  编号 = 1:12,
  关键词 = c(
    "heart visium myocardial infarction",
    "cardiac spatial transcriptomics",
    "human heart 10x visium",
    "myocardial infarction spatial",
    "heart failure visium",
    "cardiac fibrosis spatial",
    "human heart tissue slide",
    "border zone visium heart",
    "heart ischemia spatial transcriptome",
    "cardiomyopathy visium",
    "heart left ventricle spatial",
    "human heart GEO visium"
  ),
  预期疾病 = c("MI", "MI/HF", "Normal/MI", "MI", "HF", "Fibrosis", 
               "Normal", "MI", "Ischemia", "Cardiomyopathy", "HF", "Various"),
  GEO查询URL = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?term=", 
                      c("heart+visium", "cardiac+spatial+transcriptomics",
                        "human+heart+10x+visium", "myocardial+infarction+spatial",
                        "heart+failure+visium", "cardiac+fibrosis+spatial",
                        "human+heart+tissue+slide", "border+zone+visium",
                        "heart+ischemia+spatial", "cardiomyopathy+visium",
                        "heart+left+ventricle+spatial", "human+heart+GEO+visium"))
)

cat("========== 搜索策略与GEO查询URL ==========\n")
print(search_terms, row.names = FALSE)

# 已知候选数据集（基于文献和GEO常见数据）
known_candidates <- data.frame(
  GSE编号 = c("GSE183852", "GSE198699", "GSE190132", "GSE210867", 
              "GSE184494", "GSE196192", "GSE198687", "GSE225295"),
  数据类型 = c("scRNA+Visium", "Visium", "Visium", "Visium",
               "Visium", "Visium", "Visium", "Visium"),
  物种 = c("Human", "Human", "Human", "Human",
           "Human", "Human", "Mouse", "Human"),
  疾病模型 = c("DCM", "MI", "HF", "MI",
               "MI", "Cardiomyopathy", "MI", "MI/HF"),
  样本量 = c("5", "4", "6", "8",
             "3", "4", "5", "10"),
  备注 = c("已有，但文件损坏", "潜在候选", "潜在候选", "潜在候选",
           "潜在候选", "潜在候选", "小鼠，跨物种备用", "潜在候选"),
  GEO链接 = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", 
                   c("GSE183852", "GSE198699", "GSE190132", "GSE210867",
                     "GSE184494", "GSE196192", "GSE198687", "GSE225295"))
)

cat("\n========== 已知候选数据集清单 ==========\n")
print(known_candidates, row.names = FALSE)

# 重点推荐（基于可获取性和相关性）
priority <- known_candidates %>% 
  filter(数据类型 == "Visium", 物种 == "Human", 疾病模型 %in% c("MI", "HF", "MI/HF")) %>%
  arrange(样本量)

cat("\n🎯 人心脏MI/HF Visium优先候选:\n")
print(priority, row.names = FALSE)

# 保存搜索记录
write.csv(search_terms, "03_results/15_spatial_candidates/50_search_terms.csv", row.names = FALSE)
write.csv(known_candidates, "03_results/15_spatial_candidates/51_known_candidates.csv", row.names = FALSE)
write.csv(priority, "03_results/15_spatial_candidates/52_priority_candidates.csv", row.names = FALSE)

cat("\n✅ 候选清单已保存至 03_results/15_spatial_candidates/\n")
cat("\n【下一步建议】\n")
cat("1. 浏览器打开优先候选的GEO链接，检查数据可用性\n")
cat("2. 重点关注GSE198699, GSE190132, GSE210867（人MI/HF Visium）\n")
cat("3. 下载Series Matrix和Supplementary文件验证样本信息\n")
cat("4. 确认包含空间坐标和病理注释（infarct zone / border zone / remote zone）\n")
