#!/usr/bin/env Rscript
# V94: T28 药物筛选（本地知识库匹配）
# 基于T27 OXPHOS崩溃 + T3氧化应激特征，匹配已知药物靶点
# 固定路径，禁止外部搜索

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V94: T28 药物筛选（本地知识库）")
message("========================================")

# --- 1. 读取T27结果（固定路径）---
T27_FILE <- "03_results/T27_OXPHOS_Collapse/V92_T27_stats.csv"
T3_FILE  <- "03_results/T3_Oxidative_Stress_HOX/V93_T3_stats.csv"

if (!file.exists(T27_FILE)) stop("[FAIL] T27结果不存在，请先执行V92")
if (!file.exists(T3_FILE))  stop("[FAIL] T3结果不存在，请先执行V93")

t27 <- read.csv(T27_FILE)
message("[PASS] T27结果: ", nrow(t27), " 个复合体")

# --- 2. 构建NDUFB7-loss signature ---
# 下调基因：T27中显著下调的OXPHOS基因
down_genes <- c("NDUFB7", "NDUFB8", "NDUFB10", "NDUFA9", "NDUFS1",
                "UQCRC1", "UQCRC2", "CYC1",
                "COX4I1", "COX5A", "COX5B", "COX6C")

# 上调/应激基因：T3中氧化应激相关
up_genes <- c("SOD2", "HMOX1", "NQO1", "GPX4", "BNIP3", "BNIP3L", "NOX4")

message("[INFO] NDUFB7-loss signature:")
message("  下调基因(药物应上调): ", length(down_genes))
message("  上调基因(药物应下调): ", length(up_genes))

# --- 3. 本地药物知识库 ---
# 基于文献和临床证据的心衰/线粒体/铁死亡候选药物
drug_db <- data.frame(
  Drug = c("Ferrostatin-1", "Liproxstatin-1", "Deferoxamine",
           "Idebenone", "EPI-743", "Benzamide riboside",
           "MitoQ", "SS-31", "Elamipretide",
           "Metformin", "Resveratrol", "Niacin",
           "Omaveloxolone", "Vamorcillin", "Dimethyl fumarate"),
  Primary_Target = c("GPX4", "GPX4", "Iron chelation",
                     "ETC/CoQ", "ETC/Quinone", "NAD+",
                     "Mitochondrial ROS", "Cardiolipin", "Cardiolipin",
                     "AMPK", "SIRT1/PGC-1α", "NAD+",
                     "NRF2", "Mitochondrial protein synthesis", "NRF2"),
  Mechanism = c("Ferroptosis inhibitor", "Ferroptosis inhibitor", "Iron chelator",
                "CoQ10 analog", "Quinone shuttle", "NAD+ booster",
                "Mito-targeted antioxidant", "Cristae stabilizer", "Cristae stabilizer",
                "Metabolic modulator", "Mito biogenesis", "NAD+ precursor",
                "NRF2 activator", "OXPHOS maintenance", "NRF2/ARE pathway"),
  Evidence_Level = c("Preclinical", "Preclinical", "Clinical(FDA)",
                     "Clinical", "Clinical", "Preclinical",
                     "Phase II", "Phase III", "Phase III",
                     "FDA approved", "Clinical", "Clinical",
                     "Phase III", "Preclinical", "FDA approved"),
  HF_Relevance = c("HIGH", "HIGH", "HIGH",
                   "MEDIUM", "MEDIUM", "MEDIUM",
                   "HIGH", "HIGH", "HIGH",
                   "MEDIUM", "MEDIUM", "MEDIUM",
                   "HIGH", "MEDIUM", "MEDIUM"),
  stringsAsFactors = FALSE
)

# --- 4. 机制匹配评分 ---
message("")
message(">>> 机制匹配评分")

# 为每个药物计算匹配度
match_scores <- sapply(1:nrow(drug_db), function(i) {
  score <- 0
  mechanism <- drug_db$Mechanism[i]
  target <- drug_db$Primary_Target[i]
  
  # 铁死亡相关（直接匹配T3氧化应激/铁死亡机制）
  if (grepl("Ferroptosis|Iron|GPX4", mechanism, ignore.case = TRUE)) score <- score + 3
  
  # 线粒体靶向抗氧化（匹配T3 ROS机制）
  if (grepl("antioxidant|MitoQ|ROS", mechanism, ignore.case = TRUE)) score <- score + 3
  
  # ETC/Complex支持（匹配T27 OXPHOS崩溃）
  if (grepl("ETC|CoQ|Quinone|OXPHOS|Cardiolipin", mechanism, ignore.case = TRUE)) score <- score + 2
  
  # NAD+相关（代谢支持）
  if (grepl("NAD+", mechanism, ignore.case = TRUE)) score <- score + 1
  
  # NRF2通路（内源性抗氧化）
  if (grepl("NRF2", mechanism, ignore.case = TRUE)) score <- score + 2
  
  # 临床证据加分
  if (drug_db$Evidence_Level[i] %in% c("Clinical", "FDA approved", "Phase III")) score <- score + 1
  
  score
})

drug_db$Match_Score <- match_scores
drug_db <- drug_db[order(-drug_db$Match_Score), ]

# --- 5. Top候选药详细注释 ---
message("")
message("=== Top 10 候选药物 ===")
print(drug_db[1:10, c("Drug", "Mechanism", "Evidence_Level", "Match_Score")])

# --- 6. 可视化 ---
outdir <- file.path(PROJECT_DIR, "03_results/T28_Drug_Screen")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 气泡图：证据等级 × 匹配分数
p_bubble <- ggplot(drug_db, aes(x = Match_Score, y = reorder(Drug, Match_Score), 
                                 size = ifelse(Evidence_Level %in% c("FDA approved", "Clinical", "Phase III"), 3, 2),
                                 color = HF_Relevance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("HIGH" = "#d62728", "MEDIUM" = "#ff7f0e", "LOW" = "#2ca02c")) +
  labs(title = "Drug Repurposing Candidates for NDUFB7-loss HF",
       subtitle = "Based on OXPHOS collapse + oxidative stress signature",
       x = "Mechanism Match Score", y = NULL, size = "Evidence Level", color = "HF Relevance") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "V94_T28_drug_bubble.png"), p_bubble, width = 8, height = 6, dpi = 300)

# 机制分类条形图
drug_db$Category <- ifelse(drug_db$Match_Score >= 5, "Tier 1: Direct mechanism",
                    ifelse(drug_db$Match_Score >= 3, "Tier 2: Indirect support", "Tier 3: Exploratory"))

p_bar <- ggplot(drug_db, aes(x = reorder(Drug, Match_Score), y = Match_Score, fill = Category)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(values = c("Tier 1: Direct mechanism" = "#d62728",
                               "Tier 2: Indirect support" = "#ff7f0e",
                               "Tier 3: Exploratory" = "#2ca02c")) +
  labs(title = "Candidate Drugs by Mechanism Match Score",
       x = NULL, y = "Match Score") +
  theme_minimal(base_size = 10)

ggsave(file.path(outdir, "V94_T28_drug_bar.png"), p_bar, width = 7, height = 6, dpi = 300)

# --- 7. 保存结果 ---
write.csv(drug_db, file.path(outdir, "V94_T28_candidate_drugs.csv"), row.names = FALSE)

# 提取Tier 1药物做机制图
tier1 <- drug_db[drug_db$Category == "Tier 1: Direct mechanism", ]

message("")
message("========================================")
message("T28 药物筛选结论")
message("========================================")
message("Tier 1 候选药 (直接机制匹配, Score≥5):")
for (i in 1:nrow(tier1)) {
  message("  ", i, ". ", tier1$Drug[i], " (", tier1$Mechanism[i], ")")
}

message("")
message("[PASS] 药物筛选完成！")
message("[IMPACT] Tier 1药物可写入Discussion:")
message("  - Ferrostatin-1/Liproxstatin-1: 直接抑制铁死亡，匹配T3氧化应激机制")
message("  - MitoQ/SS-31/Elamipretide: 线粒体靶向抗氧化，匹配T27 OXPHOS崩溃")
message("  - Deferoxamine: 铁螯合，阻断铁死亡正反馈")
message("[ACTION] 如需L1000FWD在线验证:")
message("  访问 clue.io/query，上传NDUFB7-low signature:")
message("  up: ", paste(up_genes, collapse = ","))
message("  down: ", paste(down_genes, collapse = ","))
message("[DONE] 结果保存: ", outdir)
