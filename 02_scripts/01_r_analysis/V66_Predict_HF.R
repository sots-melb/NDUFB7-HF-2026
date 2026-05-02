suppressMessages(library(GEOquery))
suppressMessages(library(ggplot2))

matrix_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
message("▶ 正在加载并解析 GSE59867 临床特征...")
gse <- getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE, AnnotGPL = FALSE)
pdata <- pData(gse)
exprs_data <- exprs(gse)

# 提取探针
ndufb7_expr <- as.numeric(exprs_data["8034843", ])
df <- data.frame(Sample = rownames(pdata), NDUFB7 = ndufb7_expr)

# 动态特征提取函数
get_feature <- function(row, keyword) {
    for(val in row) {
        if(grepl(keyword, as.character(val), ignore.case=TRUE)) return(as.character(val))
    }
    return(NA)
}

# 提取时间点和是否进展为心衰
df$TimePoint <- apply(pdata, 1, function(x) get_feature(x, "samples collection"))
df$Progression <- apply(pdata, 1, function(x) get_feature(x, "hf progression"))

# 1. 过滤：只保留入院第一天 (Day 1 / admission) 的样本
df_day1 <- df[grepl("1st day|admission", df$TimePoint, ignore.case=TRUE), ]

# 2. 清洗：去除 "hf progression: " 前缀，保留 yes/no
df_day1$Progression <- gsub("hf progression: ", "", df_day1$Progression, ignore.case=TRUE)
df_day1 <- df_day1[df_day1$Progression %in% c("yes", "no"), ]
df_day1$Progression <- factor(df_day1$Progression, levels=c("no", "yes"), labels=c("Stable (No HF)", "HF Progression"))

message(sprintf("✅ 成功锁定 %d 名心梗患者的入院首日血液样本！", nrow(df_day1)))

# 3. 统计学检验 (Wilcoxon 秩和检验)
res <- wilcox.test(NDUFB7 ~ Progression, data = df_day1)

# 4. Logistic 回归计算风险比 (Odds Ratio)
df_day1$Prog_Bin <- ifelse(df_day1$Progression == "HF Progression", 1, 0)
# 将 NDUFB7 转换为 Z-score，计算每降低 1 个标准差带来的风险变化
df_day1$NDUFB7_Z <- scale(df_day1$NDUFB7)
log_model <- glm(Prog_Bin ~ NDUFB7_Z, data = df_day1, family = binomial)
or <- exp(-summary(log_model)$coefficients[2, "Estimate"]) # 负号代表表达量降低带来的风险倍数
p_log <- summary(log_model)$coefficients[2, "Pr(>|z|)"]

message("\n=======================================================")
message("🌟 临床转化闭环：入院首日 NDUFB7 预测心衰恶化结果 🌟")
message(sprintf(" 1. Wilcoxon 组间差异显著性: p = %.4f", res$p.value))
message(sprintf(" 2. Logistic 回归预测 p-value: p = %.4f", p_log))
message(sprintf(" 3. 临床风险: NDUFB7 每降低 1 个标准差，心衰恶化风险增加 %.2f 倍！", or))
message("=======================================================")

# 5. 绘制投稿级 Figure 6 临床预测箱线图
fig6 <- ggplot(df_day1, aes(x=Progression, y=NDUFB7, fill=Progression)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8, width=0.6, color="black", linewidth=1) +
  geom_jitter(width=0.2, size=3, color="black", alpha=0.6, shape=21, fill="white", stroke=1) +
  scale_fill_manual(values=c("Stable (No HF)"="#3C5488FF", "HF Progression"="#DC0000FF")) +
  theme_classic(base_size=16) +
  labs(title="Figure 6: Clinical Prognostic Biomarker",
       subtitle=sprintf("Low NDUFB7 at admission predicts HF progression (p=%.4f)", res$p.value),
       x="Clinical Outcome (Within 6 months post-MI)", 
       y="NDUFB7 Expression (PBMC Day 1)") +
  theme(legend.position="none", 
        plot.title=element_text(face="bold"),
        axis.text.x = element_text(face="bold", color="black"),
        axis.text.y = element_text(face="bold", color="black"))

pdf_out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/01_figures/Figure_6_Clinical_Prediction.pdf"
ggsave(pdf_out, fig6, width = 7, height = 6)
message(paste("✅ Figure 6 临床预测图已成功保存至:", pdf_out))
