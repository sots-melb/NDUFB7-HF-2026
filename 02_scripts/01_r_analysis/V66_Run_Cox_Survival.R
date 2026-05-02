# --- 自动安装生存分析依赖包 ---
for (pkg in c("survival", "survminer", "GEOquery")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos="http://cran.us.r-project.org")
}
suppressMessages(library(GEOquery))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(dplyr))

message("▶ 正在加载 GSE59867 Series Matrix...")
matrix_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"

if(file.exists(matrix_file)) {
    gse <- getGEO(filename = matrix_file, GSEMatrix = TRUE)
} else {
    message("本地未找到，正在从 GEO 极速拉取...")
    gse <- getGEO("GSE59867", GSEMatrix = TRUE)[[1]]
}

pdata <- pData(gse)
exprs_data <- exprs(gse)

# --- 1. 提取 NDUFB7 表达量 (探针: 8034843) ---
probe_id <- "8034843"
if(probe_id %in% rownames(exprs_data)) {
    ndufb7_expr <- as.numeric(exprs_data[probe_id, ])
    message("✅ 成功提取 NDUFB7 (探针 8034843) 表达量")
} else {
    stop("❌ 在矩阵中未找到探针 8034843！")
}

# --- 2. 智能提取临床生存变量 ---
message("▶ 正在解析临床随访变量 (Time & Event)...")
# 在 GEO 数据中，特征通常藏在 characteristics_ch1 等列中
# 寻找包含 survival, time, event, death, follow-up 等字眼的列
clinical_df <- data.frame(Sample = rownames(pdata), NDUFB7 = ndufb7_expr)

# 暴力提取所有特征列拼接成大字符串进行正则匹配
feature_cols <- grep("characteristics_ch1", colnames(pdata), value=TRUE)
if(length(feature_cols) == 0) feature_cols <- colnames(pdata)

# 尝试提取 Time (以月或天为单位)
time_vec <- rep(NA, nrow(pdata))
event_vec <- rep(NA, nrow(pdata))

for(i in 1:nrow(pdata)) {
    row_text <- paste(pdata[i, feature_cols], collapse=" | ")
    
    # 提取时间 (提取 follow up, time 等字眼后的数字)
    time_match <- regexpr("(follow[- ]?up time|survival time).*?([0-9.]+)", row_text, ignore.case=TRUE)
    if(time_match > 0) {
        num_str <- regmatches(row_text, regexpr("([0-9.]+)", substr(row_text, time_match, nchar(row_text))))
        time_vec[i] <- as.numeric(num_str)
    }
    
    # 提取事件 (0=存活/截尾, 1=死亡/事件)
    if(grepl("event.*?1|death.*?1|status.*?dead", row_text, ignore.case=TRUE)) {
        event_vec[i] <- 1
    } else if(grepl("event.*?0|death.*?0|status.*?alive", row_text, ignore.case=TRUE)) {
        event_vec[i] <- 0
    }
}

clinical_df$Time <- time_vec
clinical_df$Event <- event_vec

# 如果智能提取失败，打印出前 5 个样本的 metadata 供人工诊断
if(all(is.na(clinical_df$Time)) || all(is.na(clinical_df$Event))) {
    message("⚠️ 警告：无法自动识别生存时间和终点事件。请查看以下 Metadata 结构：")
    print(head(pdata[, 1:min(10, ncol(pdata))], 3))
    stop("临床变量提取失败，分析中止。")
}

# 过滤掉缺失随访数据的样本
clinical_df <- clinical_df[!is.na(clinical_df$Time) & !is.na(clinical_df$Event), ]
message(sprintf("✅ 成功提取到 %d 名患者的有效随访数据", nrow(clinical_df)))

# --- 3. 将表达量二分类 (High vs Low) ---
median_val <- median(clinical_df$NDUFB7, na.rm=TRUE)
clinical_df$Group <- ifelse(clinical_df$NDUFB7 >= median_val, "NDUFB7 High", "NDUFB7 Low")
clinical_df$Group <- factor(clinical_df$Group, levels = c("NDUFB7 High", "NDUFB7 Low")) # High作为参照组

# --- 4. 运行 Cox 比例风险模型 ---
message("▶ 正在运行 Cox 比例风险回归模型...")
cox_model <- coxph(Surv(Time, Event) ~ NDUFB7, data = clinical_df)
cox_summary <- summary(cox_model)
hr <- cox_summary$coefficients[1, "exp(coef)"]
p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]

message("\n=======================================================")
message("🌟 临床预后闭环：GSE59867 Cox 回归结果 🌟")
message(sprintf(" 1. 风险比 (Hazard Ratio, HR): %.4f", hr))
message(sprintf(" 2. 统计显著性 (p-value):      %.4e", p_val))
if(hr < 1 && p_val < 0.05) {
    message("🎉 结论完美！NDUFB7 表达量越高，死亡风险越低 (保护性因素)！")
} else if (hr > 1 && p_val < 0.05) {
    message("⚠️ 结论显著但方向相反：NDUFB7 表达越高，死亡风险越高。")
} else {
    message("⚠️ 预后趋势不显著，可能需要调整分组阈值或纳入多因素 (年龄、性别) 校正。")
}
message("=======================================================")

# --- 5. 绘制 Kaplan-Meier 生存曲线 ---
fit <- survfit(Surv(Time, Event) ~ Group, data = clinical_df)
km_plot <- ggsurvplot(
    fit, 
    data = clinical_df,
    pval = TRUE,             # 显示 Log-rank p值
    conf.int = TRUE,         # 显示置信区间
    risk.table = TRUE,       # 底部显示风险表
    palette = c("#3C5488FF", "#DC0000FF"), # 顶刊配色：高表达蓝色，低表达红色
    title = "Kaplan-Meier Survival Analysis (GSE59867)",
    xlab = "Follow-up Time",
    ylab = "Survival Probability",
    legend.title = "Expression",
    ggtheme = theme_classic(base_size = 15)
)

pdf_out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/01_figures/Figure_6_KM_Survival.pdf"
pdf(pdf_out, width = 8, height = 7, onefile = FALSE)
print(km_plot)
dev.off()

message(paste("✅ Kaplan-Meier 生存曲线已成功保存至:", pdf_out))
