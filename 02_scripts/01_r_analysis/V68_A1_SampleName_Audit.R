suppressMessages(library(dplyr))
csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_figures_tables/GSE154170_NDUFB7_expression.csv"
df <- read.csv(csv, stringsAsFactors=FALSE)
message("▶ 总样本数: ", nrow(df))
message("▶ 列名: ", paste(colnames(df), collapse=" | "))
message("▶ 全部样本名 (前16个):")
print(df[[1]])
message("▶ 样本名唯一值:")
print(unique(df[[1]]))
# 尝试自动推断分组
df$Group <- "Unknown"
df$Group[grepl("normal|healthy|control|ncm|non", df[[1]], ignore.case=TRUE)] <- "Control"
df$Group[grepl("icm|isch|infarct|mi", df[[1]], ignore.case=TRUE)] <- "ICM"
df$Group[grepl("dcm|dilated|cardio", df[[1]], ignore.case=TRUE)] <- "DCM"
df$Group[grepl("hf|fail|heart", df[[1]], ignore.case=TRUE)] <- "HF"
message("▶ 自动推断分组:")
print(table(df$Group, useNA="ifany"))
# 保存审计
out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE154170_SampleAudit.csv"
write.csv(df, out, row.names=FALSE)
message("✅ 审计保存: ", out)
