suppressMessages(library(GEOquery))
mat <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
message("▶ 加载 GSE59867...")
gse <- getGEO(filename=mat, GSEMatrix=TRUE, getGPL=FALSE, AnnotGPL=FALSE)
pdata <- pData(gse); exprs_data <- exprs(gse)
ndufb7_expr <- as.numeric(exprs_data["8034843", ])
df <- data.frame(Sample=rownames(pdata), NDUFB7=ndufb7_expr, stringsAsFactors=FALSE)
get_feature <- function(row, kw) {
  for(v in row) if(grepl(kw, as.character(v), ignore.case=TRUE)) return(as.character(v))
  return(NA)
}
df$TimePoint <- apply(pdata, 1, function(x) get_feature(x, "samples collection"))
df$Progression <- apply(pdata, 1, function(x) get_feature(x, "hf progression"))
df_day1 <- df[grepl("1st day|admission", df$TimePoint, ignore.case=TRUE), ]
df_day1$Progression <- tolower(trimws(gsub("hf progression:\\s*", "", df_day1$Progression, ignore.case=TRUE)))
message("▶ Day1 原始标签分布:")
print(table(df_day1$Progression, useNA="ifany"))
df_day1$Group <- ifelse(df_day1$Progression=="hf", "HF_Progression",
                 ifelse(df_day1$Progression=="non-hf", "No_HF_Progression",
                 ifelse(df_day1$Progression=="n/a", "Healthy_Control", NA)))
df_valid <- df_day1[!is.na(df_day1$Group), ]
df_valid$Group <- factor(df_valid$Group, levels=c("Healthy_Control","No_HF_Progression","HF_Progression"))
message("▶ 有效样本分组 (N=有效患者):")
print(table(df_valid$Group))
if(length(unique(df_valid$Group)) >= 2) {
  message("▶ 三组NDUFB7描述统计:")
  print(aggregate(NDUFB7~Group, data=df_valid, FUN=function(x) c(Mean=mean(x), SD=sd(x), Median=median(x), N=length(x))))
  kw <- kruskal.test(NDUFB7~Group, data=df_valid)
  message(sprintf("▶ Kruskal-Wallis: χ²=%.3f, p=%.4f", kw$statistic, kw$p.value))
  df_hf <- df_valid[df_valid$Group %in% c("HF_Progression","No_HF_Progression"), ]
  if(nrow(df_hf)>0 && length(unique(df_hf$Group))==2) {
    wc <- wilcox.test(NDUFB7~Group, data=df_hf)
    message(sprintf("▶ HF vs Non-HF (Wilcoxon): W=%.1f, p=%.4f", wc$statistic, wc$p.value))
  }
} else { message("❌ 分组不足") }
out_csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE59867_Day1_Descriptive.csv"
write.csv(df_valid, out_csv, row.names=FALSE)
message("✅ 保存: ", out_csv)
