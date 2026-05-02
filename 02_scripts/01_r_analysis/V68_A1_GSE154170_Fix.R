suppressMessages(library(dplyr))
csv_paths <- c(
  "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_figures_tables/GSE154170_NDUFB7_expression.csv",
  "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE154170_NDUFB7_expression.csv",
  "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE154170_Resurrected_Group.csv"
)
csv <- NULL; for(p in csv_paths) if(file.exists(p)) { csv <- p; break }
if(is.null(csv)) { message("❌ 未找到GSE154170文件"); quit(status=1) }
message("▶ 使用: ", csv)
df <- read.csv(csv, stringsAsFactors=FALSE)
message("▶ 列名: ", paste(colnames(df), collapse=", "))
samp_col <- colnames(df)[1]
df$Group <- ifelse(grepl("^RSB", df[[samp_col]], ignore.case=TRUE), "Control",
            ifelse(grepl("^RSH", df[[samp_col]], ignore.case=TRUE), "ICM",
            ifelse(grepl("^RSN", df[[samp_col]], ignore.case=TRUE), "Normal", "Unknown")))
message("▶ 分组分布:"); print(table(df$Group, useNA="ifany"))
ndufb7_col <- grep("NDUFB7", colnames(df), value=TRUE)[1]
if(!is.na(ndufb7_col)) {
  message("▶ 各组NDUFB7统计:")
  print(aggregate(as.formula(paste(ndufb7_col,"~Group")), data=df, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x),N=length(x))))
  df_icm <- df[df$Group %in% c("ICM","Control"), ]
  if(nrow(df_icm)>0 && length(unique(df_icm$Group))==2) {
    wc <- wilcox.test(as.formula(paste(ndufb7_col,"~Group")), data=df_icm)
    message(sprintf("▶ ICM vs Control: W=%.1f, p=%.4f", wc$statistic, wc$p.value))
  }
  out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE154170_ICM_Final.csv"
  write.csv(df, out, row.names=FALSE)
  message("✅ 保存: ", out)
} else { message("❌ 无NDUFB7列") }
