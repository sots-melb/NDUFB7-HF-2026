library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GDS4772: 修复NF组NA\n")
cat("=" ,rep("=", 69), "\n", sep="")

ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772_fixed.txt")

# 问题: val_cols包含"IDENTIFIER"（非样本列）
val_cols <- setdiff(names(ndufb7), c("ID_REF", " IDENTIFIER "))
cat("原始列:", val_cols, "\n")

# 排除IDENTIFIER（它是基因注释列，不是样本）
# 实际样本列是GSM开头的
sample_cols <- val_cols[grepl("^GSM", trimws(val_cols))]
cat("样本列:", length(sample_cols), "\n")
cat("  ", sample_cols, "\n")

vals <- as.numeric(unlist(ndufb7[, ..sample_cols]))
vals <- vals[!is.na(vals)]
cat("有效表达值:", length(vals), "\n")

# 重新解析分组
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/GDS4772.soft.gz"
lines <- readLines(gzfile(f))

subset_desc <- lines[grep("^!subset_description", lines)]
subset_sample <- lines[grep("^!subset_sample_id", lines)]

dcm_line <- subset_sample[grep("dilated", subset_desc)]
dcm_samples <- trimws(strsplit(gsub('!subset_sample_id = ', '', dcm_line), ",")[[1]])

nf_line <- subset_sample[grep("normal", subset_desc)]
nf_samples <- trimws(strsplit(gsub('!subset_sample_id = ', '', nf_line), ",")[[1]])

cat("\nDCM样本:", length(dcm_samples), "\n")
cat("NF样本:", length(nf_samples), "\n")

# 精确匹配
clean_cols <- trimws(sample_cols)
dcm_mask <- clean_cols %in% dcm_samples
nf_mask <- clean_cols %in% nf_samples

cat("DCM匹配:", sum(dcm_mask), "\n")
cat("NF匹配:", sum(nf_mask), "\n")

dcm_vals <- vals[dcm_mask]
nf_vals <- vals[nf_mask]

cat("\n【修正后GDS4772分层】\n")
cat("DCM (n=", length(dcm_vals), "): 均值=", round(mean(dcm_vals), 4), ", 中位数=", round(median(dcm_vals), 4), "\n")
cat("NF (n=", length(nf_vals), "): 均值=", round(mean(nf_vals), 4), ", 中位数=", round(median(nf_vals), 4), "\n")

if(length(dcm_vals) > 0 && length(nf_vals) > 0) {
  wt <- wilcox.test(dcm_vals, nf_vals)
  cat("\nWilcoxon DCM vs NF: p=", format(wt$p.value, digits=4), "\n")
  
  pooled_sd <- sqrt(((length(dcm_vals)-1)*sd(dcm_vals)^2 + (length(nf_vals)-1)*sd(nf_vals)^2) / 
                    (length(dcm_vals) + length(nf_vals) - 2))
  d <- (mean(dcm_vals) - mean(nf_vals)) / pooled_sd
  cat("Cohen's d:", round(d, 4), "\n")
}

# 保存修正版
result <- data.table(
  Group = c(rep("DCM", length(dcm_vals)), rep("NF", length(nf_vals))),
  Expression = c(dcm_vals, nf_vals),
  Sample = c(clean_cols[dcm_mask], clean_cols[nf_mask])
)
fwrite(result, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_stratified_fixed.txt", sep="\t")
cat("\n✅ 修正版分层数据已保存\n")

cat("\n完成\n")
