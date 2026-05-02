library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GDS4772: 精确病因分层\n")
cat("=" ,rep("=", 69), "\n", sep="")

# 读取NDUFB7表达
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772_fixed.txt")
val_cols <- setdiff(names(ndufb7), c("ID_REF", " IDENTIFIER "))

# 从GDS SOFT精确解析样本分组
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/GDS4772.soft.gz"
lines <- readLines(gzfile(f))

# 提取subset样本ID列表
subset_desc <- lines[grep("^!subset_description", lines)]
subset_sample <- lines[grep("^!subset_sample_id", lines)]

# 解析DCM样本
dcm_line <- subset_sample[grep("dilated", subset_desc)]
dcm_samples <- strsplit(gsub('!subset_sample_id = ', '', dcm_line), ",")[[1]]
dcm_samples <- trimws(dcm_samples)
cat("DCM样本:", length(dcm_samples), "\n")
cat("  ", paste(head(dcm_samples, 5), collapse=", "), "...\n")

# 解析Normal样本
nf_line <- subset_sample[grep("normal", subset_desc)]
nf_samples <- strsplit(gsub('!subset_sample_id = ', '', nf_line), ",")[[1]]
nf_samples <- trimws(nf_samples)
cat("Normal样本:", length(nf_samples), "\n")
cat("  ", paste(head(nf_samples, 5), collapse=", "), "...\n")

# 匹配样本ID到表达值
# val_cols格式: "GSM1053915", "GSM1053917"等（可能带空格）
clean_cols <- trimws(val_cols)

dcm_mask <- clean_cols %in% dcm_samples
nf_mask <- clean_cols %in% nf_samples

cat("\n匹配结果:\n")
cat("  DCM匹配:", sum(dcm_mask), "\n")
cat("  NF匹配:", sum(nf_mask), "\n")

# 提取表达值
vals <- as.numeric(unlist(ndufb7[, ..val_cols]))
vals <- vals[!is.na(vals)]

dcm_vals <- vals[dcm_mask]
nf_vals <- vals[nf_mask]

cat("\n【GDS4772病因分层】\n")
cat("DCM (n=", length(dcm_vals), "): 均值=", round(mean(dcm_vals), 4), ", 中位数=", round(median(dcm_vals), 4), "\n")
cat("NF (n=", length(nf_vals), "): 均值=", round(mean(nf_vals), 4), ", 中位数=", round(median(nf_vals), 4), "\n")

# 统计检验
if(length(dcm_vals) > 0 && length(nf_vals) > 0) {
  wt <- wilcox.test(dcm_vals, nf_vals)
  cat("\nWilcoxon DCM vs NF: p=", format(wt$p.value, digits=4), "\n")
  
  # 效应量
  pooled_sd <- sqrt(((length(dcm_vals)-1)*sd(dcm_vals)^2 + (length(nf_vals)-1)*sd(nf_vals)^2) / 
                    (length(dcm_vals) + length(nf_vals) - 2))
  d <- (mean(dcm_vals) - mean(nf_vals)) / pooled_sd
  cat("Cohen's d:", round(d, 4), "\n")
}

# 保存
result <- data.table(
  Group = c(rep("DCM", length(dcm_vals)), rep("NF", length(nf_vals))),
  Expression = c(dcm_vals, nf_vals),
  Sample = c(clean_cols[dcm_mask], clean_cols[nf_mask])
)
fwrite(result, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_stratified.txt", sep="\t")
cat("\n✅ 分层数据已保存\n")

cat("\n完成\n")
