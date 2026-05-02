message("========================================")
message("V69 关键数字提取")
message("========================================")

# 1. GSE183852 Condition统计
if(file.exists("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Condition_Stats.csv")) {
  df <- read.csv("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Condition_Stats.csv", stringsAsFactors=FALSE)
  message("\n[GSE183852]")
  message("  总细胞: ", nrow(df))
  message("  DCM: ", sum(df$Condition=="DCM", na.rm=TRUE), " | Donor: ", sum(df$Condition=="Donor", na.rm=TRUE))
  message("  CM_Score DCM vs Donor: 待Wilcoxon")
  message("  OXPHOS DCM vs Donor: 待Wilcoxon")
  message("  NDUFB7 DCM vs Donor: 待Wilcoxon")
}

# 2. GSE214611人类STEMI
if(file.exists("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE214611_Human_STEMI_NDUFB7.csv")) {
  df <- read.csv("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE214611_Human_STEMI_NDUFB7.csv", stringsAsFactors=FALSE)
  message("\n[GSE214611人类STEMI]")
  message("  总细胞: ", nrow(df))
  message("  NDUFB7>0比例: ", round(mean(df$NDUFB7>0), 3))
  message("  NDUFB7 Mean: ", round(mean(df$NDUFB7), 3))
}

# 3. GSE59867
if(file.exists("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE59867_Day1_Descriptive.csv")) {
  df <- read.csv("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE59867_Day1_Descriptive.csv", stringsAsFactors=FALSE)
  message("\n[GSE59867]")
  message("  总样本: ", nrow(df))
  print(aggregate(NDUFB7~Group, data=df, FUN=function(x) c(Mean=mean(x), SD=sd(x), N=length(x))))
}

# 4. KO替代
if(file.exists("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_NDUFB7_Coexpression_KO.csv")) {
  df <- read.csv("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_NDUFB7_Coexpression_KO.csv", stringsAsFactors=FALSE)
  message("\n[KO替代]")
  message("  总基因: ", nrow(df))
  message("  下调: ", sum(df$KO_Expected=="Down_regulated"))
  message("  上调: ", sum(df$KO_Expected=="Up_regulated"))
  top5 <- head(df[df$KO_Expected=="Down_regulated", c("Gene","Correlation")], 5)
  message("  Top5下调: ", paste(top5$Gene, collapse=", "))
}

# 5. WGCNA
wgcna_rds <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_WGCNA_Modules_Fix.RDS"
if(file.exists(wgcna_rds)) {
  net <- readRDS(wgcna_rds)
  message("\n[WGCNA]")
  message("  模块数: ", length(unique(net$colors)))
  if("NDUFB7" %in% names(net$colors)) {
    message("  NDUFB7模块: ", net$colors["NDUFB7"])
  }
} else {
  message("\n[WGCNA] ❌ 无产出")
}

# 6. SMR
smr_files <- list.files("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables", pattern="SMR_.*\\.smr", full.names=TRUE)
if(length(smr_files)>0) {
  message("\n[SMR]")
  for(f in smr_files) {
    message("  文件: ", basename(f))
    # 尝试读取
    tryCatch({
      df <- read.table(f, header=TRUE, stringsAsFactors=FALSE)
      message("    行数: ", nrow(df))
      if("p_SMR" %in% colnames(df)) {
        message("    SMR p范围: ", round(min(df$p_SMR, na.rm=TRUE), 4), " - ", round(max(df$p_SMR, na.rm=TRUE), 4))
      }
    }, error=function(e) message("    读取失败"))
  }
} else {
  message("\n[SMR] ❌ 无.smr产出")
}

message("\n========================================")
message("提取完成")
message("========================================")
