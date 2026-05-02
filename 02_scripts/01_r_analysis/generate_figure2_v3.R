library(data.table)
library(ggplot2)
library(gridExtra)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"

cat("生成Figure 2: 四平台跨数据集整合\n")

# 平台1: GSE57338 (Affymetrix 1.1 ST, n=313)
# 从之前保存的文件读取
gse57338_ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
gse57338_pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")

# 平台2: GDS4772 (Affymetrix 1.0 ST, n=29) - 待STEP 1提取
# 平台3: GSE116250 (RNA-seq RPKM, n=64) - 待STEP 1提取
# 平台4: GSE55296 (RNA-seq count, n=36) - 待STEP 1提取

cat("注意: 需要STEP 1完成后，手动整合四平台数据到本脚本\n")
cat("当前GSE57338数据可用，其他平台待提取\n")

# 先绘制GSE57338病因分层图
if(nrow(gse57338_ndufb7) > 0 && nrow(gse57338_pheno) > 0) {
  ndufb7_vals <- as.numeric(unlist(gse57338_ndufb7[, -1]))
  ds <- gse57338_pheno[["disease.status"]]
  
  df_57338 <- data.table(
    Expression = ndufb7_vals,
    Group = ds,
    Platform = "GSE57338\nAffymetrix 1.1 ST\nn=313"
  )
  
  p1 <- ggplot(df_57338, aes(x=Group, y=Expression, fill=Group)) +
    geom_boxplot(alpha=0.7, outlier.size=1) +
    geom_jitter(width=0.2, alpha=0.3, size=1) +
    scale_fill_manual(values=c("non-failing"="gray70", "idiopathic dilated CMP"="#377EB8", "ischemic"="#E41A1C")) +
    labs(title="GSE57338", y="log2 Expression", x="") +
    theme_bw() +
    theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1))
  
  ggsave(file.path(out_dir, "Figure2_panel_GSE57338.png"), p1, width=4, height=4, dpi=300)
  cat("✅ GSE57338 panel已保存\n")
}

cat("完成（部分）\n")
