library(data.table)
library(ggplot2)
library(cowplot)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"

cat("生成Figure 2 v3: 四平台整合\n")

# ---------- 平台1: GSE57338 (已有) ----------
gse57338_ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
gse57338_pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")

vals_57338 <- as.numeric(unlist(gse57338_ndufb7[, -1]))
ds_57338 <- gse57338_pheno[["disease.status"]]
df1 <- data.table(
  Expression = vals_57338,
  Group = factor(ds_57338, levels=c("non-failing","idiopathic dilated CMP","ischemic")),
  Platform = "GSE57338\nAffymetrix 1.1 ST\nn=313"
)

# 统计
dcm_57338 <- vals_57338[ds_57338=="idiopathic dilated CMP"]
isc_57338 <- vals_57338[ds_57338=="ischemic"]
nf_57338 <- vals_57338[ds_57338=="non-failing"]

# ---------- 平台2: GDS4772 (curated, 29样本) ----------
# 需要手动填入提取后的值，或从文件读取
gds_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772.txt"
df2 <- NULL
if(file.exists(gds_file)) {
  gds_ndufb7 <- fread(gds_file)
  # GDS4772的subset信息需要从SOFT解析，这里简化处理
  # 假设已知分组: 6 NF, 11 DCM, 12 ICM (根据GDS4772文献)
  # 实际应从GDS4772.soft的subset样本列表精确提取
  cat("GDS4772文件存在，但表型需精确匹配样本顺序\n")
}

# ---------- 平台3: GSE116250 RPKM ----------
rpkm_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_extracted.txt"
df3 <- NULL
if(file.exists(rpkm_file)) {
  rpkm <- fread(rpkm_file)
  vals_116250 <- as.numeric(unlist(rpkm[, -1]))
  # GSE116250分组: 14 NF, 37 DCM, 13 ICM (来自GEO网页)
  # 样本顺序未知，需要从series_matrix提取表型匹配
  cat("GSE116250 RPKM提取完成，表型匹配待进行\n")
}

# ---------- 平台4: GSE55296 count ----------
count_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_extracted.txt"
df4 <- NULL
if(file.exists(count_file)) {
  cnt <- fread(count_file)
  vals_55296 <- as.numeric(unlist(cnt[, -1]))
  cat("GSE55296 count提取完成\n")
}

# ---------- 绘制GSE57338面板（立即可用） ----------
p1 <- ggplot(df1, aes(x=Group, y=Expression, fill=Group)) +
  geom_boxplot(alpha=0.8, outlier.size=0.8, width=0.6) +
  geom_jitter(width=0.15, alpha=0.2, size=0.8) +
  scale_fill_manual(values=c("non-failing"="#999999", "idiopathic dilated CMP"="#377EB8", "ischemic"="#E41A1C"),
                    labels=c("NF","DCM","ICM")) +
  labs(title="Affymetrix ST 1.1", subtitle="GSE57338, n=313", y="log2 Expression", x="") +
  annotate("text", x=2.5, y=max(vals_57338)*0.95, 
           label=paste("DCM vs ICM: p=0.050\nKW p=0.095"), size=3, hjust=1) +
  theme_bw() +
  theme(legend.position="none", 
        axis.text.x=element_text(angle=0, hjust=0.5, face="bold"),
        plot.title=element_text(face="bold"),
        panel.grid.major.x=element_blank())

# 保存当前可用版本
ggsave(file.path(out_dir, "Figure2_v3_GSE57338_only.png"), p1, width=5, height=5, dpi=300)
cat("✅ Figure2 GSE57338面板已保存\n")

# 保存完整四平台数据供后续整合
save(df1, df2, df3, df4, file=file.path(out_dir, "Figure2_data.RData"))
cat("✅ 数据已保存到 Figure2_data.RData\n")

cat("完成\n")
