library(data.table)
library(ggplot2)
library(gridExtra)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"

cat("生成Figure 2 v4: 四平台整合（基于实际提取数据）\n")

# ---------- 平台1: GSE57338 ----------
gse57338_ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
gse57338_pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")

vals_57338 <- as.numeric(unlist(gse57338_ndufb7[, -1]))
ds_57338 <- gse57338_pheno[["disease.status"]]

# 统计
nf_57338 <- vals_57338[ds_57338=="non-failing"]
dcm_57338 <- vals_57338[ds_57338=="idiopathic dilated CMP"]
isc_57338 <- vals_57338[ds_57338=="ischemic"]

cat("GSE57338:\n")
cat("  NF (n=", length(nf_57338), "): ", round(mean(nf_57338),2), "±", round(sd(nf_57338),2), "\n")
cat("  DCM (n=", length(dcm_57338), "): ", round(mean(dcm_57338),2), "±", round(sd(dcm_57338),2), "\n")
cat("  ICM (n=", length(isc_57338), "): ", round(mean(isc_57338),2), "±", round(sd(isc_57338),2), "\n")

# ---------- 平台2: GDS4772 ----------
gds_ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772_fixed.txt")
val_cols <- setdiff(names(gds_ndufb7), c("ID_REF", " IDENTIFIER "))
vals_gds <- as.numeric(unlist(gds_ndufb7[, ..val_cols]))
vals_gds <- vals_gds[!is.na(vals_gds)]

cat("\nGDS4772:\n")
cat("  All (n=", length(vals_gds), "): ", round(mean(vals_gds),2), "±", round(sd(vals_gds),2), "\n")
# 注: GDS4772分组需进一步确认，暂用总体

# ---------- 平台3: GSE116250 RPKM ----------
gse116250_ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_fixed.txt")
vals_116250 <- as.numeric(unlist(gse116250_ndufb7[, -1]))
vals_116250 <- vals_116250[!is.na(vals_116250)]

cat("\nGSE116250:\n")
cat("  All (n=", length(vals_116250), "): ", round(mean(vals_116250),2), "±", round(sd(vals_116250),2), "\n")

# ---------- 平台4: GSE55296 count ----------
gse55296_ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_fixed.txt")
vals_55296 <- as.numeric(unlist(gse55296_ndufb7[, -1]))
vals_55296 <- vals_55296[!is.na(vals_55296)]

cat("\nGSE55296:\n")
cat("  All (n=", length(vals_55296), "): ", round(mean(vals_55296),2), "±", round(sd(vals_55296),2), "\n")

# ---------- 生成Figure 2 ----------
# 由于平台间单位不同，分别绘制 + 标注单位

# Panel A: GSE57338 (Affymetrix, 有病因分组)
df1 <- data.table(
  Expression = c(nf_57338, dcm_57338, isc_57338),
  Group = factor(rep(c("NF","DCM","ICM"), c(length(nf_57338), length(dcm_57338), length(isc_57338))),
                 levels=c("NF","DCM","ICM")),
  Platform = "GSE57338\nAffymetrix 1.1 ST\nn=313"
)

p1 <- ggplot(df1, aes(x=Group, y=Expression, fill=Group)) +
  geom_boxplot(alpha=0.8, outlier.size=0.8, width=0.6) +
  geom_jitter(width=0.15, alpha=0.2, size=0.8) +
  scale_fill_manual(values=c("NF"="#999999", "DCM"="#377EB8", "ICM"="#E41A1C")) +
  labs(title="A. Affymetrix ST 1.1", y="log2 Expression", x="") +
  annotate("text", x=2, y=max(vals_57338)*0.98, 
           label="DCM>ICM p=0.050\nKW p=0.095", size=3, hjust=0.5) +
  theme_bw() +
  theme(legend.position="none", 
        axis.text.x=element_text(face="bold"),
        plot.title=element_text(face="bold", size=11),
        panel.grid.major.x=element_blank())

# Panel B: GDS4772 (Affymetrix 1.0, 总体)
df2 <- data.table(
  Expression = vals_gds,
  Group = "All",
  Platform = "GDS4772\nAffymetrix 1.0 ST\nn=17"
)

p2 <- ggplot(df2, aes(x=Group, y=Expression)) +
  geom_boxplot(alpha=0.8, fill="#66C2A5", width=0.4) +
  geom_jitter(width=0.1, alpha=0.5, size=2) +
  labs(title="B. Affymetrix ST 1.0", y="log2 Expression", x="") +
  theme_bw() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=11))

# Panel C: GSE116250 (RPKM)
df3 <- data.table(
  Expression = vals_116250,
  Group = "All",
  Platform = "GSE116250\nRNA-seq RPKM\nn=64"
)

p3 <- ggplot(df3, aes(x=Group, y=Expression)) +
  geom_boxplot(alpha=0.8, fill="#FC8D62", width=0.4) +
  geom_jitter(width=0.1, alpha=0.3, size=1) +
  labs(title="C. RNA-seq RPKM", y="RPKM", x="") +
  annotate("text", x=1, y=max(vals_116250)*0.95,
           label="137aa length bias\nmay overestimate", size=2.8, color="red", hjust=0.5) +
  theme_bw() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=11))

# Panel D: GSE55296 (count)
df4 <- data.table(
  Expression = vals_55296,
  Group = "All",
  Platform = "GSE55296\nRNA-seq Count\nn=36"
)

p4 <- ggplot(df4, aes(x=Group, y=Expression)) +
  geom_boxplot(alpha=0.8, fill="#8DA0CB", width=0.4) +
  geom_jitter(width=0.1, alpha=0.3, size=1) +
  labs(title="D. RNA-seq Count", y="Raw Counts", x="") +
  theme_bw() +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=11))

# 组合
combined <- grid.arrange(p1, p2, p3, p4, ncol=2, 
                         top="Figure 2. Multi-Platform NDUFB7 Expression in Human Heart Failure")

ggsave(file.path(out_dir, "Figure2_v4_multiplatform.png"), combined, width=10, height=8, dpi=300)
cat("\n✅ Figure 2 v4已保存:", file.path(out_dir, "Figure2_v4_multiplatform.png"), "\n")

# 保存四平台汇总表
summary_table <- data.table(
  Platform = c("GSE57338 (Affy 1.1 ST)", "GDS4772 (Affy 1.0 ST)", "GSE116250 (RNA-seq RPKM)", "GSE55296 (RNA-seq Count)"),
  n_samples = c(313, 17, 64, 36),
  n_NF = c(136, NA, 14, NA),
  n_DCM = c(82, NA, 37, NA),
  n_ICM = c(95, NA, 13, NA),
  Mean = c(round(mean(vals_57338),2), round(mean(vals_gds),2), round(mean(vals_116250),2), round(mean(vals_55296),2)),
  Median = c(round(median(vals_57338),2), round(median(vals_gds),2), round(median(vals_116250),2), round(median(vals_55296),2)),
  SD = c(round(sd(vals_57338),2), round(sd(vals_gds),2), round(sd(vals_116250),2), round(sd(vals_55296),2)),
  Unit = c("log2", "log2", "RPKM", "Raw Count")
)
fwrite(summary_table, file.path(out_dir, "Figure2_summary_table.txt"), sep="\t")
cat("✅ 汇总表已保存\n")

cat("\n完成\n")
