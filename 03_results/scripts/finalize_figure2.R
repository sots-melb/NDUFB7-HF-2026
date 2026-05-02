library(data.table)
library(ggplot2)
library(cowplot)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"
cat("Figure 2定稿...\n")

# A. GSE57338
g1 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
p1 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")
df1 <- data.table(Expression=as.numeric(unlist(g1[,-1])), Group=factor(p1$disease.status, levels=c("non-failing","idiopathic dilated CMP","ischemic")))
pA <- ggplot(df1, aes(x=Group, y=Expression, fill=Group)) + geom_boxplot(alpha=0.85, width=0.6) + geom_jitter(width=0.15, alpha=0.2, size=0.8) +
  scale_fill_manual(values=c("#999999","#377EB8","#E41A1C")) + labs(title="A. Affymetrix 1.1 ST (GSE57338, n=313)", y="log2", x="") +
  annotate("text", x=2.5, y=max(df1$Expression)*0.98, label="DCM vs ICM: p=0.050\nKW p=0.095", size=3, hjust=1) +
  theme_bw() + theme(legend.position="none", axis.text.x=element_text(face="bold"), plot.title=element_text(face="bold", size=10))

# B. GDS4772
g2 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_stratified_fixed.txt")
pB <- ggplot(g2, aes(x=Group, y=Expression, fill=Group)) + geom_boxplot(alpha=0.85, width=0.5) + geom_jitter(width=0.1, alpha=0.5, size=2) +
  scale_fill_manual(values=c("NF"="#999999","DCM"="#377EB8")) + labs(title="B. Affymetrix 1.0 ST (GDS4772, n=17)", y="log2", x="") +
  annotate("text", x=1.5, y=max(g2$Expression)*0.98, label="DCM vs NF: p=0.44", size=3) +
  theme_bw() + theme(legend.position="none", plot.title=element_text(face="bold", size=10))

# C. GSE116250
g3 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_stratified.txt")
pC <- ggplot(g3, aes(x=Disease, y=Expression, fill=Disease)) + geom_boxplot(alpha=0.85, width=0.6) + geom_jitter(width=0.15, alpha=0.3, size=1) +
  scale_fill_manual(values=c("non-failing"="#999999","dilated cardiomyopathy"="#377EB8","ischemic cardiomyopathy"="#E41A1C")) +
  labs(title="C. RNA-seq RPKM (GSE116250, n=64)", y="RPKM", x="") +
  annotate("text", x=2, y=max(g3$Expression)*0.95, label="DCM/ICM > NF\nKW p=0.021\n*137aa length bias", size=2.8, color="darkred") +
  theme_bw() + theme(legend.position="none", axis.text.x=element_text(size=8), plot.title=element_text(face="bold", size=10))

# D. GSE55296
g4 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_stratified.txt")
pD <- ggplot(g4, aes(x=Group, y=Expression, fill=Group)) + geom_boxplot(alpha=0.85, width=0.6) + geom_jitter(width=0.15, alpha=0.3, size=1) +
  scale_fill_manual(values=c("Control"="#999999","Dilated"="#377EB8","Ischemic"="#E41A1C")) +
  labs(title="D. RNA-seq Count (GSE55296, n=36)", y="Raw Counts", x="") +
  annotate("text", x=2, y=max(g4$Expression)*0.95, label="ICM>DCM>Control\n*small cohort bias", size=2.8, color="darkred") +
  theme_bw() + theme(legend.position="none", axis.text.x=element_text(size=8), plot.title=element_text(face="bold", size=10))

combined <- plot_grid(pA, pB, pC, pD, ncol=2, labels=c("A","B","C","D"))
title <- ggdraw() + draw_label("Figure 2. Multi-Platform NDUFB7 Expression in Heart Failure", fontface="bold", size=14)
final <- plot_grid(title, combined, ncol=1, rel_heights=c(0.08,1))
ggsave(file.path(out_dir, "Figure2_final_v5.png"), final, width=11, height=9, dpi=300)
cat("✅ Figure 2定稿已保存\n")
