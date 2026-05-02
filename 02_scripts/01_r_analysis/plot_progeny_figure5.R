library(data.table)
library(ggplot2)

cat("生成Figure 5: PROGENy虚拟敲除通路差异\n")

# 读取结果
d <- fread("~/Projects/NDUFB7_HF_2026_04_20/03_results/virtual_ko_pathway_v4.txt")

# 计算效应方向（低表达组相对于高表达组的变化百分比）
d$Effect <- (d$Low_Mean - d$High_Mean) / d$High_Mean * 100
d$Direction <- ifelse(d$Effect < 0, "Down in NDUFB7-low", "Up in NDUFB7-low")
d$Significant <- d$Pval < 0.05

# 排序：按效应大小
d$Pathway <- factor(d$Pathway, levels=d$Pathway[order(d$Effect)])

# 绘制
p <- ggplot(d, aes(x=Pathway, y=Effect, fill=Direction, alpha=Significant)) +
  geom_bar(stat="identity", width=0.7, color="black", size=0.3) +
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) +
  geom_text(aes(label=ifelse(Significant, paste0("p=", format(Pval, digits=2)), "ns")), 
            hjust=ifelse(d$Effect > 0, -0.1, 1.1), size=3, color="black") +
  scale_fill_manual(values=c("Down in NDUFB7-low"="#377EB8", "Up in NDUFB7-low"="#E41A1C")) +
  scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.4)) +
  coord_flip() +
  labs(
    title="Figure 5. In Silico Knockdown of NDUFB7 Predicts Pathway Alterations",
    subtitle="Bottom 25% vs Top 25% NDUFB7 expressors in GSE57338 (n=313)",
    x="", 
    y="Relative Change in Pathway Activity (%)",
    caption="* p<0.05; ** p<1×10⁻⁶; ns = not significant"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face="bold", size=12),
    plot.subtitle = element_text(size=9, color="gray30"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=9),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"
ggsave(file.path(out_dir, "Figure5_virtual_ko_pathways.png"), p, width=9, height=6, dpi=300)
cat("✅ Figure 5已保存:", file.path(out_dir, "Figure5_virtual_ko_pathways.png"), "\n")

# 同时生成Xenium细胞类型参考图（Figure 6或Supplementary）
xenium <- fread("~/Projects/NDUFB7_HF_2026_04_20/03_results/xenium_cell_type_reference.txt")
xenium <- xenium[order(-Percentage)]
xenium$CellType <- factor(xenium$CellType, levels=xenium$CellType)

p2 <- ggplot(xenium[1:10], aes(x=CellType, y=Percentage, fill=CellType)) +
  geom_bar(stat="identity", width=0.7, alpha=0.85) +
  geom_text(aes(label=paste0(Percentage, "%")), vjust=-0.3, size=3) +
  scale_fill_manual(values=c(
    "Cardiomyocytes"="#E41A1C", "CD8+ T cells"="#377EB8", "Endothelial"="#4DAF4A",
    "Macrophages"="#984EA3", "Fibroblasts"="#FF7F00", "Pericytes"="#FFFF33",
    "NK"="#A65628", "Plasma"="#F781BF", "Activated endothelial"="#999999",
    "cDC2"="#66C2A5"
  )) +
  labs(
    title="Supplementary Figure. Cardiac Cell Type Composition (Xenium, n=162,638)",
    x="", y="Percentage (%)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=9),
    legend.position = "none",
    plot.title = element_text(face="bold", size=10)
  )

ggsave(file.path(out_dir, "SuppFigure_Xenium_celltypes.png"), p2, width=8, height=5, dpi=300)
cat("✅ Supplementary Figure已保存\n")

cat("\n完成\n")
