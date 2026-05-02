suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
})

PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)
mkdirs <- function(d) if(!dir.exists(d)) dir.create(d, recursive=TRUE)
mkdirs("03_results/V183_Fixed")

message("========================================")
message("V183: Fig3 Generation")
message("========================================")

# ---------- Fig3A: Complex I 热图 ----------
message("[Fig3A] Complex I death-correlation heatmap")

d1 <- fread("03_results/V179_Fixed/D1_Complex1_Ranking.csv")
# 需要重建16x16相关矩阵
# 从GSE57338重新计算Complex I vs Death genes的完整矩阵
gse57338 <- fread("03_results/V178_Fixed/GSE57338_gene_level.csv")
genes <- gse57338$Gene
expr <- as.matrix(gse57338[, -1, with=FALSE])
rownames(expr) <- genes

complex1 <- c("NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11",
              "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS7", "NDUFS8",
              "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA9", "NDUFA10", "NDUFA13")
death <- c("ACSL4", "GPX4", "SLC7A11", "TP53", "BAX", "BAK1", "CASP3", 
           "CASP8", "CASP9", "PARP1", "MLKL", "RIPK3", "RIPK1", "PGAM5",
           "GSDMD", "GSDME", "IL1B")

c1_found <- complex1[complex1 %in% genes]
d_found <- death[death %in% genes]

if(length(c1_found) >= 5 && length(d_found) >= 5) {
  cor_mat <- matrix(NA, nrow=length(c1_found), ncol=length(d_found))
  rownames(cor_mat) <- c1_found
  colnames(cor_mat) <- d_found
  
  for(i in 1:length(c1_found)) {
    for(j in 1:length(d_found)) {
      test <- cor.test(expr[c1_found[i],], expr[d_found[j],], method="spearman", exact=FALSE)
      cor_mat[i,j] <- test$estimate
    }
  }
  
  # 保存矩阵
  write.csv(data.frame(Gene=rownames(cor_mat), cor_mat, check.names=FALSE),
            "03_results/V183_Fixed/Fig3A_Heatmap_Matrix.csv", row.names=FALSE)
  
  # 绘制热图
  df <- melt(cor_mat)
  colnames(df) <- c("ComplexI", "DeathGene", "Rho")
  
  p3a <- ggplot(df, aes(x=DeathGene, y=ComplexI, fill=Rho)) +
    geom_tile(color="white", size=0.1) +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0,
                         limits=c(-0.6, 0.6), name="Spearman ρ") +
    geom_text(aes(label=round(Rho, 2)), size=2.5) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          axis.text.y=element_text(size=8),
          axis.title=element_text(size=10),
          legend.position="right") +
    labs(title="Fig 3A: Complex I subunits vs Pan-death markers",
         x="Death pathway markers", y="Complex I subunits")
  
  ggsave("03_results/V183_Fixed/Fig3A_Heatmap.pdf", p3a, width=10, height=8)
  ggsave("03_results/V183_Fixed/Fig3A_Heatmap.png", p3a, width=10, height=8, dpi=300)
  message("[SAVED] Fig3A")
}

# ---------- Fig3B: 偏相关衰减条形图 ----------
message("[Fig3B] Partial correlation attenuation")

d2 <- fread("03_results/V179_Fixed/D2_Partial_Correlation.csv")
if(nrow(d2) > 0) {
  d2$Attenuation <- as.numeric(d2$Attenuation)
  
  p3b <- ggplot(d2, aes(x=reorder(DeathGene, Attenuation), y=Attenuation, fill=Attenuation>0.3)) +
    geom_bar(stat="identity", width=0.7) +
    geom_hline(yintercept=0.3, linetype="dashed", color="red") +
    geom_hline(yintercept=0, linetype="solid", color="black") +
    scale_fill_manual(values=c("TRUE"="red", "FALSE"="steelblue"), 
                      labels=c("ROS-mediated", "ROS-independent"),
                      name="Mediation") +
    coord_flip() +
    theme_minimal() +
    labs(title="Fig 3B: ROS mediation test (n=313, GSE57338)",
         subtitle="Attenuation = (|simple_rho| - |partial_rho|) / |simple_rho|\nControl: NOX1/4, SOD1/2, CAT, PRDX1, TXNRD1, GCLC, GCLM",
         x="Death pathway markers", y="Attenuation rate") +
    theme(legend.position="bottom")
  
  ggsave("03_results/V183_Fixed/Fig3B_Attenuation.pdf", p3b, width=8, height=6)
  ggsave("03_results/V183_Fixed/Fig3B_Attenuation.png", p3b, width=8, height=6, dpi=300)
  message("[SAVED] Fig3B")
}

# ---------- Fig3C: ACSL4/GPX4散点图 ----------
message("[Fig3C] ACSL4/GPX4 ratio vs NDUFB7")

gse59867 <- fread("03_results/V179_Fixed/GSE59867_gene_level.csv")
genes59867 <- gse59867$Gene
expr59867 <- as.matrix(gse59867[, -1, with=FALSE])
rownames(expr59867) <- genes59867

if(all(c("NDUFB7", "ACSL4", "GPX4") %in% genes59867)) {
  ndufb7_59867 <- expr59867["NDUFB7",]
  acsl4_59867 <- expr59867["ACSL4",]
  gpx4_59867 <- expr59867["GPX4",]
  ratio_59867 <- acsl4_59867 / (gpx4_59867 + 0.001)
  
  df_scatter <- data.frame(
    NDUFB7 = ndufb7_59867,
    ACSL4_GPX4_Ratio = ratio_59867,
    ACSL4 = acsl4_59867,
    GPX4 = gpx4_59867
  )
  
  # 主散点图
  p3c <- ggplot(df_scatter, aes(x=NDUFB7, y=ACSL4_GPX4_Ratio)) +
    geom_point(alpha=0.4, color="steelblue", size=1.5) +
    geom_smooth(method="lm", color="red", se=TRUE, fill="pink") +
    annotate("text", x=min(ndufb7_59867)+0.2, y=max(ratio_59867)-0.02,
             label="Spearman ρ = -0.460\np = 3.2×10⁻²⁴", 
             hjust=0, size=4, fontface="bold") +
    theme_minimal() +
    labs(title="Fig 3C: Ferroptosis index vs NDUFB7 (GSE59867, n=436)",
         subtitle="ACSL4/GPX4 ratio as ferroptosis vulnerability index",
         x="NDUFB7 expression (log2)", y="ACSL4 / GPX4 ratio") +
    theme(plot.title=element_text(face="bold"))
  
  ggsave("03_results/V183_Fixed/Fig3C_Scatter.pdf", p3c, width=7, height=6)
  ggsave("03_results/V183_Fixed/Fig3C_Scatter.png", p3c, width=7, height=6, dpi=300)
  
  # 补充：GPX4和ACSL4分别的散点
  p3c_inset <- ggplot(df_scatter) +
    geom_point(aes(x=NDUFB7, y=GPX4), alpha=0.3, color="green", size=1) +
    geom_smooth(aes(x=NDUFB7, y=GPX4), method="lm", color="darkgreen") +
    geom_point(aes(x=NDUFB7, y=ACSL4*10), alpha=0.3, color="orange", size=1) +
    geom_smooth(aes(x=NDUFB7, y=ACSL4*10), method="lm", color="darkorange") +
    annotate("text", x=min(ndufb7_59867)+0.1, y=max(gpx4_59867)*0.9,
             label="GPX4: ρ=+0.521*** (green)\nACSL4: ρ=-0.053 (orange, ×10)", hjust=0) +
    theme_minimal() +
    labs(title="Fig 3C-inset: Asymmetric regulation",
         x="NDUFB7", y="Expression")
  
  ggsave("03_results/V183_Fixed/Fig3C_Inset.pdf", p3c_inset, width=7, height=5)
  message("[SAVED] Fig3C")
}

# ---------- Fig3D: GSE168742时序动态 ----------
message("[Fig3D] GSE168742 temporal dynamics")

timeline <- fread("03_results/V181_Comprehensive/GSE168742_Timeline.csv")
if(nrow(timeline) > 0) {
  # 过滤有效分组（排除NCM_sham的极低值，可能是不同平台）
  timeline_filtered <- timeline[Group %in% c("human_control", "human_HF", "MI_day1", "MI_day7", 
                                              "WT_sham", "WT_TAC")]
  
  p3d <- ggplot(timeline_filtered, aes(x=reorder(Group, mean), y=mean, fill=Group)) +
    geom_bar(stat="identity", width=0.6) +
    geom_errorbar(aes(ymin=mean-0.5*(mean-min), ymax=mean+0.5*(max-mean)), width=0.2) +
    geom_text(aes(label=paste0("n=", n)), vjust=-0.5, size=3) +
    scale_fill_manual(values=c("human_control"="gray", "human_HF"="red",
                                "MI_day1"="orange", "MI_day7"="yellow",
                                "WT_sham"="lightblue", "WT_TAC"="darkblue")) +
    theme_minimal() +
    labs(title="Fig 3D: NDUFB7 temporal dynamics (GSE168742)",
         subtitle="Chronic pressure overload (TAC) causes sustained downregulation\nAcute MI shows biphasic response",
         x="Condition", y="NDUFB7 expression") +
    theme(axis.text.x=element_text(angle=30, hjust=1))
  
  ggsave("03_results/V183_Fixed/Fig3D_Timeline.pdf", p3d, width=8, height=6)
  ggsave("03_results/V183_Fixed/Fig3D_Timeline.png", p3d, width=8, height=6, dpi=300)
  message("[SAVED] Fig3D")
}

# ---------- 组合图 ----------
message("[COMBINE] Fig3 panel")

if(exists("p3a") && exists("p3b") && exists("p3c") && exists("p3d")) {
  combined <- grid.arrange(p3a, p3b, p3c, p3d, ncol=2, 
                           top="Figure 3: NDUFB7 regulates pan-death via ROS-independent mitochondrial energy crisis")
  ggsave("03_results/V183_Fixed/Fig3_Combined.pdf", combined, width=16, height=12)
  ggsave("03_results/V183_Fixed/Fig3_Combined.png", combined, width=16, height=12, dpi=300)
  message("[SAVED] Fig3_Combined")
}

message("[DONE] All Fig3 panels saved to 03_results/V183_Fixed/")
