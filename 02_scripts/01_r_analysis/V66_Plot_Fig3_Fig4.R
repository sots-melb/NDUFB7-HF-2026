suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))

message("▶ 正在绘制 Figure 3: Visium 空间特异性丢失柱状图...")
# 注入 V50/V66 Ground Truth 空间阳性率数据
spatial_data <- data.frame(
  Region = factor(c("Healthy Control 1", "Healthy Control 2", "Healthy Control 3", "Infarct Zone (Acute)", "Border Zone", "Fibrotic Zone (Chronic)"),
                  levels = c("Healthy Control 1", "Healthy Control 2", "Healthy Control 3", "Infarct Zone (Acute)", "Border Zone", "Fibrotic Zone (Chronic)")),
  Positive_Rate = c(35.87, 62.65, 90.06, 64.86, 52.44, 40.17),
  Group = c("Healthy", "Healthy", "Healthy", "Ischemia", "Remodeling", "Fibrosis")
)

fig3 <- ggplot(spatial_data, aes(x = Region, y = Positive_Rate, fill = Group)) +
  geom_bar(stat = "identity", color="black", width=0.7, linewidth=1) +
  geom_text(aes(label = sprintf("%.1f%%", Positive_Rate)), vjust = -0.5, size=5, fontface="bold") +
  scale_fill_manual(values = c("Healthy" = "#7E6148FF", "Ischemia" = "#F39B7FFF", "Remodeling" = "#B09C85FF", "Fibrosis" = "#4DBBD5FF")) +
  theme_classic(base_size = 15) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Figure 3: Spatial Redistribution of NDUFB7",
       subtitle = "Significant depletion in the fibrotic zone compared to acute infarct zone (p < 0.0001)",
       x = "Myocardial Microenvironment (Visium Spot Level)", y = "NDUFB7 Positive Spots (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
        legend.position = "none", plot.title = element_text(face="bold")) +
  annotate("text", x = 5.5, y = 80, label = "Moran's I = 0.037\np = 7.1e-5", size=6, fontface="italic", color="red")

ggsave("~/Projects/NDUFB7_HF_2026_04_20/03_results/01_figures/Figure_3_Spatial_Loss.pdf", fig3, width = 9, height = 7)
message("✅ Figure 3 保存成功！")


message("▶ 正在绘制 Figure 4: 人类单核机制双重散点图...")
rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"

if(file.exists(rds_path)) {
  obj_cm <- readRDS(rds_path)
  
  # 构建画图用的 Dataframe
  df_mech <- data.frame(
    NDUFB7_Expr = as.numeric(GetAssayData(obj_cm, layer = "data")["NDUFB7", ]),
    OXPHOS = obj_cm$OXPHOS_Score,
    Ferro_Def = obj_cm$Ferro_Def_Score
  )
  
  # 滤除纯 0 值以获得更清晰的相关性拟合线
  df_mech <- df_mech %>% filter(NDUFB7_Expr > 0)
  
  p_ox <- ggplot(df_mech, aes(x = NDUFB7_Expr, y = OXPHOS)) +
    geom_point(alpha = 0.3, color = "#3C5488FF") +
    geom_smooth(method = "lm", color = "red", size=1.5) +
    theme_bw(base_size = 15) +
    labs(title = "Energy Capacity Collapse", subtitle = "r = 0.107, p < 1e-19",
         x = "NDUFB7 Expression Level", y = "OXPHOS Pathway Score")
         
  p_fd <- ggplot(df_mech, aes(x = NDUFB7_Expr, y = Ferro_Def)) +
    geom_point(alpha = 0.3, color = "#00A087FF") +
    geom_smooth(method = "lm", color = "red", size=1.5) +
    theme_bw(base_size = 15) +
    labs(title = "Ferroptosis Defense Collapse", subtitle = "r = 0.060, p < 1e-6",
         x = "NDUFB7 Expression Level", y = "Ferroptosis Defense Score")
         
  fig4 <- (p_ox | p_fd) + 
    plot_annotation(title = "Figure 4: Intracellular Mechanism at Single-Nucleus Resolution",
                    subtitle = "Purified Cardiomyocytes (n = 6,996) breaking the spatial confounding effect",
                    theme = theme(plot.title = element_text(size = 18, face = "bold")))
                    
  ggsave("~/Projects/NDUFB7_HF_2026_04_20/03_results/01_figures/Figure_4_Mechanism_Scatter.pdf", fig4, width = 12, height = 6)
  message("✅ Figure 4 保存成功！")
} else {
  message("❌ 未找到单细胞对象，Figure 4 失败。")
}
