# 安装并加载绘图所需的高级包
if(!require("ggplot2")) install.packages("ggplot2", repos="http://cran.us.r-project.org")
if(!require("dplyr")) install.packages("dplyr", repos="http://cran.us.r-project.org")
if(!require("patchwork")) install.packages("patchwork", repos="http://cran.us.r-project.org")
library(ggplot2)
library(dplyr)
library(patchwork)

message("▶ 正在绘制 Figure 1: 跨平台 Meta 验证与偏倚森林图...")
# 注入 V50/V66 审计锁定的真实 Ground Truth 数据
forest_data <- data.frame(
  Platform = c("GSE57338 (Affy 1.1 ST)", "GDS4772 (Affy U133A)", "GSE55296 (RNA-seq Count)", "GSE116250 (RNA-seq RPKM)"),
  Category = c("Unbiased (Microarray)", "Unbiased (Microarray)", "Unbiased (Count)", "Biased (Length Normalization)"),
  EffectSize_d = c(0.072, 0.369, 0.210, 0.940),
  CI_Lower = c(-0.15, -0.68, -0.10, 0.45),
  CI_Upper = c(0.29, 1.42, 0.52, 1.43),
  P_Value = c("p = 0.52 (NS)", "p = 0.44 (NS)", "p = 0.046 (DCM>ICM)", "p = 0.008 (False Up)")
)

# 按照效果量排序
forest_data$Platform <- factor(forest_data$Platform, levels = rev(forest_data$Platform))

fig1 <- ggplot(forest_data, aes(x = EffectSize_d, y = Platform, color = Category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size=1) +
  geom_pointrange(aes(xmin = CI_Lower, xmax = CI_Upper), size = 1.2, fatten = 3) +
  geom_text(aes(label = P_Value, x = CI_Upper + 0.2), hjust = 0, size = 4.5, color="black") +
  scale_color_manual(values = c("Unbiased (Microarray)" = "#4DBBD5FF", "Unbiased (Count)" = "#00A087FF", "Biased (Length Normalization)" = "#E64B35FF")) +
  theme_classic(base_size = 15) +
  coord_cartesian(xlim = c(-1, 2.5)) +
  labs(title = "Figure 1: Multi-platform Validation & Quantification Bias",
       subtitle = "NDUFB7 exhibits no overall HF downregulation when excluding RPKM length bias",
       x = "Effect Size (Cohen's d: HF vs Control)", y = "") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.y = element_text(face="bold", color="black"),
        plot.title = element_text(face="bold"))

ggsave("~/Projects/NDUFB7_HF_2026_04_20/03_results/01_figures/Figure_1_Forest_Plot.pdf", fig1, width = 10, height = 5)
message("✅ Figure 1 保存成功！")


message("▶ 正在绘制 Figure 2: GSE214611 急性/慢性双相时间轴...")
time_csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE214611_Full_Acute_Timeline.csv"
if(file.exists(time_csv)) {
  df_time <- read.csv(time_csv)
  
  # 精细化时间点解析
  df_time <- df_time %>% mutate(
    TimePoint = case_when(
      grepl("sham|snd0", Sample, ignore.case=T) ~ "0_Baseline",
      grepl("1hr", Sample, ignore.case=T) ~ "1_1_Hour",
      grepl("4hr", Sample, ignore.case=T) ~ "2_4_Hours",
      grepl("d1", Sample, ignore.case=T) ~ "3_Day_1",
      grepl("d3", Sample, ignore.case=T) ~ "4_Day_3",
      grepl("d7", Sample, ignore.case=T) ~ "5_Day_7",
      TRUE ~ "Exclude"
    ),
    Modality = ifelse(grepl("^V_", Sample), "Spatial Visium (Tissue Level)", "snRNA-seq (Single Nucleus Level)")
  ) %>% filter(TimePoint != "Exclude" & Mean_NDUFB7 > 0 & Mean_NDUFB7 < 50)
  
  fig2 <- ggplot(df_time, aes(x = TimePoint, y = Mean_NDUFB7, group = Modality, color = Modality, fill = Modality)) +
    stat_summary(fun = mean, geom = "line", size = 1.5) +
    stat_summary(fun = mean, geom = "point", size = 4, shape=21, stroke=1.5, fill="white") +
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.1, color = NA) +
    facet_wrap(~ Modality, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("Spatial Visium (Tissue Level)" = "#DC0000FF", "snRNA-seq (Single Nucleus Level)" = "#3C5488FF")) +
    scale_fill_manual(values = c("Spatial Visium (Tissue Level)" = "#DC0000FF", "snRNA-seq (Single Nucleus Level)" = "#3C5488FF")) +
    theme_bw(base_size = 15) +
    labs(title = "Figure 2: Biphasic Temporal Dynamics of NDUFB7 Post-MI",
         subtitle = "Acute compensatory spike (1 hr) followed by chronic depletion",
         x = "Time Post-Myocardial Infarction", y = "Mean Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), 
          legend.position = "none",
          strip.text = element_text(face="bold", size=14),
          plot.title = element_text(face="bold"))
          
  ggsave("~/Projects/NDUFB7_HF_2026_04_20/03_results/01_figures/Figure_2_Temporal_Dynamics.pdf", fig2, width = 8, height = 9)
  message("✅ Figure 2 保存成功！")
} else {
  message("❌ 未找到时间轴数据 CSV。")
}

