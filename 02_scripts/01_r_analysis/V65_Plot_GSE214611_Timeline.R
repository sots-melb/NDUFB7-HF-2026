if(!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if(!require("dplyr", quietly = TRUE)) install.packages("dplyr")
library(ggplot2)
library(dplyr)

message("▶ 开始处理时间轴数据...")
in_csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE214611_Full_Acute_Timeline.csv"
df <- read.csv(in_csv)

# 过滤异常庞大/可疑的空斑点 (如sn4hr 6百万细胞明显是未过滤的空滴)
df <- df %>% filter(Cells > 500 & Cells < 100000)

# 解析时间点与数据模态
df <- df %>% mutate(
  Modality = ifelse(grepl("^V_", Sample), "Spatial Visium", "snRNA-seq"),
  TimePoint = case_when(
    grepl("sham|snd0", Sample, ignore.case=T) ~ "0_Baseline",
    grepl("1hr", Sample, ignore.case=T) ~ "1_1_Hour",
    grepl("4hr", Sample, ignore.case=T) ~ "2_4_Hours",
    grepl("d1", Sample, ignore.case=T) ~ "3_Day_1",
    grepl("d3", Sample, ignore.case=T) ~ "4_Day_3",
    grepl("d7", Sample, ignore.case=T) ~ "5_Day_7",
    TRUE ~ "Exclude"
  )
) %>% filter(TimePoint != "Exclude")

# 绘图：分为两个 Panel 分别展示 Visium 和 snRNA
p <- ggplot(df, aes(x = TimePoint, y = Mean_NDUFB7, group = Modality, color = Modality)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ Modality, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("Spatial Visium" = "#D55E00", "snRNA-seq" = "#0072B2")) +
  labs(title = "Biphasic Evolution of NDUFB7 in Acute Myocardial Infarction",
       subtitle = "Acute compensation (1 hr) followed by chronic depletion",
       x = "Post-MI Timeline", y = "Mean Expression (Normalized)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

out_pdf <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/01_figures/GSE214611_Biphasic_Timeline.pdf"
ggsave(out_pdf, p, width = 7, height = 8)
message(paste("✅ 完美！双相时间折线图已保存至:", out_pdf))
