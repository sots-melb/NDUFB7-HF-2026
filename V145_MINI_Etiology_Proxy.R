#!/usr/bin/env Rscript
library(ggplot2)
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
df <- read.csv("03_results/V141_GSE57338_Bulk_Fix/V141_bulk_ferroptosis.csv", stringsAsFactors = FALSE)

# 用NDUFB7三分位作为"病因严重程度"代理（低=ICM-like, 高=NF-like）
df$Severity <- cut(df$NDUFB7, breaks = quantile(df$NDUFB7, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                   labels = c("Severe (Low)", "Moderate", "Mild (High)"), include.lowest = TRUE)

kw <- kruskal.test(ferroptosis_score ~ Severity, data = df)
p <- ggplot(df, aes(x = Severity, y = ferroptosis_score, fill = Severity)) +
  geom_boxplot(alpha = 0.8) + geom_jitter(width = 0.1, alpha = 0.3) +
  scale_fill_manual(values = c("Severe (Low)"="#D55E00", "Moderate"="#E69F00", "Mild (High)"="#0072B2")) +
  labs(title = "Ferroptosis by NDUFB7 Severity Tertiles (GSE57338)",
       subtitle = paste0("Kruskal-Wallis p = ", signif(kw$p.value, 2))) +
  theme_minimal() + theme(legend.position = "none")
ggsave("03_results/V145_GSE57338_Etiology/V145_Severity_Tertile_Proxy.png", p, width = 6, height = 5, dpi = 300, bg = "white")
cat("[MINI] Severity proxy saved. True etiology labels require manual GEO inspection.\n")
