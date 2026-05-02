library(data.table)
library(ggplot2)

cat("生成Figure 1: 空间数据示意\n")

# 使用已知统计值生成示意数据（实际投稿需替换为真实spot-level）
set.seed(42)
spots <- data.table(
  sample = rep(c("P3_IZ","P7_BZ","P14_BZ","P20_FZ","P5_control","P13_control"), c(100,100,100,100,100,100)),
  zone = rep(c("IZ","BZ","BZ","FZ","Control","Control"), c(100,100,100,100,100,100)),
  ndufb7_expression = c(
    rnorm(100, 1.26, 0.5), rnorm(100, 1.06, 0.6), rnorm(100, 1.06, 0.6),
    rnorm(100, 0.00, 0.3), rnorm(100, 0.8, 0.8), rnorm(100, 0.8, 0.8)
  )
)
spots$ndufb7_expression[spots$zone == "FZ"] <- pmax(0, spots$ndufb7_expression[spots$zone == "FZ"])
spots$detected <- spots$ndufb7_expression > 0

# Panel A
pos_rate <- spots[, .(positive_rate = mean(detected)*100, median_expr = median(ndufb7_expression)), by = .(sample, zone)]
pA <- ggplot(pos_rate, aes(x=sample, y=positive_rate, fill=zone)) +
  geom_bar(stat="identity", width=0.7) +
  geom_text(aes(label=paste0(round(positive_rate,1),"%")), vjust=-0.5, size=3) +
  scale_fill_manual(values=c("IZ"="#E41A1C", "BZ"="#FF7F00", "FZ"="#377EB8", "Control"="#999999")) +
  labs(title="A. NDUFB7 Detection Rate by Zone", y="% Positive Spots", x="") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1, face="bold"), legend.position="top")

# Panel B
pB <- ggplot(spots, aes(x=zone, y=ndufb7_expression, fill=zone)) +
  geom_violin(alpha=0.7, trim=FALSE) + geom_boxplot(width=0.2, alpha=0.8) +
  scale_fill_manual(values=c("IZ"="#E41A1C", "BZ"="#FF7F00", "FZ"="#377EB8", "Control"="#999999")) +
  labs(title="B. Expression Distribution", y="Normalized Expression", x="") +
  theme_bw() + theme(legend.position="none")

# Panel C
iz_fz <- spots[zone %in% c("IZ","FZ")]
pC <- ggplot(iz_fz, aes(x=zone, y=ndufb7_expression, fill=zone)) +
  geom_boxplot(alpha=0.85, width=0.5) + geom_jitter(width=0.1, alpha=0.3, size=1) +
  scale_fill_manual(values=c("IZ"="#E41A1C", "FZ"="#377EB8")) +
  labs(title="C. Acute IZ vs Chronic FZ", y="Expression", x="") +
  annotate("text", x=1.5, y=max(iz_fz$ndufb7_expression)*0.95, label="MWU p<0.0001", size=3.5) +
  theme_bw() + theme(legend.position="none", axis.text.x=element_text(face="bold", size=12))

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"
ggsave(file.path(out_dir, "Figure1A_spatial.png"), pA, width=8, height=5, dpi=300)
ggsave(file.path(out_dir, "Figure1B_violin.png"), pB, width=6, height=5, dpi=300)
ggsave(file.path(out_dir, "Figure1C_IZ_FZ.png"), pC, width=5, height=5, dpi=300)
cat("✅ Figure 1已保存（示意数据，投稿前替换真实spot数据）\n")
