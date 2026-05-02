library(data.table)
library(ggplot2)

cat("生成Figure 3...\n")

cor_df <- fread("~/Projects/NDUFB7_HF_2026_04_20/03_results/nduf_family_correlation_long.txt")

b_genes <- c("NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11")
ref_genes <- c("NDUFS4","NDUFA9","NDUFV1","NDUFAF2","NDUFC2")
all_genes <- c(b_genes, ref_genes)

cor_df$Gene1 <- factor(cor_df$Gene1, levels=all_genes)
cor_df$Gene2 <- factor(cor_df$Gene2, levels=all_genes)

# 热图
p1 <- ggplot(cor_df, aes(x=Gene2, y=Gene1, fill=Correlation)) +
  geom_tile(color="white", size=0.1) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, limit=c(0,1)) +
  geom_text(aes(label=round(Correlation,2)), size=2.5, 
            color=ifelse(cor_df$Correlation > 0.65, "white", "black")) +
  labs(title="A. NDUF Family Co-expression Network", x="", y="") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8, face="bold"),
        axis.text.y=element_text(size=8, face="bold"),
        panel.grid=element_blank(),
        plot.title=element_text(face="bold", size=10))

# NDUFB7特异性
ndufb7_data <- cor_df[Gene1 == "NDUFB7" & Gene2 != "NDUFB7"]
ndufb7_data$Type <- ifelse(ndufb7_data$Gene2 %in% b_genes, "B-family", "Reference")
ndufb7_data$Gene2 <- factor(ndufb7_data$Gene2, levels=ndufb7_data$Gene2[order(ndufb7_data$Correlation, decreasing=TRUE)])

p2 <- ggplot(ndufb7_data, aes(x=Gene2, y=Correlation, fill=Type)) +
  geom_bar(stat="identity", width=0.7, alpha=0.85) +
  geom_hline(yintercept=median(ndufb7_data[Type=="B-family"]$Correlation), linetype="dashed", color="#E41A1C") +
  scale_fill_manual(values=c("B-family"="#E41A1C", "Reference"="#377EB8")) +
  labs(title="B. NDUFB7 Co-expression Specificity", x="", y="Pearson r") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1, face="bold"), legend.position="top")

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"
ggsave(file.path(out_dir, "Figure3A_heatmap.png"), p1, width=10, height=9, dpi=300)
ggsave(file.path(out_dir, "Figure3B_specificity.png"), p2, width=8, height=5, dpi=300)
cat("✅ Figure 3已保存\n")
