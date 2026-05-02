library(data.table)
library(ggplot2)

eqtlgen <- data.table(
  Dataset = "eQTLGen\n(Whole Blood, n=31,684)",
  Beta = 0.15, SE = 0.025, Pval = 2.6e-9,
  Direction = "A ↑ NDUFB7", SNP = "rs11085898"
)
gtex <- data.table(
  Dataset = "GTEx v11\n(Heart LV, n=387)",
  Beta = -0.0803, SE = 0.0166, Pval = 6.6e-6,
  Direction = "T ↓ NDUFB7", SNP = "rs8103021"
)
combined <- rbind(eqtlgen, gtex)

pA <- ggplot(combined, aes(x=Beta, y=Dataset, color=Direction)) +
  geom_vline(xintercept=0, linetype="dashed", color="gray50") +
  geom_errorbarh(aes(xmin=Beta-1.96*SE, xmax=Beta+1.96*SE), height=0.2, size=1) +
  geom_point(size=4) +
  geom_text(aes(label=paste0("p=", format(Pval, digits=2))), hjust=-0.3, size=3.5) +
  scale_color_manual(values=c("A ↑ NDUFB7"="#E41A1C", "T ↓ NDUFB7"="#377EB8")) +
  labs(title="A. Tissue-Specific eQTL Effect Directions", x="Effect Size (slope/beta)", y="") +
  theme_bw() + theme(legend.position="top", plot.title=element_text(face="bold", size=11))

snp_pos <- data.table(SNP=c("rs11085898","rs8103021"), Pos=c(14679428,14641837), Dataset=c("eQTLGen","GTEx"))
pB <- ggplot(snp_pos, aes(x=Pos, y=Dataset, color=Dataset)) +
  geom_segment(aes(x=14641837, xend=14679428, y="GTEx", yend="eQTLGen"), color="gray50", linetype="dashed") +
  geom_point(size=5) + geom_text(aes(label=SNP), vjust=-1, size=3.5) +
  scale_color_manual(values=c("eQTLGen"="#E41A1C", "GTEx"="#377EB8")) +
  labs(title="B. Lead SNP Positions on chr19", x="Genomic Position (bp)", y="") +
  theme_bw() + theme(legend.position="none")

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"
ggsave(file.path(out_dir, "Figure4A_forest.png"), pA, width=7, height=4, dpi=300)
ggsave(file.path(out_dir, "Figure4B_position.png"), pB, width=7, height=3, dpi=300)
cat("✅ Figure 4已保存\n")
