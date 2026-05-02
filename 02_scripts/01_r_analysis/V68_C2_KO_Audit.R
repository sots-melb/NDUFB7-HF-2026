ko <- read.csv("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_NDUFB7_Coexpression_KO.csv", stringsAsFactors=FALSE)
message("▶ 总基因数: ", nrow(ko))
message("▶ 正相关(KO下调): ", sum(ko$KO_Expected=="Down_regulated"))
message("▶ 负相关(KO上调): ", sum(ko$KO_Expected=="Up_regulated"))
message("▶ 无变化: ", sum(ko$KO_Expected=="Unchanged"))

# Top正相关（KO后最可能下调）
top_pos <- head(ko[ko$KO_Expected=="Down_regulated", ], 20)
message("▶ Top 20 KO后预期下调基因:")
print(top_pos[, c("Gene","Correlation")])

# Top负相关（KO后最可能上调）——如果存在
top_neg <- head(ko[ko$KO_Expected=="Up_regulated", ], 20)
if(nrow(top_neg)>0) {
  message("▶ Top 20 KO后预期上调基因:")
  print(top_neg[, c("Gene","Correlation")])
} else {
  message("⚠️ 无显著负相关基因（|r|>0.1）")
}

# 铁死亡基因子集
ferro <- c("FTL","FTH1","SLC7A11","GPX4","SAT1","ACSL4","NFE2L2","KEAP1","LPCAT3","PTGS2","ALOX15","GLS2")
ferro_df <- ko[ko$Gene %in% ferro, ]
message("▶ 铁死亡基因完整表:")
print(ferro_df[order(-ferro_df$Correlation), ])

# 通路富集模拟（基于Top 200 KO下调基因）
top200 <- head(ko[ko$KO_Expected=="Down_regulated", ], 200)
message("▶ Top 200下调基因中已知线粒体基因:")
mito_genes <- grep("^NDU|^COX|^ATP|^UQCR|^SDH|^CYB|^MT-", top200$Gene, value=TRUE)
message("  命中: ", paste(mito_genes, collapse=", "))
