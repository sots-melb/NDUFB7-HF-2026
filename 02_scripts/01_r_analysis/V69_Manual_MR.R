suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

data_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/smr_input"
gwas_file <- file.path(data_dir, "HERMES_HF_SMR.ma")
eqtl_file <- file.path(data_dir, "eQTLGen_NDUFB7.ma")

if(file.exists(gwas_file) && file.exists(eqtl_file)) {
  gwas <- read.table(gwas_file, header=TRUE, stringsAsFactors=FALSE)
  eqtl <- read.table(eqtl_file, header=TRUE, stringsAsFactors=FALSE)
  
  merged <- inner_join(eqtl, gwas, by="SNP", suffix=c("_exp", "_out"))
  
  if(nrow(merged) > 0) {
    # 对齐方向并计算 Wald Ratio
    merged <- merged %>% 
      mutate(b_out_aligned = ifelse(A1_exp == A1_out, b_out, ifelse(A1_exp == A2_out, -b_out, NA))) %>%
      filter(!is.na(b_out_aligned))
      
    top <- merged[which.min(merged$p_exp), ]
    wr_beta <- top$b_out_aligned / top$b_exp
    wr_se <- abs(wr_beta) * sqrt((top$se_out/top$b_out_aligned)^2 + (top$se_exp/top$b_exp)^2)
    wr_p <- 2 * pnorm(abs(wr_beta)/wr_se, lower.tail=FALSE)
    
    or <- exp(wr_beta)
    lci <- exp(wr_beta - 1.96*wr_se)
    uci <- exp(wr_beta + 1.96*wr_se)
    
    message(sprintf("✅ 手动 MR 计算成功！OR: %.3f (95%% CI: %.3f - %.3f), p-value: %.2e", or, lci, uci, wr_p))
  } else {
    message("⚠️ 未找到重叠的 SNP。")
  }
} else {
  message("⚠️ 缺少 .ma 数据文件。")
}
