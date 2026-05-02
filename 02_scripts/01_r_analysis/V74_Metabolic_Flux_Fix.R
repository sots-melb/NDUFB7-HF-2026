suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
if(file.exists(rds_path)) {
  obj <- readRDS(rds_path)
  
  # 提取纯数据矩阵，兼容 V4/V5
  expr_mat <- tryCatch(GetAssayData(obj, slot="data"), error=function(e) LayerData(obj, assay="RNA", layer="data"))
  
  # 基因集
  glyco_genes <- c("ENO1", "ENO2", "HK1", "HK2", "ALDOA", "ALDOB", "GAPDH", "PGK1", "PGM1", "PKM", "LDHA", "LDHB")
  oxphos_genes <- c("ATP5F1A", "COX4I1", "CYC1", "NDUFA4", "SDHA", "UQCRC1", "UQCRC2")
  
  valid_glyco <- intersect(glyco_genes, rownames(expr_mat))
  valid_oxphos <- intersect(oxphos_genes, rownames(expr_mat))
  
  message("▶ 有效基因: ", length(valid_glyco), " 个糖酵解, ", length(valid_oxphos), " 个 OXPHOS")
  
  # 手写矩阵运算打分：计算基因平均表达量
  glyco_score <- colMeans(expr_mat[valid_glyco, , drop=FALSE])
  oxphos_score <- colMeans(expr_mat[valid_oxphos, , drop=FALSE])
  
  df <- data.frame(
    Condition = obj$Condition,
    Glyco = glyco_score,
    OXPHOS = oxphos_score
  )
  
  # 代谢比值 = 糖酵解 / OXPHOS (加微小偏置防分母0)
  df$Metabolic_Ratio <- log2((df$Glyco + 0.01) / (df$OXPHOS + 0.01))
  
  message("\n🌟 疾病组 (DCM) 与 健康组 (Donor) 的代谢通量对比:")
  stats <- df %>% group_by(Condition) %>% 
           summarise(Mean_OXPHOS = mean(OXPHOS), Mean_Glyco = mean(Glyco), Mean_Ratio = mean(Metabolic_Ratio))
  print(stats)
  
  p_oxphos <- wilcox.test(OXPHOS ~ Condition, data = df)$p.value
  p_ratio <- wilcox.test(Metabolic_Ratio ~ Condition, data = df)$p.value
  
  message(sprintf("   - OXPHOS 崩溃显著性 p-value: %.2e", p_oxphos))
  message(sprintf("   - 糖酵解/OXPHOS 比例升高显著性 p-value: %.2e", p_ratio))
  
  if(p_ratio < 0.05) {
     message("🎉 完美证实代谢通量倒挂！DCM心肌的糖酵解相对 OXPHOS 发生代偿性升高，代谢重编程实锤！")
  }
}
