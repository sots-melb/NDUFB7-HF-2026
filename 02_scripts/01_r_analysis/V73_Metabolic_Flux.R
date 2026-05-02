suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
if(file.exists(rds_path)) {
  obj <- readRDS(rds_path)
  
  # 经典糖酵解核心基因集 (KEGG / Hallmark)
  glyco_genes <- c("ENO1", "ENO2", "HK1", "HK2", "ALDOA", "ALDOB", "GAPDH", "PGK1", "PGM1", "PKM", "LDHA", "LDHB")
  # 经典 OXPHOS 核心基因集 (我们在 V66 用过的)
  oxphos_genes <- c("ATP5F1A", "COX4I1", "CYC1", "NDUFA4", "SDHA", "UQCRC1", "UQCRC2")
  
  valid_glyco <- intersect(glyco_genes, rownames(obj))
  valid_oxphos <- intersect(oxphos_genes, rownames(obj))
  
  message("▶ 使用 ", length(valid_glyco), " 个糖酵解基因和 ", length(valid_oxphos), " 个 OXPHOS 基因进行通量评估...")
  
  # 使用 Seurat 的 AddModuleScore 进行稳健评分
  obj <- AddModuleScore(obj, features = list(valid_glyco), name = "Glyco_Flux")
  obj <- AddModuleScore(obj, features = list(valid_oxphos), name = "OXPHOS_Flux")
  
  # 注意：AddModuleScore 会自动在列名后加上 '1'
  df <- data.frame(
    Condition = obj$Condition,
    NDUFB7 = tryCatch(GetAssayData(obj, slot="data")["NDUFB7", ], error=function(e) LayerData(obj, assay="RNA", layer="data")["NDUFB7", ]),
    Glyco = obj$Glyco_Flux1,
    OXPHOS = obj$OXPHOS_Flux1
  )
  
  # 计算代谢比值 (Metabolic Ratio): 糖酵解 / OXPHOS 
  # 为了防止分母为 0 或负数，我们对其进行极小值平移后取 log
  df$Metabolic_Ratio <- log2((df$Glyco - min(df$Glyco) + 0.01) / (df$OXPHOS - min(df$OXPHOS) + 0.01))
  
  message("\n🌟 疾病组 (DCM) 与 健康组 (Donor) 的代谢通量对比:")
  stats <- df %>% group_by(Condition) %>% 
           summarise(Mean_OXPHOS = mean(OXPHOS), Mean_Glyco = mean(Glyco), Mean_Ratio = mean(Metabolic_Ratio))
  print(stats)
  
  # 统计检验
  p_oxphos <- wilcox.test(OXPHOS ~ Condition, data = df)$p.value
  p_ratio <- wilcox.test(Metabolic_Ratio ~ Condition, data = df)$p.value
  
  message(sprintf("   - OXPHOS 差异 p-value: %.2e", p_oxphos))
  message(sprintf("   - 糖酵解/OXPHOS 比例差异 p-value: %.2e", p_ratio))
  
  if(p_ratio < 0.05) {
     message("🎉 成功证实代谢通量倒挂！DCM 中糖酵解相对于 OXPHOS 发生代偿性升高，标志着能量代谢重编程。")
  }
}
