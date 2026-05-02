suppressMessages(library(WGCNA))
csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE57338_expression_matrix_annotated.csv"
if(file.exists(csv)) {
  df <- read.csv(csv, stringsAsFactors=FALSE, check.names=FALSE)
  genes <- df[,1]
  
  # 强制转换为数值矩阵，忽略无法转换的字符并转为 NA
  expr <- suppressWarnings(apply(df[, -1], 2, as.numeric))
  rownames(expr) <- genes
  
  # 去除全 NA 的行或列
  expr <- expr[rowSums(is.na(expr)) < ncol(expr)*0.5, colSums(is.na(expr)) < nrow(expr)*0.5]
  
  datExpr <- t(expr)
  message("▶ 强制数值化完成！WGCNA 输入维度: ", ncol(datExpr), " 样本 x ", nrow(datExpr), " 基因")
  
  gsg <- goodSamplesGenes(datExpr, verbose=0)
  datExpr <- datExpr[, gsg$goodGenes]
  message("▶ 质控过滤后基因数: ", ncol(datExpr))
  
  # 快速构建网络 (使用预设软阈值 6 加快速度)
  net <- blockwiseModules(datExpr, power=6, TOMType="unsigned", minModuleSize=30, numericLabels=TRUE, verbose=0)
  
  if("NDUFB7" %in% colnames(datExpr)) {
    mod <- net$colors["NDUFB7"]
    message("✅ WGCNA 修复成功！NDUFB7 位于模块: ", mod, " (包含基因数: ", sum(net$colors == mod), ")")
  }
}
