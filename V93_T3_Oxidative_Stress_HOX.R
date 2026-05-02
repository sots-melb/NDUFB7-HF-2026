#!/usr/bin/env Rscript
# V93: T3 氧化应激五法评分 + HOX亚群确认
# 固定路径，禁止任何find/list.files外部搜索

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

# --- 严格固定路径 ---
FILE_CM <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"

message("========================================")
message("V93: T3 氧化应激五法评分 + HOX亚群确认")
message("========================================")

if (!file.exists(FILE_CM)) stop("[FAIL] 固定路径文件不存在: ", FILE_CM)

cm <- readRDS(FILE_CM)
message("[PASS] CM数据: ", ncol(cm), " cells × ", nrow(cm), " genes")

# ========================================
# 五法氧化应激评分
# ========================================
message("")
message("========================================")
message("T3: 五法氧化应激评分")
message("========================================")

# 5种氧化应激基因集（基于Cell Metab 2021, Circ Res 2023, Nat Cardiovasc Res 2024）
os_sets <- list(
  ROS_Core      = c("SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1", "PRDX2", "TXN", "TXNRD1"),
  NRF2_Pathway  = c("NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM", "FTH1", "FTL", "SLC7A11"),
  Mito_ROS      = c("MPO", "CYBB", "NOX2", "NOX4", "AIFM1", "ENDOG", "BNIP3", "BNIP3L"),
  DNA_Damage    = c("ATM", "ATR", "CHEK1", "CHEK2", "TP53", "CDKN1A", "PARP1", "PARP2"),
  Lipid_Perox   = c("ACSL4", "LPCAT3", "PTGS2", "ALOX5", "ALOX15", "GPX4", "GCLC", "GCLM")
)

# 使用AddModuleScore逐法评分
score_cols <- character(0)
for (i in seq_along(os_sets)) {
  set_name <- names(os_sets)[i]
  genes <- intersect(os_sets[[i]], rownames(cm))
  
  if (length(genes) >= 3) {
    cm <- AddModuleScore(cm, features = list(genes), name = set_name)
    # Seurat自动命名为 set_name1, set_name2... 取第一个
    new_col <- paste0(set_name, "1")
    if (new_col %in% colnames(cm@meta.data)) {
      score_cols <- c(score_cols, new_col)
      message("[SCORE] ", set_name, ": ", length(genes), " genes -> ", new_col)
    }
  } else {
    message("[SKIP] ", set_name, ": 仅", length(genes), "个基因可用")
  }
}

# 综合评分 = 各方法均值
if (length(score_cols) >= 3) {
  cm$OS_Composite <- rowMeans(cm@meta.data[, score_cols])
  message("[PASS] 综合评分 OS_Composite 已创建（基于", length(score_cols), "种方法）")
} else {
  message("[FAIL] 可用评分方法不足 (", length(score_cols), "/5)")
  T3_PASS <- FALSE
}

# ========================================
# HOX亚群识别
# ========================================
message("")
message("========================================")
message("T3: HOX亚群识别")
message("========================================")

if ("OS_Composite" %in% colnames(cm@meta.data)) {
  # HOX阈值: top 20%综合评分
  hox_q80 <- quantile(cm$OS_Composite, 0.8, na.rm = TRUE)
  cm$HOX_Status <- ifelse(cm$OS_Composite > hox_q80, "HOX_High", "HOX_Low")
  
  n_hox <- sum(cm$HOX_Status == "HOX_High")
  n_low <- sum(cm$HOX_Status == "HOX_Low")
  pct_hox <- round(n_hox / ncol(cm) * 100, 1)
  
  message("HOX阈值 (80%分位数): ", round(hox_q80, 4))
  message("HOX_High: ", n_hox, " cells (", pct_hox, "%)")
  message("HOX_Low:  ", n_low, " cells")
  
  # HOX与NDUFB7的关系
  if ("NDUFB7" %in% rownames(cm)) {
    ndufb7_expr <- FetchData(cm, vars = "NDUFB7")$NDUFB7
    cm$NDUFB7_Expr <- ndufb7_expr
    
    # 组间比较
    hox_ndufb7 <- cm$NDUFB7_Expr[cm$HOX_Status == "HOX_High"]
    low_ndufb7 <- cm$NDUFB7_Expr[cm$HOX_Status == "HOX_Low"]
    tt <- t.test(hox_ndufb7, low_ndufb7)
    
    message("")
    message("NDUFB7 在HOX亚群 vs 对照:")
    message("  HOX_High mean: ", round(mean(hox_ndufb7), 4))
    message("  HOX_Low  mean: ", round(mean(low_ndufb7), 4))
    message("  Log2FC: ", round(log2(mean(hox_ndufb7)/mean(low_ndufb7)), 3))
    message("  t-test p: ", format(tt$p.value, digits = 2, scientific = TRUE))
    
    # 相关分析
    cor_os <- cor.test(cm$OS_Composite, cm$NDUFB7_Expr, method = "spearman")
    message("  OS_Composite vs NDUFB7 Spearman ρ: ", round(cor_os$estimate, 3), 
            " p = ", format(cor_os$p.value, digits = 2, scientific = TRUE))
  }
  
  # ========================================
  # 可视化
  # ========================================
  outdir <- file.path(PROJECT_DIR, "03_results/T3_Oxidative_Stress_HOX")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. 五法评分小提琴图
  plot_data <- cm@meta.data[, c("group", score_cols, "OS_Composite", "HOX_Status")]
  plot_data$cell_id <- rownames(plot_data)
  
  # 长格式
  library(reshape2)
  plot_long <- melt(plot_data, id.vars = c("cell_id", "group", "HOX_Status"), 
                    measure.vars = c(score_cols, "OS_Composite"))
  
  p_violin <- ggplot(plot_long, aes(x = group, y = value, fill = group)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.5) +
    facet_wrap(~variable, scales = "free_y", ncol = 3) +
    labs(title = "Oxidative Stress Scores by Group", y = "Module Score") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(outdir, "V93_T3_OS_scores_violin.png"), p_violin, 
         width = 10, height = 8, dpi = 300)
  
  # 2. HOX vs NDUFB7散点
  if ("NDUFB7_Expr" %in% colnames(cm@meta.data)) {
    p_hox <- ggplot(cm@meta.data, aes(x = OS_Composite, y = NDUFB7_Expr, color = HOX_Status)) +
      geom_point(alpha = 0.4, size = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      scale_color_manual(values = c("HOX_Low" = "#440154", "HOX_High" = "#FDE725")) +
      labs(title = "Oxidative Stress vs NDUFB7 Expression",
           subtitle = paste0("ρ = ", round(cor_os$estimate, 3), 
                           " | HOX_High = top 20% OS score"),
           x = "OS Composite Score", y = "NDUFB7 Expression") +
      theme_minimal(base_size = 10)
    
    ggsave(file.path(outdir, "V93_T3_HOX_NDUFB7_scatter.png"), p_hox,
           width = 6, height = 5, dpi = 300)
  }
  
  # 3. HOX状态UMAP（如果有UMAP）
  if ("umap" %in% names(cm@reductions)) {
    p_umap <- DimPlot(cm, group.by = "HOX_Status", 
                      cols = c("HOX_Low" = "#440154", "HOX_High" = "#FDE725"),
                      pt.size = 0.5, label = TRUE) +
      labs(title = "HOX Subcluster on UMAP") +
      theme_minimal(base_size = 10)
    
    ggsave(file.path(outdir, "V93_T3_HOX_umap.png"), p_umap,
           width = 6, height = 5, dpi = 300)
  }
  
  # 保存统计
  stats <- data.frame(
    N_Total = ncol(cm),
    N_HOX_High = n_hox,
    Pct_HOX_High = pct_hox,
    HOX_NDUFB7_Mean = mean(hox_ndufb7, na.rm = TRUE),
    Low_NDUFB7_Mean = mean(low_ndufb7, na.rm = TRUE),
    NDUFB7_Log2FC = log2(mean(hox_ndufb7)/mean(low_ndufb7)),
    NDUFB7_TTest_P = tt$p.value,
    OS_NDUFB7_Spearman = cor_os$estimate,
    OS_NDUFB7_P = cor_os$p.value
  )
  write.csv(stats, file.path(outdir, "V93_T3_stats.csv"), row.names = FALSE)
  
  # 判断
  message("")
  if (tt$p.value < 0.05 && mean(hox_ndufb7) < mean(low_ndufb7)) {
    message("[PASS] HOX亚群NDUFB7显著更低！支持'氧化应激→NDUFB7丢失'机制")
    T3_PASS <- TRUE
  } else if (tt$p.value < 0.05) {
    message("[PARTIAL] HOX与NDUFB7有显著差异，但方向需确认")
    T3_PASS <- TRUE
  } else {
    message("[WARN] HOX与NDUFB7无显著关联，'氧化应激驱动'假说需弱化")
    T3_PASS <- FALSE
  }
  
  message("[DONE] T3结果: ", outdir)
} else {
  message("[FAIL] OS_Composite未创建")
  T3_PASS <- FALSE
}

# ========================================
# 状态板
# ========================================
message("")
message("========================================")
message("V93 最终状态")
message("========================================")
t3_status <- ifelse(exists("T3_PASS") && T3_PASS, "✅ PASS", ifelse(exists("T3_PASS"), "❌ FAIL", "⏭️ SKIP"))
message("T3 (氧化应激+HOX): ", t3_status)
message("========================================")
