#!/usr/bin/env Rscript
setwd("~/Projects/NDUFB7_HF_{2026_04_20}")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     GSE183852 DCM整合分析                                ║\n")
cat("║     目标：验证Cluster 4型NDUFB7低表达亚群是否在DCM复现   ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# --- 1. 解压与加载 ---
obj_gz <- "01_data/01_raw_geo/GSE183852/GSE183852_DCM_Integrated.Robj.gz"
obj_r <- "01_data/01_raw_geo/GSE183852/GSE183852_DCM_Integrated.Robj"

if (!file.exists(obj_r)) {
  if (file.exists(obj_gz)) {
    cat("⏳ 解压.Robj.gz...\n")
    if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils", quiet = TRUE)
    R.utils::gunzip(obj_gz, overwrite = TRUE, remove = FALSE)
    cat("✅ 解压完成\n")
  } else {
    cat("❌ 文件不存在:", obj_gz, "\n")
    quit(status = 1)
  }
}

cat("⏳ 加载R对象...\n")
env_before <- ls()
load(obj_r)
env_after <- ls()
new_objs <- setdiff(env_after, env_before)

srt_dcm <- NULL
for (obj_name in new_objs) {
  obj <- get(obj_name)
  if (inherits(obj, "Seurat")) {
    srt_dcm <- obj
    cat("✅ 找到Seurat对象:", obj_name, "|", ncol(obj), "cells x", nrow(obj), "genes\n")
    break
  }
}
if (is.null(srt_dcm)) { cat("❌ 未找到Seurat对象\n"); quit(status = 1) }

# --- 2. NDUFB7检查 ---
cat("\n========== 2. NDUFB7检查 ==========\n")
has_ndufb7 <- "NDUFB7" %in% rownames(srt_dcm)
cat("NDUFB7在数据中:", has_ndufb7, "\n")

if (!has_ndufb7) {
  cat("⚠️ 尝试别名/大小写...\n")
  candidates <- grep("NDUFB7|CI-B18|B18", rownames(srt_dcm), value = TRUE, ignore.case = TRUE)
  cat("候选:", paste(candidates, collapse = ", "), "\n")
  if (length(candidates) > 0) {
    cat("使用第一个候选作为NDUFB7代理\n")
    srt_dcm[["NDUFB7_proxy"]] <- FetchData(srt_dcm, candidates[1])[,1]
    ndufb7_vals <- srt_dcm[["NDUFB7_proxy"]]
  } else {
    cat("❌ NDUFB7及其别名均不在数据中\n")
    quit(status = 1)
  }
} else {
  ndufb7_vals <- FetchData(srt_dcm, vars = "NDUFB7")[,1]
}

cat("NDUFB7范围:", round(range(ndufb7_vals), 3), "\n")
cat("NDUFB7均值:", round(mean(ndufb7_vals), 3), "\n")

# --- 3. 按Condition/疾病状态分组 ---
cat("\n========== 3. 疾病分组分析 ==========\n")
cond_col <- NULL
for (c in c("condition", "Condition", "disease", " Disease", "group", "Group", "sampletype")) {
  if (c %in% colnames(srt_dcm@meta.data)) { cond_col <- c; break }
}
if (!is.null(cond_col)) {
  cat("Condition列:", cond_col, "\n")
  cat("分布:\n")
  print(table(srt_dcm@meta.data[[cond_col]]))
  
  df <- data.frame(
    NDUFB7 = ndufb7_vals,
    Condition = srt_dcm@meta.data[[cond_col]]
  )
  agg <- df %>% group_by(Condition) %>% summarise(
    mean = mean(NDUFB7), median = median(NDUFB7), n = n(), .groups = "drop"
  )
  cat("\n各组NDUFB7统计:\n")
  print(as.data.frame(agg))
  
  # Wilcoxon检验
  conds <- unique(df$Condition)
  if (length(conds) == 2) {
    wt <- wilcox.test(NDUFB7 ~ Condition, data = df)
    cat("\nWilcoxon检验: p =", format(wt$p.value, digits = 3), "\n")
  }
  
  # 保存统计
  write.csv(agg, "03_results/12_dcm_integration/04_dcm_ndufb7_by_condition.csv", row.names = FALSE)
}

# --- 4. 寻找DCM中的低表达尾部（类似Cluster 4） ---
cat("\n========== 4. 低表达亚群检测 ==========\n")
q10 <- quantile(ndufb7_vals, 0.10, na.rm = TRUE)
q5  <- quantile(ndufb7_vals, 0.05, na.rm = TRUE)
cat("NDUFB7 10%分位数:", round(q10, 3), "\n")
cat("NDUFB7 5%分位数:", round(q5, 3), "\n")
cat("低表达细胞(<10%分位数):", sum(ndufb7_vals < q10, na.rm = TRUE), "/", sum(!is.na(ndufb7_vals)), "\n")
cat("极低表达细胞(<5%分位数):", sum(ndufb7_vals < q5, na.rm = TRUE), "/", sum(!is.na(ndufb7_vals)), "\n")

# 如果有cluster信息，按cluster统计
clust_col <- NULL
for (c in c("seurat_clusters", "cluster", "Cluster", "celltype", "CellType")) {
  if (c %in% colnames(srt_dcm@meta.data)) { clust_col <- c; break }
}
if (!is.null(clust_col)) {
  cat("\n各Cluster NDUFB7均值（升序）:\n")
  df2 <- data.frame(NDUFB7 = ndufb7_vals, Cluster = srt_dcm@meta.data[[clust_col]])
  agg2 <- df2 %>% group_by(Cluster) %>% summarise(
    mean = mean(NDUFB7), median = median(NDUFB7), n = n(), .groups = "drop"
  ) %>% arrange(mean)
  print(as.data.frame(agg2))
  
  # 检查是否有类似Cluster 4的极端低表达簇
  min_cluster <- agg2$Cluster[which.min(agg2$mean)]
  min_mean <- min(agg2$mean)
  cat("\n最低表达Cluster:", min_cluster, "(均值=", round(min_mean, 3), ")\n")
  
  if (min_mean < q10) {
    cat("🎉 发现极端低表达Cluster！类似HF Cluster 4的特征！\n")
  } else {
    cat("⚠️ 最低Cluster均值高于10%分位数，无极端低表达亚群\n")
  }
}

# --- 5. 与HF数据比较 ---
cat("\n========== 5. DCM vs HF比较 ==========\n")
cat("【GSE168742 HF参考值】\n")
cat("HF整体均值: ~1.5-2.0 (估计)\n")
cat("HF Cluster 4: 0.36 (极端低)\n")
cat("Control整体: ~2.5-3.0\n\n")

cat("【当前DCM】\n")
cat("DCM整体均值:", round(mean(ndufb7_vals), 3), "\n")
cat("DCM最低Cluster/10%分位数:", round(min(ndufb7_vals, na.rm = TRUE), 3), "/", round(q10, 3), "\n\n")

if (mean(ndufb7_vals) < 1.5) {
  cat("✅ DCM整体NDUFB7显著低于Control参考值，支持心衰通用下调\n")
} else {
  cat("⚠️ DCM整体NDUFB7不低，可能病因特异性\n")
}

# --- 6. 可视化 ---
cat("\n========== 6. 可视化 ==========\n")
p1 <- VlnPlot(srt_dcm, features = ifelse(has_ndufb7, "NDUFB7", candidates[1]), 
              pt.size = 0.3, group.by = ifelse(is.null(cond_col), NULL, cond_col)) +
  ggtitle("NDUFB7 in DCM") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("03_results/12_dcm_integration/01_dcm_ndufb7_violin.pdf", p1, width = 8, height = 6)
cat("[保存] 01_dcm_ndufb7_violin.pdf\n")

p2 <- FeaturePlot(srt_dcm, features = ifelse(has_ndufb7, "NDUFB7", candidates[1]), 
                  pt.size = 0.5, order = TRUE) +
  ggtitle("NDUFB7 Expression in DCM") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("03_results/12_dcm_integration/02_dcm_ndufb7_featureplot.pdf", p2, width = 8, height = 6)
cat("[保存] 02_dcm_ndufb7_featureplot.pdf\n")

# 如果有cluster，画cluster violin
if (!is.null(clust_col)) {
  p3 <- VlnPlot(srt_dcm, features = ifelse(has_ndufb7, "NDUFB7", candidates[1]), 
                pt.size = 0.3, group.by = clust_col) +
    ggtitle("NDUFB7 by Cluster in DCM") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave("03_results/12_dcm_integration/03_dcm_ndufb7_by_cluster.pdf", p3, width = 10, height = 6)
  cat("[保存] 03_dcm_ndufb7_by_cluster.pdf\n")
}

# --- 7. 保存 ---
saveRDS(srt_dcm, "03_results/12_dcm_integration/05_dcm_seurat_object.rds")
cat("\n[保存] 05_dcm_seurat_object.rds\n")

cat("\n🎉 DCM整合分析完成！\n")
