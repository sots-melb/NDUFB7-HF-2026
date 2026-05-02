#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     hdWGCNA 线粒体共表达模块分析                         ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# --- 0. 安装hdWGCNA（若未安装）---
cat("========== 0. 环境检查 ==========\n")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
}

if (!requireNamespace("hdWGCNA", quietly = TRUE)) {
  cat("⏳ hdWGCNA未安装，正在安装（约5-10分钟）...\n")
  tryCatch({
    remotes::install_github("smorabit/hdWGCNA", ref = "dev", quiet = TRUE)
    cat("✅ hdWGCNA安装成功\n")
  }, error = function(e) {
    cat("❌ 自动安装失败:", conditionMessage(e), "\n")
    cat("请手动运行: remotes::install_github('smorabit/hdWGCNA', ref='dev')\n")
    quit(status = 1)
  })
}

library(hdWGCNA)
cat("✅ hdWGCNA加载成功，版本:", packageVersion("hdWGCNA"), "\n")

# --- 1. 加载数据 ---
cat("\n========== 1. 加载Seurat对象 ==========\n")
srt <- readRDS("03_results/09_single_cell/22_srt_with_pseudotime.rds")
cat("✅ 对象:", ncol(srt), "cells x", nrow(srt), "genes\n")

# --- 2. 定义线粒体基因集 ---
cat("\n========== 2. 定义线粒体基因集 ==========\n")

# 复合体I核心亚基（N模块 + Q模块 + B模块）
cI_n <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NDUFAB1")
cI_q <- c("NDUFA1","NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10",
          "NDUFA11","NDUFA12","NDUFA13","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5",
          "NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11")
cI_core <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6")

# 复合体II-V标志
cII <- c("SDHA","SDHB","SDHC","SDHD")
cIII <- c("UQCRC1","UQCRC2","CYC1","UQCRB","UQCRQ","UQCR10","UQCR11")
cIV <- c("COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A1","COX7B","COX8A")
cV <- c("ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5MC1","ATP5MC2","ATP5MC3")
mito_biogenesis <- c("PPARGC1A","PGC1A","NRF1","NRF2","TFAM","TFB1M","TFB2M")
mito_dynamics <- c("MFN1","MFN2","OPA1","DNM1L","FIS1")

all_mito <- unique(c(cI_n, cI_q, cI_core, cII, cIII, cIV, cV, mito_biogenesis, mito_dynamics))
cat("定义线粒体基因集:", length(all_mito), "个\n")

# 检查可用性
avail_mito <- all_mito[all_mito %in% rownames(srt)]
cat("数据中可用:", length(avail_mito), "个\n")
cat("包含NDUFB7:", "NDUFB7" %in% avail_mito, "\n")
cat("包含PGC-1α(PPARGC1A):", "PPARGC1A" %in% avail_mito, "\n")
cat("包含NRF1:", "NRF1" %in% avail_mito, "\n")
cat("包含TFAM:", "TFAM" %in% avail_mito, "\n")

if (length(avail_mito) < 30) {
  cat("⚠️ 线粒体基因不足30个，扩展至全基因组但聚焦线粒体功能注释\n")
  # 备用：使用全基因组，后续通过模块富集筛选
  avail_mito <- rownames(srt)
}

# --- 3. hdWGCNA设置 ---
cat("\n========== 3. hdWGCNA网络构建 ==========\n")

srt <- SetupForWGCNA(
  srt,
  gene_select = "custom",
  gene_list = avail_mito,
  wgcna_name = "mito_modules"
)
cat("✅ Setup完成\n")

# 按condition构建metacells（减少噪声）
srt <- MetacellsByGroups(
  srt,
  group.by = "condition",
  reduction = "umap",
  k = 15,  # 503 cells较小，k不宜过大
  max_shared = 5,
  ident.group = "condition"
)
cat("✅ Metacells构建完成\n")

srt <- SetDatExpr(
  srt,
  group_name = c("HF", "Control"),
  group.by = "condition",
  assay = "RNA",
  slot = "data"
)
cat("✅ 表达矩阵设置完成\n")

# 测试软阈值
cat("\n========== 4. 测试软阈值 ==========\n")
srt <- TestSoftPowers(srt, networkType = "signed")
cat("✅ 软阈值测试完成\n")

# 构建网络（核心步骤，约10-20分钟）
cat("\n========== 5. 构建共表达网络（约10-20分钟）==========\n")
srt <- ConstructNetwork(
  srt,
  soft_power = 6,  # 若TestSoftPowers建议不同，可调整
  setDatExpr = FALSE,
  tom_outdir = "03_results/09_single_cell/hdwgcna",
  tom_name = "mito_tom"
)
cat("✅ 网络构建完成\n")

# --- 6. 模块检测与特征计算 ---
cat("\n========== 6. 模块检测 ==========\n")
srt <- ModuleEigengenes(srt)
cat("✅ 模块特征基因计算完成\n")

srt <- ModuleConnectivity(srt)
cat("✅ 模块连接度计算完成\n")

# 获取模块分配
modules <- GetModules(srt)
cat("检测到的模块数:", length(unique(modules$module)) - 1, "\n")  # 排除grey
cat("模块大小统计:\n")
print(table(modules$module))

# --- 7. NDUFB7模块归属 ---
cat("\n========== 7. NDUFB7模块归属 ==========\n")
if ("NDUFB7" %in% modules$gene_name) {
  ndufb7_module <- modules$module[modules$gene_name == "NDUFB7"]
  cat("🎯 NDUFB7所属模块:", ndufb7_module, "\n")
  
  # 同模块基因
  same_module <- modules$gene_name[modules$module == ndufb7_module]
  cat("同模块基因数:", length(same_module), "\n")
  cat("同模块基因（前20）:", paste(head(same_module, 20), collapse = ", "), "\n")
  
  # 检查PGC-1α/NRF1/TFAM是否同模块
  buddies <- c("PPARGC1A","PGC1A","NRF1","NRF2","TFAM","NDUFB8","NDUFB10")
  for (g in buddies) {
    if (g %in% same_module) {
      cat("✅", g, "与NDUFB7同模块！\n")
    }
  }
  
  write.csv(data.frame(
    gene = same_module,
    module = ndufb7_module,
    kME = modules$kME[modules$module == ndufb7_module]
  ), "03_results/09_single_cell/hdwgcna/23_ndufb7_module_genes.csv", row.names = FALSE)
  cat("[保存] 23_ndufb7_module_genes.csv\n")
} else {
  cat("⚠️ NDUFB7未进入任何模块（可能在grey模块）\n")
}

# --- 8. 模块与Condition关联 ---
cat("\n========== 8. 模块与Condition关联 ==========\n")
srt <- ModuleExprScore(srt, n_genes = 25)
srt <- SetModuleScore(srt, n_genes = 25)

# 检验每个模块与condition的关联
me_df <- GetMEs(srt)
me_df$condition <- srt$condition[match(rownames(me_df), colnames(srt))]

module_stats <- data.frame()
for (mod in colnames(me_df)[grepl("^ME", colnames(me_df))]) {
  wt <- wilcox.test(me_df[[mod]] ~ me_df$condition)
  hf_mean <- mean(me_df[[mod]][me_df$condition == "HF"], na.rm = TRUE)
  ctrl_mean <- mean(me_df[[mod]][me_df$condition == "Control"], na.rm = TRUE)
  module_stats <- rbind(module_stats, data.frame(
    module = mod,
    HF_mean = round(hf_mean, 4),
    Control_mean = round(ctrl_mean, 4),
    log2FC = round(log2((hf_mean + 0.001) / (ctrl_mean + 0.001)), 4),
    p_value = format(wt$p.value, digits = 4, scientific = TRUE),
    significant = wt$p.value < 0.05
  ))
}
module_stats <- module_stats[order(module_stats$p_value), ]
print(module_stats)
write.csv(module_stats, "03_results/09_single_cell/hdwgcna/24_module_condition_stats.csv", row.names = FALSE)
cat("[保存] 24_module_condition_stats.csv\n")

# --- 9. 可视化 ---
cat("\n========== 9. 可视化 ==========\n")

# 模块UMAP
tryCatch({
  srt <- RunModuleUMAP(srt, n_hubs = 10)
  p_umap <- ModuleUMAPPlot(srt, label_hubs = TRUE)
  ggsave("03_results/09_single_cell/hdwgcna/25_module_umap.pdf", p_umap, width = 10, height = 8)
  cat("[保存] 25_module_umap.pdf\n")
}, error = function(e) {
  cat("⚠️ ModuleUMAP失败:", conditionMessage(e), "\n")
})

# 模块特征基因热图
tryCatch({
  p_me <- PlotModuleEigengenes(srt, features = "condition")
  ggsave("03_results/09_single_cell/hdwgcna/26_module_eigengenes.pdf", p_me, width = 10, height = 8)
  cat("[保存] 26_module_eigengenes.pdf\n")
}, error = function(e) {
  cat("⚠️ Eigengenes热图失败:", conditionMessage(e), "\n")
})

# Hub基因网络（NDUFB7模块）
if ("NDUFB7" %in% modules$gene_name) {
  ndufb7_mod <- modules$module[modules$gene_name == "NDUFB7"]
  tryCatch({
    p_hub <- HubGeneNetworkPlot(srt, mods = ndufb7_mod)
    ggsave("03_results/09_single_cell/hdwgcna/27_ndufb7_hub_network.pdf", p_hub, width = 10, height = 10)
    cat("[保存] 27_ndufb7_hub_network.pdf\n")
  }, error = function(e) {
    cat("⚠️ Hub网络失败:", conditionMessage(e), "\n")
  })
}

# --- 10. 保存 ---
cat("\n========== 10. 保存结果 ==========\n")
saveRDS(srt, "03_results/09_single_cell/hdwgcna/28_srt_hdwgcna.rds")
saveRDS(modules, "03_results/09_single_cell/hdwgcna/29_modules_table.rds")
cat("[保存] 28_srt_hdwgcna.rds + 29_modules_table.rds\n")

cat("\n🎉 hdWGCNA分析完成！\n")
cat("\n【关键产出】\n")
cat("  23_ndufb7_module_genes.csv —— NDUFB7同模块基因列表\n")
cat("  24_module_condition_stats.csv —— 模块与HF/Control关联\n")
cat("  25_module_umap.pdf —— 模块UMAP\n")
cat("  27_ndufb7_hub_network.pdf —— NDUFB7模块hub网络\n")
