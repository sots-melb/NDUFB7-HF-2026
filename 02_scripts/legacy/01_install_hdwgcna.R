#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     hdWGCNA 安装与模块分析（重试版）                     ║\n")
cat("║     时间: 2026-04-21 16:50（Rate Limit已重置）           ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# --- 0. 安装依赖 ---
cat("========== 0. 安装依赖 ==========\n")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
}

# hdWGCNA依赖WGCNA等包
tryCatch({
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    BiocManager::install("WGCNA", ask = FALSE, update = FALSE)
  }
}, error = function(e) {
  cat("⚠️ WGCNA安装警告:", conditionMessage(e), "\n")
})

# --- 1. 安装hdWGCNA ---
cat("\n========== 1. 安装hdWGCNA ==========\n")
install_success <- FALSE
max_try <- 3

for (i in 1:max_try) {
  cat("⏳ 第", i, "次尝试安装 hdWGCNA...\n")
  tryCatch({
    remotes::install_github("smorabit/hdWGCNA", ref = "dev", 
                          force = FALSE, quiet = TRUE,
                          build_vignettes = FALSE)
    install_success <- TRUE
    cat("✅ hdWGCNA安装成功（第", i, "次尝试）\n")
    break
  }, error = function(e) {
    cat("❌ 第", i, "次失败:", conditionMessage(e), "\n")
    if (i < max_try) cat("等待10秒后重试...\n")
    Sys.sleep(10)
  })
}

if (!install_success) {
  cat("\n❌ hdWGCNA安装最终失败\n")
  cat("备选方案1: 使用本地简化模块分析结果（已完成的64基因模块）\n")
  cat("备选方案2: 明日使用GitHub Token认证后重试\n")
  cat("备选方案3: 下载源码本地安装\n")
  
  # 写入失败标记
  writeLines("FAILED", "03_results/09_single_cell/hdwgcna/INSTALL_STATUS.txt")
  quit(status = 0)  # 不退出，继续执行后续备选分析
}

library(hdWGCNA)
cat("✅ hdWGCNA加载成功，版本:", packageVersion("hdWGCNA"), "\n")
writeLines("SUCCESS", "03_results/09_single_cell/hdwgcna/INSTALL_STATUS.txt")

# --- 2. 加载数据 ---
cat("\n========== 2. 加载数据 ==========\n")
srt <- readRDS("03_results/09_single_cell/34_srt_with_module.rds")
cat("✅ 加载对象:", ncol(srt), "cells x", nrow(srt), "genes\n")

# --- 3. hdWGCNA分析 ---
cat("\n========== 3. hdWGCNA模块分析 ==========\n")

# 设置默认assay为RNA
DefaultAssay(srt) <- "RNA"

# 使用线粒体基因子集（避免全基因组计算过慢）
mito_genes <- c(
  grep("^NDUF", rownames(srt), value = TRUE),
  grep("^COX", rownames(srt), value = TRUE),
  grep("^ATP5", rownames(srt), value = TRUE),
  grep("^UQCRC", rownames(srt), value = TRUE),
  grep("^SDH", rownames(srt), value = TRUE),
  grep("^CYC1", rownames(srt), value = TRUE),
  "NDUFB7", "PPARGC1A", "NRF1", "TFAM", "ESRRA"
)
mito_genes <- unique(mito_genes[mito_genes %in% rownames(srt)])
cat("线粒体候选基因:", length(mito_genes), "个\n")

# 构建hdWGCNA对象
tryCatch({
  srt <- SetupForWGCNA(
    srt,
    gene_select = "custom",
    gene_list = mito_genes,
    wgcna_name = "mito"
  )
  
  # 构建metacells
  srt <- MetacellsByGroups(
    srt,
    group.by = c("condition", "seurat_clusters"),
    k = 10,
    max_shared = 5,
    ident.group = "condition"
  )
  cat("✅ Metacells构建完成\n")
  
  # 设置表达矩阵
  srt <- SetDatExpr(
    srt,
    group_name = "HF",
    group.by = "condition",
    assay = "RNA",
    slot = "data"
  )
  
  # 测试软阈值
  srt <- TestSoftPowers(srt, networkType = "signed")
  cat("✅ 软阈值测试完成\n")
  
  # 构建网络（使用简化参数加速）
  srt <- ConstructNetwork(
    srt,
    soft_power = 6,
    setDatExpr = FALSE,
    tom_outdir = "03_results/09_single_cell/hdwgcna/TOM",
    overwrite_tom = TRUE
  )
  cat("✅ 网络构建完成\n")
  
  # 模块特征基因
  srt <- ModuleEigengenes(srt)
  cat("✅ 模块特征基因计算完成\n")
  
  # 模块连通性
  srt <- ModuleConnectivity(srt)
  cat("✅ 模块连通性计算完成\n")
  
  # 获取NDUFB7所在模块
  mods <- GetModules(srt)
  ndufb7_module <- mods$module[mods$gene_name == "NDUFB7"]
  cat("\n🎯 NDUFB7所属模块:", ndufb7_module, "\n")
  
  # 同模块基因
  if (length(ndufb7_module) > 0 && !is.na(ndufb7_module)) {
    same_mod <- mods$gene_name[mods$module == ndufb7_module & mods$gene_name != "NDUFB7"]
    cat("同模块基因数:", length(same_mod), "\n")
    cat("Top 10 hub基因:", head(same_mod[order(mods$kWithin[mods$module == ndufb7_module & mods$gene_name != "NDUFB7"], decreasing = TRUE)], 10), "\n")
    
    # 检查PGC-1α/NRF1/TFAM
    check_genes <- c("PPARGC1A", "NRF1", "TFAM", "ESRRA")
    for (g in check_genes) {
      if (g %in% mods$gene_name) {
        g_mod <- mods$module[mods$gene_name == g]
        cat("  ", g, "模块:", g_mod, ifelse(g_mod == ndufb7_module, "✅同模块", "❌不同模块"), "\n")
      }
    }
    
    # 保存模块基因
    mod_df <- mods[mods$module == ndufb7_module, c("gene_name", "module", "kWithin")]
    write.csv(mod_df, "03_results/09_single_cell/hdwgcna/24_hdwgcna_ndufb7_module.csv", row.names = FALSE)
  }
  
  # 保存结果
  saveRDS(srt, "03_results/09_single_cell/hdwgcna/25_srt_hdwgcna.rds")
  cat("\n✅ hdWGCNA分析完成，结果已保存\n")
  
}, error = function(e) {
  cat("\n❌ hdWGCNA分析失败:", conditionMessage(e), "\n")
  cat("建议: 使用简化版模块分析结果作为Figure 4\n")
})

cat("\n🎉 脚本执行完毕！\n")
