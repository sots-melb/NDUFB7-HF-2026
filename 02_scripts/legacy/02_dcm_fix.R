#!/usr/bin/env Rscript
setwd("~/Projects/NDUFB7_HF_{2026_04_20}")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     GSE183852 DCM数据修复与整合                          ║\n")
cat("║     解决0字节问题，重新下载或修复                        ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# --- 0. 检查文件状态 ---
gz_file <- "01_data/01_raw_geo/GSE183852/GSE183852_DCM_Integrated.Robj.gz"
obj_file <- "01_data/01_raw_geo/GSE183852/GSE183852_DCM_Integrated.Robj"

cat("========== 0. 文件状态检查 ==========\n")
if (file.exists(gz_file)) {
  gz_info <- file.info(gz_file)
  cat("GZ文件存在，大小:", gz_info$size, "bytes (", round(gz_info$size/1024/1024, 2), "MB)\n")
  
  # 检查是否为有效gzip
  is_valid <- system(paste("file", gz_file), intern = TRUE)
  cat("文件类型:", is_valid, "\n")
  
  if (gz_info$size < 1000) {
    cat("⚠️ GZ文件过小（<1KB），可能已损坏，需要重新下载\n")
    need_download <- TRUE
  } else {
    # 尝试读取gzip头
    con <- gzfile(gz_file, "rb")
    header <- tryCatch({
      readBin(con, raw(), n = 10)
    }, error = function(e) NULL)
    close(con)
    
    if (is.null(header) || length(header) == 0) {
      cat("❌ GZ文件无法读取，需要重新下载\n")
      need_download <- TRUE
    } else {
      cat("✅ GZ文件头可读，尝试重新解压...\n")
      need_download <- FALSE
    }
  }
} else {
  cat("❌ GZ文件不存在，需要下载\n")
  need_download <- TRUE
}

# --- 1. 重新下载（如需要） ---
if (need_download) {
  cat("\n========== 1. 重新下载GSE183852 ==========\n")
  
  # GEO FTP路径
  ftp_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz"
  https_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz"
  
  cat("⏳ 尝试FTP下载...\n")
  dl_result <- tryCatch({
    download.file(ftp_url, gz_file, mode = "wb", timeout = 300)
    TRUE
  }, error = function(e) {
    cat("FTP失败:", conditionMessage(e), "\n")
    FALSE
  })
  
  if (!dl_result || !file.exists(gz_file) || file.info(gz_file)$size < 1000) {
    cat("⏳ FTP失败，尝试HTTPS...\n")
    dl_result <- tryCatch({
      download.file(https_url, gz_file, mode = "wb", timeout = 300, method = "libcurl")
      TRUE
    }, error = function(e) {
      cat("HTTPS失败:", conditionMessage(e), "\n")
      FALSE
    })
  }
  
  if (!file.exists(gz_file) || file.info(gz_file)$size < 1000) {
    cat("\n❌ 自动下载失败，请手动下载:\n")
    cat("  1) 浏览器访问: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183852\n")
    cat("  2) 找到 'GSE183852_DCM_Integrated.Robj.gz' (约13.5MB)\n")
    cat("  3) 保存到:", gz_file, "\n")
    cat("  4) 重新运行此脚本\n")
    
    # 写入失败记录
    writeLines("DOWNLOAD_FAILED", "03_results/12_dcm_integration/STATUS.txt")
    quit(status = 0)
  }
  
  cat("✅ 下载完成，大小:", file.info(gz_file)$size, "bytes\n")
}

# --- 2. 解压（使用R.utils作为备选） ---
cat("\n========== 2. 解压文件 ==========\n")

# 删除旧的0字节文件
if (file.exists(obj_file) && file.info(obj_file)$size == 0) {
  cat("删除旧的0字节文件...\n")
  file.remove(obj_file)
}

# 尝试系统gunzip
unzip_success <- FALSE
tryCatch({
  system(paste("gunzip -c", gz_file, ">", obj_file))
  if (file.exists(obj_file) && file.info(obj_file)$size > 1000) {
    cat("✅ 系统gunzip成功，大小:", file.info(obj_file)$size, "\n")
    unzip_success <- TRUE
  }
}, error = function(e) {
  cat("系统gunzip失败:", conditionMessage(e), "\n")
})

# 如果系统gunzip失败，使用R.utils
if (!unzip_success) {
  cat("⏳ 尝试R.utils::gunzip...\n")
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    install.packages("R.utils", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(R.utils)
  
  tryCatch({
    R.utils::gunzip(gz_file, destname = obj_file, overwrite = TRUE, remove = FALSE)
    if (file.exists(obj_file) && file.info(obj_file)$size > 1000) {
      cat("✅ R.utils解压成功，大小:", file.info(obj_file)$size, "\n")
      unzip_success <- TRUE
    }
  }, error = function(e) {
    cat("R.utils解压失败:", conditionMessage(e), "\n")
  })
}

if (!unzip_success) {
  cat("\n❌ 解压失败，文件可能已损坏\n")
  writeLines("UNZIP_FAILED", "03_results/12_dcm_integration/STATUS.txt")
  quit(status = 0)
}

# --- 3. 加载与分析 ---
cat("\n========== 3. 加载DCM对象 ==========\n")
env_before <- ls()
load(obj_file)
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

if (is.null(srt_dcm)) {
  cat("❌ 未找到Seurat对象，可用对象:", paste(new_objs, collapse = ", "), "\n")
  writeLines("LOAD_FAILED", "03_results/12_dcm_integration/STATUS.txt")
  quit(status = 0)
}

# --- 4. NDUFB7分析 ---
cat("\n========== 4. DCM中NDUFB7分析 ==========\n")

if ("NDUFB7" %in% rownames(srt_dcm)) {
  # 基础统计
  ndufb7_exp <- FetchData(srt_dcm, vars = "NDUFB7")$NDUFB7
  cat("NDUFB7表达范围:", round(range(ndufb7_exp), 3), "\n")
  cat("NDUFB7均值:", round(mean(ndufb7_exp), 3), "\n")
  
  # 按condition分组（如果有）
  if ("condition" %in% colnames(srt_dcm@meta.data)) {
    library(dplyr)
    ndufb7_df <- data.frame(
      cell = colnames(srt_dcm),
      NDUFB7 = ndufb7_exp,
      condition = srt_dcm$condition
    )
    sum_df <- ndufb7_df %>% group_by(condition) %>% 
      summarise(mean = mean(NDUFB7), median = median(NDUFB7), n = n())
    cat("\n按Condition分组:\n")
    print(sum_df)
    
    # 检验
    conds <- unique(srt_dcm$condition)
    if (length(conds) >= 2) {
      wt <- wilcox.test(NDUFB7 ~ condition, data = ndufb7_df)
      cat("\nWilcoxon检验: p =", format(wt$p.value, digits = 3), "\n")
    }
  }
  
  # 查找极低表达亚群（类似Cluster 4）
  low_thresh <- quantile(ndufb7_exp, 0.05)  # 底部5%
  low_cells <- colnames(srt_dcm)[ndufb7_exp <= low_thresh]
  cat("\n底部5%低表达细胞:", length(low_cells), "个 (阈值:", round(low_thresh, 3), ")\n")
  
  # 如果有cluster信息
  if ("seurat_clusters" %in% colnames(srt_dcm@meta.data) || "cluster" %in% colnames(srt_dcm@meta.data)) {
    clust_col <- ifelse("seurat_clusters" %in% colnames(srt_dcm@meta.data), "seurat_clusters", "cluster")
    clust_ndufb7 <- data.frame(
      cluster = srt_dcm@meta.data[[clust_col]],
      NDUFB7 = ndufb7_exp
    ) %>% group_by(cluster) %>% 
      summarise(mean = mean(NDUFB7), median = median(NDUFB7), n = n(), 
                pct_low = mean(NDUFB7 <= low_thresh)) %>%
      arrange(mean)
    cat("\n各Cluster NDUFB7统计:\n")
    print(clust_ndufb7)
    
    # 判断是否有类似Cluster 4的极端低表达亚群
    min_mean <- min(clust_ndufb7$mean)
    if (min_mean < 0.5) {
      cat("\n🎯 发现极端低表达Cluster（均值<0.5），类似HF Cluster 4特征！\n")
    } else {
      cat("\n⚠️ 未发现极端低表达Cluster（最低均值:", round(min_mean, 3), "）\n")
    }
  }
  
  # 可视化
  library(ggplot2)
  p <- FeaturePlot(srt_dcm, features = "NDUFB7", cols = c("grey", "red"), 
                   pt.size = 0.5) + ggtitle("NDUFB7 in DCM")
  ggsave("03_results/12_dcm_integration/40_dcm_ndufb7_umap.pdf", p, width = 6, height = 5)
  cat("\n[保存] 40_dcm_ndufb7_umap.pdf\n")
  
  # Violin图
  if ("condition" %in% colnames(srt_dcm@meta.data)) {
    p2 <- VlnPlot(srt_dcm, features = "NDUFB7", group.by = "condition", pt.size = 0) + 
      ggtitle("NDUFB7 by Condition (DCM)")
    ggsave("03_results/12_dcm_integration/41_dcm_ndufb7_violin.pdf", p2, width = 5, height = 4)
    cat("[保存] 41_dcm_ndufb7_violin.pdf\n")
  }
  
} else {
  cat("❌ DCM对象中无NDUFB7基因\n")
}

# --- 5. 保存 ---
cat("\n========== 5. 保存结果 ==========\n")
saveRDS(srt_dcm, "03_results/12_dcm_integration/42_srt_dcm_processed.rds")
writeLines("SUCCESS", "03_results/12_dcm_integration/STATUS.txt")
cat("✅ DCM整合完成\n")
