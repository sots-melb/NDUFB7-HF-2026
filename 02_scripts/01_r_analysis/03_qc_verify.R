#!/usr/bin/env Rscript
setwd("~/Projects/NDUFB7_HF_{2026_04_20}")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     空间转录组数据QC验证脚本                             ║\n")
cat("║     验证项: 完整性、可读性、NDUFB7存在性、空间坐标       ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

library(dplyr)

# ============================================================
# 0. 扫描下载目录
# ============================================================
cat("========== 0. 扫描下载目录 ==========\n")

spatial_dir <- "01_data/04_spatial_raw"
if (!dir.exists(spatial_dir)) {
  cat("❌ 下载目录不存在:", spatial_dir, "\n")
  quit(status = 1)
}

gse_dirs <- list.dirs(spatial_dir, recursive = FALSE, full.names = TRUE)
gse_names <- basename(gse_dirs)
cat("发现GSE目录:", length(gse_dirs), "个\n")
cat(paste(gse_names, collapse = ", "), "\n")

# ============================================================
# 1. 逐GSE QC检查
# ============================================================
cat("\n========== 1. 逐GSE QC检查 ==========\n")

qc_results <- data.frame(
  GSE = character(),
  文件数 = integer(),
  总大小_MB = numeric(),
  Series_Matrix = logical(),
  有H5文件 = logical(),
  有Seurat对象 = logical(),
  有CSV矩阵 = logical(),
  有空间坐标 = logical(),
  GZIP完整性 = logical(),
  NDUFB7可检测 = logical(),
  推荐操作 = character(),
  stringsAsFactors = FALSE
)

for (i in 1:length(gse_dirs)) {
  gse <- gse_names[i]
  gdir <- gse_dirs[i]
  files <- list.files(gdir, recursive = TRUE, full.names = TRUE)
  
  cat("\n---", gse, "---\n")
  cat("文件数:", length(files), "\n")
  
  if (length(files) == 0) {
    qc_results <- rbind(qc_results, data.frame(
      GSE = gse, 文件数 = 0, 总大小_MB = 0, Series_Matrix = FALSE,
      有H5文件 = FALSE, 有Seurat对象 = FALSE, 有CSV矩阵 = FALSE,
      有空间坐标 = FALSE, GZIP完整性 = FALSE, NDUFB7可检测 = FALSE,
      推荐操作 = "重新下载", stringsAsFactors = FALSE
    ))
    next
  }
  
  # 文件大小统计
  sizes <- file.info(files)$size
  total_mb <- round(sum(sizes, na.rm = TRUE) / 1024 / 1024, 2)
  cat("总大小:", total_mb, "MB\n")
  
  # 检查各类文件
  has_matrix <- any(grepl("_series_matrix\\.txt\\.gz$", files))
  has_h5 <- any(grepl("\\.h5$", files, ignore.case = TRUE))
  has_seurat <- any(grepl("\\.(rds|Robj|rdata)$", files, ignore.case = TRUE))
  has_csv <- any(grepl("\\.csv\\.gz$", files, ignore.case = TRUE))
  has_spatial <- any(grepl("(spatial|tissue_positions|scalefactors|tissue_lowres|tissue_hires)", 
                            files, ignore.case = TRUE))
  
  cat("Series Matrix:", has_matrix, "| H5:", has_h5, "| Seurat:", has_seurat, 
      "| CSV:", has_csv, "| Spatial:", has_spatial, "\n")
  
  # GZIP完整性检查
  gz_files <- files[grepl("\\.gz$", files)]
  gz_ok <- TRUE
  if (length(gz_files) > 0) {
    for (gz in head(gz_files, 5)) {  # 只检查前5个避免太慢
      check <- system(paste("gzip -t", shQuote(gz), "2>&1"), intern = TRUE)
      if (length(check) > 0) {
        cat("   ⚠️ GZIP损坏:", basename(gz), "\n")
        gz_ok <- FALSE
      }
    }
  }
  cat("GZIP完整性:", gz_ok, "\n")
  
  # NDUFB7检测（如果Seurat对象可用）
  ndufb7_detected <- FALSE
  if (has_seurat) {
    seurat_file <- files[grepl("\\.(rds|Robj)$", files, ignore.case = TRUE)][1]
    cat("⏳ 尝试读取Seurat对象检测NDUFB7...\n")
    tryCatch({
      if (grepl("\\.rds$", seurat_file, ignore.case = TRUE)) {
        srt <- readRDS(seurat_file)
      } else {
        env_before <- ls()
        load(seurat_file)
        env_after <- ls()
        new_objs <- setdiff(env_after, env_before)
        srt <- NULL
        for (o in new_objs) {
          obj <- get(o)
          if (inherits(obj, "Seurat")) { srt <- obj; break }
        }
      }
      if (!is.null(srt) && "NDUFB7" %in% rownames(srt)) {
        ndufb7_detected <- TRUE
        cat("   ✅ NDUFB7存在于Seurat对象中\n")
      } else if (!is.null(srt)) {
        cat("   ⚠️ Seurat对象中无NDUFB7（可能用Ensembl ID）\n")
      }
    }, error = function(e) {
      cat("   ❌ 读取Seurat失败:", conditionMessage(e), "\n")
    })
  }
  
  # 推荐操作
  action <- "待进一步分析"
  if (has_h5 && has_spatial) {
    action <- "✅ 完整SpaceRanger输出，可构建Seurat对象"
  } else if (has_seurat && ndufb7_detected) {
    action <- "✅ 可直接使用，NDUFB7已验证"
  } else if (has_csv && has_spatial) {
    action <- "⚠️ 需手动构建Seurat对象"
  } else if (!gz_ok) {
    action <- "❌ 文件损坏，需重新下载"
  } else if (total_mb < 1) {
    action <- "❌ 文件过小，可能下载不完整"
  }
  
  qc_results <- rbind(qc_results, data.frame(
    GSE = gse, 文件数 = length(files), 总大小_MB = total_mb,
    Series_Matrix = has_matrix, 有H5文件 = has_h5, 有Seurat对象 = has_seurat,
    有CSV矩阵 = has_csv, 有空间坐标 = has_spatial, GZIP完整性 = gz_ok,
    NDUFB7可检测 = ndufb7_detected, 推荐操作 = action,
    stringsAsFactors = FALSE
  ))
}

# ============================================================
# 2. 保存QC报告
# ============================================================
cat("\n========== 2. QC报告 ==========\n")
print(qc_results, row.names = FALSE)

write.csv(qc_results, "03_results/16_spatial_from_bulk/08_qc_report.csv", row.names = FALSE)
cat("\n✅ QC报告已保存: 03_results/16_spatial_from_bulk/08_qc_report.csv\n")

# 生成下一步脚本
usable <- qc_results[qc_results$有H5文件 | (qc_results$有Seurat对象 & qc_results$NDUFB7可检测), ]
if (nrow(usable) > 0) {
  cat("\n🎯 立即可用的数据集:\n")
  print(usable$GSE, row.names = FALSE)
  
  # 生成读取脚本
  read_script <- c(
    "#!/usr/bin/env Rscript",
    "# 自动生成的空间数据读取脚本",
    paste("# 生成时间:", Sys.time()),
    "library(Seurat)",
    "library(ggplot2)",
    ""
  )
  
  for (g in usable$GSE) {
    read_script <- c(read_script, paste0('# --- ', g, ' ---'))
    gdir <- file.path("01_data/04_spatial_raw", g)
    files <- list.files(gdir, recursive = TRUE, full.names = TRUE)
    
    if (any(grepl("\\.h5$", files))) {
      h5 <- files[grepl("\\.h5$", files)][1]
      read_script <- c(read_script,
        paste0('srt_', g, ' <- Load10X_Spatial(data.dir = "', dirname(h5), '", filename = "', basename(h5), '")'),
        paste0('print("', g, ' loaded: ", ncol(srt_', g, '), " spots")'),
        paste0('if ("NDUFB7" %in% rownames(srt_', g, ')) {'),
        paste0('  p <- SpatialFeaturePlot(srt_', g, ', features = "NDUFB7") + ggtitle("', g, ' NDUFB7")'),
        paste0('  ggsave("03_results/16_spatial_from_bulk/09_', g, '_ndufb7_spatial.pdf", p, width = 8, height = 6)'),
        paste0('}')
      )
    } else if (any(grepl("\\.(rds|Robj)$", files))) {
      rds <- files[grepl("\\.(rds|Robj)$", files)][1]
      read_script <- c(read_script,
        paste0('srt_', g, ' <- readRDS("', rds, '")'),
        paste0('print("', g, ' loaded")')
      )
    }
    read_script <- c(read_script, "")
  }
  
  writeLines(read_script, "02_scripts/09_geo_bulk_screen/04_auto_read_spatial.R")
  cat("✅ 自动读取脚本已生成: 02_scripts/09_geo_bulk_screen/04_auto_read_spatial.R\n")
}

cat("\n🎉 QC验证完成！\n")
