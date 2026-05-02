suppressMessages(library(Seurat))

data_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE290094"
h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE, recursive = TRUE)

if(length(h5_files) > 0) {
  message("\n▶ 尝试加载第一个小鼠空间样本...")
  tryCatch({
    obj <- Load10X_Spatial(data.dir = dirname(h5_files[1]), filename = basename(h5_files[1]))
    if("Ndufb7" %in% rownames(obj)) {
      expr <- GetAssayData(obj, slot="data")["Ndufb7", ]
      message("🎉 成功跨物种定位！小鼠 Ndufb7 表达 Spot 比例: ", round(mean(expr > 0) * 100, 2), "%")
    } else {
      message("⚠️ 未找到 Ndufb7。")
    }
  }, error = function(e) message("⚠️ 空间图像数据不完整: ", conditionMessage(e)))
}
