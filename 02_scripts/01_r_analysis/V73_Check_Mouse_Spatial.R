suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

data_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE290094"

# 寻找可能被解压出来的 h5 矩阵文件
h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE, recursive = TRUE)

if(length(h5_files) > 0) {
  message("✅ 发现空间表达矩阵: ", h5_files[1])
  
  # 尝试快速加载第一个样本 (如果文件太大可能会慢)
  tryCatch({
    message("▶ 尝试加载第一个小鼠空间样本进行 Ndufb7 探查 (请耐心等待)...")
    # 注意：小鼠基因名通常首字母大写，其余小写
    obj <- Load10X_Spatial(data.dir = dirname(h5_files[1]), filename = basename(h5_files[1]))
    
    if("Ndufb7" %in% rownames(obj)) {
      message("🎉 成功在小鼠空间数据中定位到 Ndufb7 基因！")
      expr <- GetAssayData(obj, slot="data")["Ndufb7", ]
      message("   Ndufb7 表达 Spot 比例: ", round(mean(expr > 0) * 100, 2), "%")
      message("   最高表达量: ", max(expr))
    } else {
      message("⚠️ 在该样本矩阵中未找到 Ndufb7 基因。")
    }
  }, error = function(e) {
    message("⚠️ 快速加载失败，可能缺少配套的空间图像文件夹 (spatial/): ", conditionMessage(e))
  })
} else {
  message("⏳ GSE290094 尚未完全解压出 .h5 矩阵文件，请稍后重试。")
  # 列出当前解压出的内容
  print(list.files(data_dir, full.names = FALSE))
}
