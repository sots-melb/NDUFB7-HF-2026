suppressMessages(library(hdf5r))
h5 <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE214611_extracted/GSM6613090_V_Human_STEMI/filtered_feature_bc_matrix.h5"
file <- H5File$new(h5, mode="r")
message("▶ HDF5结构:")
print(file$ls(recursive=FALSE))
# 10x标准路径
mat <- file[["matrix"]]
message("▶ matrix组结构:")
print(mat$ls(recursive=FALSE))
# 读取
barcodes <- mat[["barcodes"]][]; message("▶ Barcodes: ", length(barcodes))
features <- mat[["features"]]; message("▶ Features组:")
print(features$ls(recursive=FALSE))
gene_names <- features[["name"]][]; message("▶ Genes: ", length(gene_names))
data <- mat[["data"]][]; message("▶ Data长度: ", length(data))
indices <- mat[["indices"]][]; message("▶ Indices长度: ", length(indices))
indptr <- mat[["indptr"]][]; message("▶ Indptr长度: ", length(indptr))
shape <- mat[["shape"]][]; message("▶ Shape: ", paste(shape, collapse=" x "))
# 查找NDUFB7
ndufb7_idx <- which(gene_names=="NDUFB7")
message("▶ NDUFB7索引: ", paste(ndufb7_idx, collapse=", "))
if(length(ndufb7_idx)>0) {
  # 提取该基因的表达
  gene_idx <- ndufb7_idx[1] - 1  # 0-based for CSR
  expr_vals <- c()
  for(i in 1:(length(indptr)-1)) {
    start <- indptr[i] + 1
    end <- indptr[i+1]
    if(end >= start) {
      cols <- indices[start:end] + 1
      vals <- data[start:end]
      if(gene_idx %in% (cols-1)) {
        expr_vals <- c(expr_vals, vals[cols==(gene_idx+1)])
      }
    }
  }
  message("▶ NDUFB7表达细胞数: ", length(expr_vals))
  message("▶ NDUFB7表达摘要: "); print(summary(expr_vals))
}
file$close_all()
message("✅ HDF5审计完成")
