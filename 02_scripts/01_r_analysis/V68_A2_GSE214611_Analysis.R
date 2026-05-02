suppressMessages(library(hdf5r))
suppressMessages(library(Matrix))

h5 <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE214611_extracted/GSM6613090_V_Human_STEMI/filtered_feature_bc_matrix.h5"
file <- H5File$new(h5, mode="r")
mat <- file[["matrix"]]

barcodes <- mat[["barcodes"]][]
features <- mat[["features"]]
gene_names <- features[["name"]][]
data <- mat[["data"]][]
indices <- mat[["indices"]][]
indptr <- mat[["indptr"]][]
shape <- mat[["shape"]][]
file$close_all()

# 构建稀疏矩阵
sparse_mat <- sparseMatrix(i=indices+1, p=indptr, x=data, dims=shape)
rownames(sparse_mat) <- gene_names
colnames(sparse_mat) <- barcodes

message("▶ 矩阵: ", nrow(sparse_mat), " x ", ncol(sparse_mat))
message("▶ 非零元素: ", length(data), " (密度: ", round(length(data)/prod(shape)*100, 2), "%)")

# NDUFB7表达
ndufb7_expr <- as.numeric(sparse_mat["NDUFB7", ])
df <- data.frame(Barcode=barcodes, NDUFB7=ndufb7_expr, stringsAsFactors=FALSE)
message("▶ NDUFB7表达统计:")
print(summary(df$NDUFB7))
message("▶ NDUFB7>0 比例: ", round(mean(df$NDUFB7>0), 4))

# 分位数分层
df$Level <- ifelse(df$NDUFB7==0, "Zero",
            ifelse(df$NDUFB7<=1, "Low",
            ifelse(df$NDUFB7<=3, "Medium", "High")))
message("▶ NDUFB7分层:")
print(table(df$Level))

# 铁死亡基因共表达（在NDUFB7>0的细胞中）
ferro_genes <- c("FTL","FTH1","SLC7A11","GPX4","SAT1","ACSL4","NFE2L2","KEAP1","LPCAT3","PTGS2","ALOX15","GLS2")
ferro_in <- ferro_genes[ferro_genes %in% rownames(sparse_mat)]
message("▶ 铁死亡基因命中: ", paste(ferro_in, collapse=", "))
if(length(ferro_in)>0) {
  ferro_mat <- as.matrix(sparse_mat[ferro_in, df$NDUFB7>0])
  cor_vec <- apply(ferro_mat, 1, function(g) cor(g, df$NDUFB7[df$NDUFB7>0], use="complete.obs"))
  ferro_cor <- data.frame(Gene=names(cor_vec), Correlation=cor_vec, stringsAsFactors=FALSE)
  ferro_cor <- ferro_cor[order(-abs(ferro_cor$Correlation)), ]
  message("▶ 铁死亡基因与NDUFB7相关（仅在NDUFB7>0细胞）:")
  print(ferro_cor)
}

# 线粒体复合体I基因
complexI <- c("NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11","NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NDUFV3","ND1","ND2","ND3","ND4","ND4L","ND5","ND6")
complexI_in <- complexI[complexI %in% rownames(sparse_mat)]
if(length(complexI_in)>0) {
  ci_mat <- as.matrix(sparse_mat[complexI_in, df$NDUFB7>0])
  ci_expr <- colMeans(ci_mat)
  message("▶ Complex I平均表达 vs NDUFB7相关: r=", round(cor(ci_expr, df$NDUFB7[df$NDUFB7>0]), 3))
}

out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE214611_Human_STEMI_NDUFB7.csv"
write.csv(df, out, row.names=FALSE)
message("✅ 保存: ", out)
