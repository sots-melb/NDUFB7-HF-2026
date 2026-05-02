suppressMessages(library(hdf5r))
suppressMessages(library(dplyr))

base <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE214611_FULL_EXTRACT"
samples <- list.dirs(base, full.names=FALSE, recursive=FALSE)
samples <- samples[samples!=""]
message("▶ 发现样本: ", paste(samples, collapse=", "))

results <- data.frame(Sample=character(), TimePoint=character(), N_Cells=integer(), NDUFB7_Mean=numeric(), NDUFB7_Pct=numeric(), ComplexI_Mean=numeric(), ComplexI_Cor=numeric(), FTH1_Cor=numeric(), stringsAsFactors=FALSE)

complexI <- c("NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11","NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NDUFV3")

for(s in samples) {
  h5_path <- file.path(base, s, "filtered_feature_bc_matrix.h5")
  if(!file.exists(h5_path)) next
  message("▶ 处理: ", s)
  tryCatch({
    file <- H5File$new(h5_path, mode="r")
    mat <- file[["matrix"]]
    barcodes <- mat[["barcodes"]][]
    features <- mat[["features"]]
    gene_names <- features[["name"]][]
    data <- mat[["data"]][]
    indices <- mat[["indices"]][]
    indptr <- mat[["indptr"]][]
    shape <- mat[["shape"]][]
    file$close_all()
    
    sparse_mat <- sparseMatrix(i=indices+1, p=indptr, x=data, dims=shape)
    rownames(sparse_mat) <- gene_names
    ndufb7 <- as.numeric(sparse_mat["NDUFB7", ])
    n_cells <- length(ndufb7)
    ndufb7_mean <- mean(ndufb7)
    ndufb7_pct <- mean(ndufb7>0)
    
    # Complex I
    ci_in <- complexI[complexI %in% rownames(sparse_mat)]
    ci_cor <- NA
    if(length(ci_in)>0 && mean(ndufb7)>0) {
      ci_mat <- as.matrix(sparse_mat[ci_in, ])
      ci_mean <- colMeans(ci_mat)
      ci_cor <- cor(ci_mean, ndufb7, use="complete.obs")
    }
    
    # FTH1
    fth1_cor <- NA
    if("FTH1" %in% rownames(sparse_mat)) {
      fth1 <- as.numeric(sparse_mat["FTH1", ])
      fth1_cor <- cor(fth1, ndufb7, use="complete.obs")
    }
    
    # 推断时间点
    time <- "Unknown"
    if(grepl("sham", s, ignore.case=TRUE)) time <- "Sham"
    else if(grepl("1hr", s, ignore.case=TRUE)) time <- "1hr"
    else if(grepl("4hr", s, ignore.case=TRUE)) time <- "4hr"
    else if(grepl("d1", s, ignore.case=TRUE)) time <- "Day1"
    else if(grepl("d3", s, ignore.case=TRUE)) time <- "Day3"
    else if(grepl("d7", s, ignore.case=TRUE)) time <- "Day7"
    else if(grepl("Human", s, ignore.case=TRUE)) time <- "Human_STEMI"
    
    results <- rbind(results, data.frame(Sample=s, TimePoint=time, N_Cells=n_cells, NDUFB7_Mean=ndufb7_mean, NDUFB7_Pct=ndufb7_pct, ComplexI_Mean=ifelse(length(ci_in)>0, mean(colMeans(as.matrix(sparse_mat[ci_in, ]))), NA), ComplexI_Cor=ci_cor, FTH1_Cor=fth1_cor, stringsAsFactors=FALSE))
  }, error=function(e) message("  ❌ ", s, ": ", conditionMessage(e)))
}

results <- results[order(results$TimePoint), ]
message("▶ 时间序列汇总:")
print(results)
out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE214611_TimeSeries.csv"
write.csv(results, out, row.names=FALSE)
message("✅ 保存: ", out)
