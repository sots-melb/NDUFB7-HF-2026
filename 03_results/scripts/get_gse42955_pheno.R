library(GEOquery)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE42955"

cat("下载GSE42955完整GEO对象...\n")
gse <- getGEO("GSE42955", GSEMatrix = FALSE)

cat("提取样本表型...\n")
pd <- pData(phenoData(gse[[1]]))
cat("样本数:", nrow(pd), "\n")
cat("可用列名:", names(pd), "\n")

# 保存完整表型
fwrite(as.data.table(pd), file=file.path(out_dir, "GSE42955_full_phenotype.txt"), sep="\t")

# 查找病因相关列
disease_cols <- names(pd)[grepl("disease|diagnosis|etiology|condition|status", names(pd), ignore.case=TRUE)]
cat("\n病因相关列:", disease_cols, "\n")

if(length(disease_cols) > 0) {
  for(col in disease_cols) {
    cat("\n---", col, "---\n")
    print(table(pd[[col]]))
  }
} else {
  cat("\n未找到病因列，显示所有列样本值:\n")
  print(head(pd[, 1:min(5, ncol(pd))]))
}

cat("\n完成\n")
