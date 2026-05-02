#!/usr/bin/env Rscript
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V133_GSE57338")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V135: GSE57338表达矩阵结构审计")
message("========================================")

rds_file <- list.files(c("01_data","03_results"), pattern="GSE57338.*gene_level.*\\.rds$|GSE57338.*\\.rds$", full.names=TRUE, recursive=TRUE)[1]
if (is.na(rds_file) || !file.exists(rds_file)) {
  message("[FAIL] GSE57338 RDS not found"); quit(status=1)
}

message("[LOAD] ", basename(rds_file))
obj <- readRDS(rds_file)
message("[INFO] Class: ", paste(class(obj), collapse="/"))

# 深度结构探测
cat("\n=== 对象结构 ===\n")
str(obj, max.level=2)

# 如果是list，检查names
if (is.list(obj) && !is.null(names(obj))) {
  message("\nList names: ", paste(names(obj), collapse=", "))
  
  for (nm in names(obj)) {
    item <- obj[[nm]]
    message("\n--- ", nm, " ---")
    message("  Class: ", paste(class(item), collapse="/"))
    if (is.matrix(item) || is.data.frame(item)) {
      message("  Dim: ", nrow(item), " × ", ncol(item))
      message("  Row names (head): ", paste(head(rownames(item), 3), collapse=", "))
      message("  Col names (head): ", paste(head(colnames(item), 3), collapse=", "))
    }
  }
}

# 尝试提取任何包含数字的矩阵
cat("\n=== 自动提取尝试 ===\n")
for (nm in names(obj)) {
  item <- obj[[nm]]
  if (is.matrix(item) || is.data.frame(item)) {
    if (nrow(item) > 100 && ncol(item) > 10) {
      message("[CANDIDATE] ", nm, ": ", nrow(item), " × ", ncol(item))
      # 尝试找NDUFB7
      if ("NDUFB7" %in% rownames(item)) {
        message("[FOUND] NDUFB7 in ", nm)
        expr <- as.numeric(item["NDUFB7", ])
        message("  n=", length(expr), " mean=", round(mean(expr),2))
        write.csv(data.frame(sample=colnames(item), NDUFB7=expr), 
                  file.path(outdir, "V135_NDUFB7_from_list.csv"), row.names=FALSE)
      }
    }
  }
}

# 如果是ExpressionSet
if (inherits(obj, "ExpressionSet")) {
  message("\n[ExpressionSet detected]")
  message("  exprs dim: ", paste(dim(Biobase::exprs(obj)), collapse=" × "))
  message("  featureNames (head): ", paste(head(Biobase::featureNames(obj)), 3), collapse=", ")
  if ("NDUFB7" %in% Biobase::featureNames(obj)) {
    expr <- as.numeric(Biobase::exprs(obj)["NDUFB7", ])
    message("[FOUND] NDUFB7: n=", length(expr))
  }
}

message("\n[DONE] V135: ", outdir)
message("[ACTION] 根据结构输出，定制下一步提取策略")
