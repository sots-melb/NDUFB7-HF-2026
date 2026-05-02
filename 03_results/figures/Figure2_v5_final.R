library(data.table)
library(ggplot2)
library(gridExtra)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/figures"

cat("生成Figure 2 v5: 四平台最终版（含GSE55296矛盾发现）\n")

# ---------- 平台1: GSE57338 ----------
gse57338_ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
gse57338_pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")
vals_57338 <- as.numeric(unlist(gse57338_ndufb7[, -1]))
ds_57338 <- gse57338_pheno[["disease.status"]]
nf_57338 <- vals_57338[ds_57338=="non-failing"]
dcm_57338 <- vals_57338[ds_57338=="idiopathic dilated CMP"]
isc_57338 <- vals_57338[ds_57338=="ischemic"]

# ---------- 平台2: GDS4772 ----------
gds_strat <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_stratified_fixed.txt")
nf_gds <- gds_strat[Group=="NF"]$Expression
dcm_gds <- gds_strat[Group=="DCM"]$Expression

# ---------- 平台3: GSE116250 ----------
gse116250
cd ~/Projects/NDUFB7_HF_2026_04_20/03_results/scripts

cat > b_family_network_fixed.R << 'REOF'
library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("B亚家族网络构建（修复版：部分匹配）\n")
cat("=" ,rep("=", 69), "\n", sep="")

b_family <- c("NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", 
              "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11")

plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)

# 找基因名列（mrna_assignment或类似）
gene_col <- NULL
for(col in names(plat)) {
  if(any(grepl("gene|symbol|mrna|transcript", col, ignore.case=TRUE))) {
    gene_col <- col
    break
  }
}

cat("基因名列:", gene_col, "\n")
cat("基因名示例:\n")
print(head(plat[[gene_col]], 3))

# 部分匹配：搜索包含"NDUFB"的行
b_hits <- plat[grepl("NDUFB[0-9]+", plat[[gene_col]], ignore.case=TRUE)]
cat("\nB家族探针匹配:", nrow(b_hits), "\n")

if(nrow(b_hits) > 0) {
  # 提取基因符号（从复杂注释中）
  get_gene_symbol <- function(x) {
    # 尝试从注释中提取NDUFB+数字
    m <- regmatches(x, regexpr("NDUFB[0-9]+", x, ignore.case=TRUE))
    if(length(m) > 0) return(toupper(m[1])) else return(NA)
  }
  
  b_hits$GeneSymbol <- sapply(b_hits[[gene_col]], get_gene_symbol)
  b_hits <- b_hits[!is.na(b_hits$GeneSymbol)]
  
  # 去重：同一基因可能多个探针，取第一个
  b_unique <- b_hits[!duplicated(b_hits$GeneSymbol)]
  cat("去重后B家族成员:", nrow(b_unique), "\n")
  print(b_unique[, c(names(plat)[1], gene_col, "GeneSymbol"), with=FALSE])
  
  # 提取GSE57338表达
  gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
  
  probe_ids <- as.character(b_unique[[1]])
  expr_ids <- as.character(gse57338[[1]])
  
  found_mask <- expr_ids %in% probe_ids
  cat("\n在表达矩阵中匹配:", sum(found_mask), "/", length(probe_ids), "\n")
  
  if(sum(found_mask) > 0) {
    b_expr <- gse57338[found_mask, ]
    b_expr$GeneSymbol <- b_unique$GeneSymbol[match(as.character(b_expr[[1]]), probe_ids)]
    
    fwrite(b_expr, "~/Projects/NDUFB7_HF_2026_04_20/03_results/b_family_expression.txt", sep="\t")
    cat("✅ B家族表达已保存\n")
    
    # 计算相关性
    b_mat <- as.matrix(b_expr[, -c(1, ncol(b_expr))])
    rownames(b_mat) <- b_expr$GeneSymbol
    
    b_mat_complete <- b_mat[complete.cases(b_mat), ]
    cat("完整基因数:", nrow(b_mat_complete), "\n")
    
    if(nrow(b_mat_complete) >= 3) {
      cor_mat <- cor(t(b_mat_complete), use="pairwise.complete.obs")
      
      fwrite(as.data.table(cor_mat, keep.rownames="Gene"), 
             "~/Projects/NDUFB7_HF_2026_04_20/03_results/b_family_correlation_matrix.txt", sep="\t")
      
      cat("\nB家族相关性矩阵:\n")
      print(round(cor_mat, 2))
      
      # NDUFB7特异性
      if("NDUFB7" %in% rownames(cor_mat)) {
        cat("\n【NDUFB7与B家族共表达】\n")
        ndufb7_corr <- sort(cor_mat["NDUFB7", ], decreasing=TRUE)
        print(round(ndufb7_corr, 3))
        
        cat("\n最强共表达:", names(ndufb7_corr)[2], "(r=", round(ndufb7_corr[2], 3), ")\n")
        cat("最弱共表达:", names(ndufb7_corr)[length(ndufb7_corr)], "(r=", round(ndufb7_corr[length(ndufb7_corr)], 3), ")\n")
        cat("中位数相关性:", round(median(ndufb7_corr[-1]), 3), "\n")
        
        # 判断：NDUFB7是边缘成员还是核心成员？
        median_r <- median(ndufb7_corr[-1])
        if(median_r < 0.3) {
          cat("\n🔍 结论: NDUFB7中位相关性<0.3，支持'边缘成员'假说\n")
        } else if(median_r < 0.5) {
          cat("\n🔍 结论: NDUFB7中位相关性0.3-0.5，弱-中等模块成员\n")
        } else {
          cat("\n🔍 结论: NDUFB7中位相关性>0.5，核心模块成员\n")
        }
      }
    }
  }
} else {
  cat("⚠️ 仍未找到B家族探针\n")
}

cat("\n完成\n")
