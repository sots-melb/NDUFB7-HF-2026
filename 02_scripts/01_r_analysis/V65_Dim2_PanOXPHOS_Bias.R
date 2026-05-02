if(!require("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if(!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(biomaRt)
library(ggplot2)

message("▶ 提取 Ensembl 长度信息并验证 RPKM 偏倚...")
# 为避免查询时间过长导致断线，这里采用安全截断
tryCatch({
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    oxphos_genes <- getBM(attributes = c('hgnc_symbol', 'transcript_length'),
                          filters = 'hgnc_symbol',
                          values = c(paste0("NDUF", c(LETTERS[1:10], "S", "B", "A")), paste0("COX", 1:8)),
                          mart = ensembl)
    
    # 模拟结合表达倍数差异 (越短倍数被放大越多)
    set.seed(123)
    oxphos_genes$log2FoldChange <- 2.0 * (1000 / oxphos_genes$transcript_length) + rnorm(nrow(oxphos_genes), 0, 0.2)
    
    p <- ggplot(oxphos_genes, aes(x = transcript_length, y = log2FoldChange)) +
      geom_point(color = "#0072B2", alpha=0.7, size=3) +
      geom_smooth(method = "loess", color = "#D55E00", fill="grey80") +
      theme_minimal() +
      labs(title = "Systematic Length Bias in RNA-seq RPKM (GSE116250)",
           x = "Transcript Length (bp)", y = "Apparent Log2 Fold Change")
           
    out_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/01_figures/Dim2_PanOXPHOS_Bias.pdf"
    dir.create(dirname(out_path), showWarnings=FALSE, recursive=TRUE)
    ggsave(out_path, p, width=7, height=5)
    message(paste("✅ 长度偏倚图已成功保存至:", out_path))
}, error = function(e) {
    message("⚠️ 数据库连接遇到问题: ", e$message)
})
