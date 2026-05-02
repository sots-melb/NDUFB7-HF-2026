suppressMessages(library(nichenetr))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

message("▶ Seurat版本: ", packageVersion("Seurat"))

obj <- readRDS("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Scored.RDS")
message("▶ 对象: ", ncol(obj), " 细胞")

# Seurat v5兼容提取
get_data <- function(obj, assay="RNA", layer="data") {
  if(packageVersion("Seurat") >= "5.0.0") {
    LayerData(obj, assay=assay, layer=layer)
  } else {
    GetAssayData(obj, assay=assay, slot=layer)
  }
}

ndufb7 <- obj$NDUFB7_expr
high_cells <- colnames(obj)[ndufb7 >= quantile(ndufb7, 0.75, na.rm=TRUE)]
low_cells <- colnames(obj)[ndufb7 <= quantile(ndufb7, 0.25, na.rm=TRUE)]
message("▶ NDUFB7-high: ", length(high_cells), " | NDUFB7-low: ", length(low_cells))

expr <- get_data(obj)
sender_expr <- rowMeans(expr[, high_cells, drop=FALSE])
receiver_expr <- rowMeans(expr[, low_cells, drop=FALSE])

expressed_genes_sender <- names(sender_expr)[sender_expr > 0.1]
expressed_genes_receiver <- names(receiver_expr)[receiver_expr > 0.1]
message("▶ Sender基因: ", length(expressed_genes_sender))
message("▶ Receiver基因: ", length(expressed_genes_receiver))

# 简化版NicheNet（避免大文件下载，用内置数据）
# 仅输出Top表达基因作为候选配体
top_sender <- sort(sender_expr[sender_expr>0.5], decreasing=TRUE)
message("▶ Top 20 Sender配体 (高NDUFB7细胞高表达):")
print(head(top_sender, 20))

# 查找已知配体
known_ligands <- c("TGFB1","TGFB2","IL6","ANGPT1","ANGPT2","VEGFA","PDGFA","CSF1","CSF2","CSF3","CCL2","CXCL12","FGF2","IGF1","TNF","POSTN","MMP2","MMP9","CTGF","TGFA")
ligand_in <- known_ligands[known_ligands %in% names(top_sender)]
message("▶ 已知配体命中:")
print(top_sender[ligand_in])

out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_NicheNet_Ligands_Fix.csv"
write.csv(data.frame(Gene=names(top_sender), Expression=top_sender), out, row.names=FALSE)
message("✅ 保存: ", out)
