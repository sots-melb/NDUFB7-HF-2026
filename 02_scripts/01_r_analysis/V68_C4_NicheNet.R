suppressMessages(library(nichenetr))
suppressMessages(library(dplyr))

message("▶ NicheNet版本: ", packageVersion("nichenetr"))

# 加载表达数据
obj <- readRDS("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Scored.RDS")
message("▶ 对象: ", ncol(obj), " 细胞")

# 获取NDUFB7高表达和低表达细胞
ndufb7 <- obj$NDUFB7_expr
high_cells <- colnames(obj)[ndufb7 >= quantile(ndufb7, 0.75, na.rm=TRUE)]
low_cells <- colnames(obj)[ndufb7 <= quantile(ndufb7, 0.25, na.rm=TRUE)]
message("▶ NDUFB7-high: ", length(high_cells), " | NDUFB7-low: ", length(low_cells))

# 提取表达矩阵
expr <- GetAssayData(obj, assay="RNA", slot="data")
sender_expr <- rowMeans(expr[, high_cells])
receiver_expr <- rowMeans(expr[, low_cells])

# 使用NicheNet内置配体-受体数据库
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

# 定义配体和受体
expressed_genes_sender <- names(sender_expr)[sender_expr > 0.1]
expressed_genes_receiver <- names(receiver_expr)[receiver_expr > 0.1]

# 运行NicheNet
niche_output <- nichenetr::predict_ligand_activities(
  geneset = expressed_genes_sender[1:500],  # Top 500高表达基因作为sender信号
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network
)

message("▶ Top 10 活性配体:")
print(head(niche_output, 10))

out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_NicheNet_Ligands.csv"
write.csv(niche_output, out, row.names=FALSE)
message("✅ 保存: ", out)
