message("启动 SMR 格式审计与转换...")
# 模拟加载已有 eQTLGen 或 GTEx 数据
# SMR .ma 要求列: SNP, A1(effect), A2(other), freq, b, se, p, N
# 此处为处理骨架，直接基于 V62 审计结果(S07, S08)生成符合格式的文件
output_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/smr_formatted/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 如果你本地有原始全量 eQTLGen txt，在这里执行 dplyr::select 与 rename
# write.table(ma_data, file=paste0(output_dir, "NDUFB7_eQTL.ma"), quote=F, row.names=F, col.names=T)
message("✅ SMR .ma 格式审查与转换脚本准备就绪，待关联原始大文件。")
