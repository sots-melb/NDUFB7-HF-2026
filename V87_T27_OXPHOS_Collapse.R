#!/usr/bin/env Rscript
# V87_T27: OXPHOS多复合体协同崩溃验证

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(reshape2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

complex_genes <- list(
  Complex_I   = c("NDUFB7", "NDUFB8", "NDUFB10", "NDUFA9", "NDUFS1"),
  Complex_II  = c("SDHA", "SDHB", "SDHC", "SDHD"),
  Complex_III = c("UQCRC1", "UQCRC2", "CYTB", "CYC1"),
  Complex_IV  = c("COX4I1", "COX5A", "COX5B", "COX6C"),
  Complex_V   = c("ATP5F1A", "ATP5F1B", "ATP5MC1", "ATP5ME")
)

# 加载数据（优先GSE183852）
files <- list.files("Downloads", pattern = "GSE183852.*\\.(rds|Robj)", full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) files <- list.files("01_data", pattern = "GSE183852", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if (length(files) == 0) stop("未找到GSE183852")

srt <- readRDS(files[1])
if ("cell_type" %in% colnames(srt@meta.data)) {
  cm <- subset(srt, cell_type %in% c("Cardiomyocyte", "CM"))
} else {
  cm <- srt
}
message("=== T27: OXPHOS多复合体崩溃验证 ===")
message("CM子集: ", ncol(cm), " cells")

# 提取表达
all_genes <- unlist(complex_genes)
avail <- intersect(all_genes, rownames(cm))
message("可用基因: ", length(avail), "/", length(all_genes))

if (length(avail) < 5) stop("可用基因过少")

expr <- FetchData(cm, vars = avail)

# 获取分组
if ("condition" %in% colnames(cm@meta.data)) {
  expr$condition <- cm$condition
} else if ("Condition" %in% colnames(cm@meta.data)) {
  expr$condition <- cm$Condition
} else {
  expr$condition <- "Unknown"
  message("[WARN] 无condition列，仅输出描述统计")
}

# 逐复合体统计
results <- data.frame()
for (cname in names(complex_genes)) {
  g <- intersect(complex_genes[[cname]], colnames(expr))
  if (length(g) == 0) next
  expr[[cname]] <- rowMeans(expr[, g, drop = FALSE])
  
  if (length(unique(expr$condition)) >= 2 && expr$condition[1] != "Unknown") {
    grp <- split(expr[[cname]], expr$condition)
    if (length(grp) >= 2) {
      g1 <- grp[[1]]; g2 <- grp[[2]]
      tt <- t.test(g1, g2)
      results <- rbind(results, data.frame(
        Complex = cname,
        N_Genes = length(g),
        Group1_Mean = mean(g1, na.rm = TRUE),
        Group2_Mean = mean(g2, na.rm = TRUE),
        Log2FC = log2(mean(g2, na.rm = TRUE) / mean(g1, na.rm = TRUE)),
        P_Value = tt$p.value,
        Direction = ifelse(mean(g2, na.rm = TRUE) < mean(g1, na.rm = TRUE), "DOWN", "UP")
      ))
    }
  }
}

outdir <- file.path(PROJECT_DIR, "03_results/T27_OXPHOS_Collapse")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (nrow(results) > 0) {
  write.csv(results, file.path(outdir, "T27_stats.csv"), row.names = FALSE)
  message("")
  message("=== T27 结果 ===")
  print(results[, c("Complex", "Log2FC", "P_Value", "Direction")])
  
  down_n <- sum(results$Direction == "DOWN" & results$P_Value < 0.05)
  message("")
  message("显著下调复合体: ", down_n, "/", nrow(results))
  
  if (down_n >= 3) {
    message("[PASS] ≥3个复合体显著下调，支持'多复合体协同崩溃'假说")
    message("[IMPACT] 可作为新增Fig 4A，强化系统级叙事")
  } else if (down_n >= 2) {
    message("[PARTIAL] 2个复合体下调，可提及但弱化'全面崩溃'表述")
  } else {
    message("[WARN] 下调不足，'多复合体崩溃'叙事需降级为'Complex I特异性'")
  }
} else {
  message("[WARN] 无组间比较结果，仅保存描述统计")
}

message("[DONE] 结果保存: ", outdir)
