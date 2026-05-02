#!/usr/bin/env Rscript
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)

message("========================================")
message("V124: 2026-05-01 夜间审计")
message("========================================")

# 检查核心产出
checks <- c(
  "03_results/V120_Fig3C/V120_Fig3C_Monocle3_stepwise.png",
  "03_results/V123_Fig4/V123_Fig4_mechanism_integration.png",
  "03_results/V113A_CellChat/V113A_heatmap_fixed.png",
  "03_results/V122_GSE121893/V122_GSE121893_bimodal_density.png",
  "03_results/V121_CellChat_Interpretation/V121_cluster_ndufb7_levels.csv",
  "03_results/V113B_Monocle3_Stepwise/V113B_stepwise_stats.csv"
)

message("\n=== 核心产出检查 ===")
for (f in checks) {
  if (file.exists(f)) {
    sz <- round(file.size(f)/1024, 1)
    message("[✅] ", basename(f), " (", sz, " KB)")
  } else {
    message("[❌] MISSING: ", f)
  }
}

# 叙事状态
message("\n=== 叙事状态 ===")
message("[PIVOTED] 阈值假说 → 阶梯式下降")
message("[FROZEN] Monocle3 breakpoint: pseudotime=6.14, p=4.3e-04")
message("[FROZEN] OXPHOS collapse: Complex I/III/IV下调, p<1e-6")
message("[FROZEN] Ferroptosis index: ACSL4/GPX4 ratio 25.96 vs 24.15")
message("[FROZEN] CellChat: 3742 pairs, ANXA1-FPR1 dominant")

# 明日优先级
message("\n=== 明日优先级 (2026-05-02) ===")
message("P0: Moran's I重建 — 需要GEO原始Visium数据或Space Ranger输出")
message("P0: GitHub SSH配置 — 代码公开是Reviewer硬性要求")
message("P1: Results撰写 — Paragraph 3-4现在可写（素材已齐备）")
message("P1: Methods补全 — Monocle3/CellChat/ACSL4-GPX4方法学段落")
message("P2: GSE57338探针映射 — 解锁Bulk层面全部基因提取")
message("P2: Fig 2优化 — 依赖Moran's I真实值")

message("\n[DONE] V124 Night Audit")
message("[RECOMMENDATION] 今晚止损，明日P0-1优先")
