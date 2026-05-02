# V69 最终资产审计 — 生成论文可用的统计摘要

message("========================================")
message("V69 资产审计 | ", Sys.time())
message("========================================")

assets <- list()

# S01: Visium空间
assets[["S01_Visium_FZ"]] <- list(
  Region="Fibrotic_Zone", N_Spots=19525, FZ_Pct=40.2, IZ_Pct=64.9, P=1.5e-10,
  Interpretation="NDUFB7在FZ特异性丢失", Figure="Fig 2A"
)

# S02: GSE183852单核
assets[["S02_LV_SingleNuclei"]] <- list(
  Dataset="GSE183852", N_Nuclei=220752, 
  CM_Score_DCM=0.558, CM_Score_Donor=0.849, CM_P="<0.0001",
  OXPHOS_DCM=0.160, OXPHOS_Donor=0.283, OXPHOS_P="<0.0001",
  NDUFB7_DCM=0.130, NDUFB7_Donor=0.201, NDUFB7_P="<0.0001",
  Ferro_DCM=0.0036, Ferro_Donor=0.0054, Ferro_P=0.3666,
  CM_Subset_NDUFB7_P=0.541, CM_Subset_OXPHOS_P="7.9e-05",
  Interpretation="CM丢失为主因，OXPHOS障碍，铁死亡不显著",
  Figure="Fig 4A-B"
)

# S03: GSE109816 LA单细胞
assets[["S03_LA_SingleCell"]] <- list(
  Dataset="GSE109816", N_Cells=9994,
  NCM_Median=4, CM_Median=1, P=4.4e-82,
  Interpretation="LA中NCM高表达，CM低表达",
  Figure="Supp"
)

# S04: GSE168742
assets[["S04_GSE168742"]] <- list(
  Control_CM=462, HF_CM=187, P=1.5e-10,
  Interpretation="第三人类数据集验证",
  Figure="Fig 4C"
)

# S05: 小鼠TAC
assets[["S05_Mouse_TAC"]] <- list(
  Dataset="GSE315590", N_Cells=10164,
  W8_Mean=203, W16_Mean=46,
  Interpretation="小鼠8W↑→16W↓双相模式",
  Figure="Fig 3B"
)

# S06: GSE57338铁死亡
assets[["S06_Bulk_Ferroptosis"]] <- list(
  Dataset="GSE57338", N=313,
  R=-0.444, P="1.4e-16", N_Genes=5,
  Interpretation="Bulk水平铁死亡显著，单细胞被稀释",
  Figure="Fig 4D"
)

# S07: PROGENy阴性对照
assets[["S07_PROGENy"]] <- list(
  Method="PROGENy", N_Random=1000, P=0.0000,
  Interpretation="NDUFB7效应特异性确认",
  Figure="Supp"
)

# S08: eQTL方向
assets[["S08_eQTL_Direction"]] <- list(
  eQTLGen="T allele", GTEx="T allele", Direction="一致/相反待确认",
  Figure="Fig 6A"
)

# S09: MR因果
assets[["S09_MR"]] <- list(
  eQTLGen_Beta=0.208, eQTLGen_P=0.117,
  GTEx_Beta=0.081, GTEx_P=0.465,
  Interpretation="MR不显著，SMR待补充",
  Figure="Fig 6A"
)

# S10: GSE59867描述性
assets[["S10_GSE59867_PBMC"]] <- list(
  Day1_N=111, HF=9, NonHF=8, Healthy=94,
  HF_Mean=8.074, NonHF_Mean=7.880, Healthy_Mean=8.052,
  KW_P=0.0301, HF_vs_NonHF_P=0.0745,
  Interpretation="炎症代偿标志物，非独立病因",
  Figure="Fig 6B/Discussion"
)

# S11: GSE154170
assets[["S11_GSE154170"]] <- list(
  N=16, Mean=100.423, Group="All_Control",
  Interpretation="表达基线参考，降级为基线",
  Figure="Supp"
)

# S12: GPL11532
assets[["S12_GPL11532"]] <- list(
  Correct="8034843", Deprecated="7896736",
  Interpretation="探针映射验证",
  Figure="Methods"
)

# S14: Moran's I
assets[["S14_Moran"]] <- list(
  Control=0.015, FZ=0.119, IZ=0.041,
  Interpretation="FZ显著聚集",
  Figure="Fig 2B"
)

# S15: 五平台森林图
assets[["S15_Forest"]] <- list(
  N=431, P_Range="0.44-0.52",
  Interpretation="Affymetrix无差异，count DCM>ICM",
  Figure="Fig 3A"
)

# S16: GSE214611人类STEMI（新增）
assets[["S16_GSE214611_Human"]] <- list(
  N_Cells=1551, NDUFB7_Pct=0.481, NDUFB7_Median=0, NDUFB7_Mean=0.874,
  ComplexI_Cor=0.597, FTH1_Cor=0.487, GPX4_Cor=0.363,
  Interpretation="急性期NDUFB7低表达，与Complex I强协同",
  Figure="Fig 3C"
)

# S17: KO替代（新增）
assets[["S17_KO_Alternative"]] <- list(
  N_Genes=26240, Down_Regulated=28, Up_Regulated=0,
  Top_Genes="MB, MYL3, MT-CO1, TNNC1, S100A1",
  Ferro_FTH1=0.100, Ferro_GPX4=0.051,
  Interpretation="KO后线粒体+收缩基因下调，铁死亡基因弱正相关",
  Figure="Fig 5B"
)

# 打印
for(name in names(assets)) {
  message("\n--- ", name, " ---")
  for(k in names(assets[[name]])) {
    message("  ", k, ": ", assets[[name]][[k]])
  }
}

# 保存为CSV
out_df <- do.call(rbind, lapply(names(assets), function(n) {
  as.data.frame(assets[[n]], stringsAsFactors=FALSE)
}))
write.csv(out_df, "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V69_Final_Asset_Summary.csv", row.names=FALSE)
message("\n✅ 资产汇总已保存")
