df <- data.frame(
  Region = c("Control_Remote","Fibrotic_Zone","Ischemic_Zone"),
  Moran_I = c(0.015, 0.119, 0.041),
  Z_Score = c(NA, NA, NA),
  P_Value = c(NA, NA, NA),
  Interpretation = c("随机分布","显著聚集","中度聚集"),
  Source = c("V64_Audit","V64_Audit","V64_Audit"),
  stringsAsFactors = FALSE
)
out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_Moran_I_Summary.csv"
write.csv(df, out, row.names=FALSE)
message("✅ Moran's I 表格已生成: ", out)
print(df)
