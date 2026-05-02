if(!require("meta", quietly = TRUE)) install.packages("meta")
library(meta)
message("▶ 计算平台异质性...")
# 五平台数据输入
data <- data.frame(
  Study = c("GSE57338(Affy)", "GDS4772(Affy)", "GSE55296(Count)", "GSE116250(RPKM)"),
  d = c(0.072, 0.369, 0.218, 0.940),
  se = c(0.11, 0.45, 0.15, 0.35),
  Is_RPKM = c(0, 0, 0, 1)
)
m_gen <- metagen(TE = d, seTE = se, studlab = Study, data = data, sm = "SMD")
m_reg <- metareg(m_gen, ~ Is_RPKM)
print(m_reg)
message("✅ Meta 回归完成，查看 Is_RPKM 的 p-value。")
