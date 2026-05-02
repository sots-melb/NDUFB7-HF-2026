if(!require("nichenetr", quietly=TRUE)) {
  message("▶ 安装 nichenetr...")
  install.packages("nichenetr", repos="https://cloud.r-project.org")
}
if(!require("nichenetr", quietly=TRUE)) {
  message("▶ 尝试GitHub...")
  if(!require("remotes", quietly=TRUE)) install.packages("remotes")
  remotes::install_github("saeyslab/nichenetr")
}
library(nichenetr)
message("✅ nichenetr 版本: ", packageVersion("nichenetr"))
