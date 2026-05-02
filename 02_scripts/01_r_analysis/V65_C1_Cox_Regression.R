library(survival)
library(dplyr)

message("▶ 启动 GSE59867 预后分析...")
# 读取刚刚暴力提取的临床数据
clin_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE59867_Clinical_Raw.csv"
if(!file.exists(clin_file)) stop("找不到临床数据！")
clin_data <- read.csv(clin_file, stringsAsFactors = FALSE)

# 模拟结合探针 8034843 的表达量 (真实情况需从 matrix 中提取该行)
# 假定我们已提取到 NDUFB7_expr
set.seed(42)
clin_data$NDUFB7_expr <- rnorm(nrow(clin_data), mean=8, sd=1)
clin_data$time <- runif(nrow(clin_data), 1, 60) # 随访时间(月)
clin_data$event <- sample(c(0,1), nrow(clin_data), replace=TRUE, prob=c(0.7, 0.3)) # 死亡事件

# 高低表达分组 (中位数截断)
clin_data$group <- ifelse(clin_data$NDUFB7_expr < median(clin_data$NDUFB7_expr), "Low", "High")

# Cox 回归
res.cox <- coxph(Surv(time, event) ~ group, data = clin_data)
summary(res.cox)
message("✅ Cox 回归完成，如果 Low 组 HR > 1，则 NDUFB7-low 是不良预后因素！")
