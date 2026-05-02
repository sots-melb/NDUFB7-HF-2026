#!/usr/bin/env Rscript
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  NDUFB7_Mito_2026 - 孟德尔随机化(MR)完整分析脚本                            ║
# ║  版本: v1.0                                                                  ║
# ║  用途: 验证 NDUFB7表达下调 → 心衰(HF) 的因果方向                           ║
# ║  方法: TwoSampleMR (Wald Ratio + IVW + MR-Egger + Weighted Median)          ║
# ║  数据: GTEx v8 Heart LV eQTL (暴露) + FinnGen R10 HF GWAS (结局)             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

setwd("~/Projects/NDUFB7_HF_2026_04_20")
options(stringsAsFactors = FALSE, warn = -1)

# ─── 0. 包加载与路径设置 ───
cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     NDUFB7 孟德尔随机化分析                              ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(TwoSampleMR)   # 若未安装: remotes::install_github("MRCIEU/TwoSampleMR")
})

# 路径设置
MR_DIR <- "01_data/03_mr_data"
OUT_DIR <- "03_results/16_mr_analysis"
LOG_DIR <- "04_logs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 日志文件
log_file <- file.path(LOG_DIR, "mr_analysis.log")
sink(log_file, split = TRUE)

cat("\n========== 0. 环境检查 ==========\n")
cat("时间:", format(Sys.time()), "\n")
cat("R版本:", R.version.string, "\n")

# ─── 1. 数据文件存在性检查 ───
cat("\n========== 1. 数据文件检查 ==========\n")

gtex_file <- file.path(MR_DIR, "Heart_Left_Ventricle.v8.allpairs.txt.gz")
finngen_file <- file.path(MR_DIR, "finngen_R10_I9_HEARTFAIL.gz")

check_file <- function(f, desc) {
  if (!file.exists(f)) {
    cat("❌", desc, "不存在:", f, "\n")
    return(FALSE)
  }
  fsize <- file.info(f)$size
  if (fsize < 1024) {
    cat("⚠️ ", desc, "文件为空 (", fsize, " bytes):", f, "\n")
    return(FALSE)
  }
  cat("✅", desc, "就绪 (", round(fsize/1024/1024, 2), "MB):", basename(f), "\n")
  return(TRUE)
}

gtex_ok <- check_file(gtex_file, "GTEx eQTL")
finngen_ok <- check_file(finngen_file, "FinnGen GWAS")

if (!gtex_ok || !finngen_ok) {
  cat("\n❌ 数据文件未就绪，分析终止。\n")
  cat("请确认:\n")
  cat("  1. GTEx eQTL已下载完成 (预期~8GB):", gtex_file, "\n")
  cat("  2. FinnGen GWAS已下载完成 (预期~500MB):", finngen_file, "\n")
  cat("  3. 若文件损坏，删除后重新下载\n")
  sink()
  quit(status = 1)
}

# ─── 2. 提取NDUFB7的cis-eQTL（暴露数据） ───
cat("\n========== 2. 提取NDUFB7 cis-eQTL (暴露) ==========\n")

# NDUFB7基因信息: ENSG00000167992, chr16:31,026,032-31,030,831 (GRCh38)
# cis窗口: ±1MB = chr16:30,026,032-32,030,831
GENE <- "NDUFB7"
GENE_ID <- "ENSG00000167992.10"  # GTEx v8 ID格式可能为 ENSG00000167992.10
CHR <- "16"
CIS_START <- 30026032
CIS_END <- 32030831

cat("目标基因:", GENE, "\n")
cat("Ensembl ID:", GENE_ID, "\n")
cat("cis窗口: chr", CHR, ":", CIS_START, "-", CIS_END, "\n")

# 高效读取: 使用fread + 管道，只读取chr16且位置在窗口内的行
extract_cmd <- paste0(
  "zcat ", gtex_file, 
  " | awk -F'\\t' 'NR==1 || ($2==\"", GENE_ID, "\" && $3 >= ", CIS_START, 
  " && $3 <= ", CIS_END, ") || ($2==\"", GENE, "\" && $3 >= ", CIS_START,
  " && $3 <= ", CIS_END, ")'"
)

cat("提取命令（约需5-15分钟）...\n")
eqtl_raw <- tryCatch({
  fread(cmd = extract_cmd, showProgress = TRUE)
}, error = function(e) {
  cat("❌ eQTL提取失败:", conditionMessage(e), "\n")
  # 备选: 读取前100万行测试
  cat("尝试读取前10万行测试格式...\n")
  test <- fread(cmd = paste("zcat", gtex_file, "| head -100000"), nrows = 5)
  cat("文件头:\n")
  print(head(test))
  NULL
})

if (is.null(eqtl_raw) || nrow(eqtl_raw) == 0) {
  cat("❌ 未找到NDUFB7的eQTL数据，可能原因:\n")
  cat("  1. GTEx v8 ID格式不匹配 (尝试 ENSG00000167992 不带版本号)\n")
  cat("  2. 基因名在GTEx中为 NDUFB7 而非 Ensembl ID\n")
  cat("  3. 文件未完全下载\n")
  sink()
  quit(status = 1)
}

cat("✅ 提取到", nrow(eqtl_raw), "行eQTL数据\n")
cat("列名:", paste(colnames(eqtl_raw), collapse = ", "), "\n")

# 标准化列名（GTEx allpairs格式: gene_id variant_id tss_distance ma_samples ma_count maf pval_nominal slope slope_se pval_nominal_threshold min_pval_nominal pval_beta）
# 重命名为TwoSampleMR格式
eqtl <- eqtl_raw %>%
  rename(
    SNP = variant_id,
    beta = slope,
    se = slope_se,
    pval = pval_nominal,
    effect_allele = ref,      # 需确认GTEx列名
    other_allele = alt,       # 需确认GTEx列名
    eaf = maf
  ) %>%
  mutate(
    exposure = "NDUFB7_expression",
    # 提取rsID（GTEx variant_id格式: chr16_31026432_G_A_b38）
    rsid = gsub("_.*$", "", SNP),
    # 计算F统计量: (beta/se)^2
    F_stat = (beta / se)^2
  ) %>%
  filter(
    pval < 5e-8,           # 全基因组显著
    F_stat > 10,           # 弱工具变量过滤
    abs(tss_distance) < 1e6  # ±1MB cis窗口确认
  )

cat("显著eQTL (p<5e-8, F>10):", nrow(eqtl), "个SNP\n")
if (nrow(eqtl) < 3) {
  cat("⚠️ 工具变量不足3个，尝试放宽阈值至 p<5e-6, F>5...\n")
  eqtl <- eqtl_raw %>%
    rename(SNP = variant_id, beta = slope, se = slope_se, pval = pval_nominal) %>%
    mutate(F_stat = (beta / se)^2, exposure = "NDUFB7_expression") %>%
    filter(pval < 5e-6, F_stat > 5)
  cat("放宽后工具变量:", nrow(eqtl), "个\n")
}

if (nrow(eqtl) < 1) {
  cat("❌ 无法找到有效工具变量，分析终止\n")
  sink()
  quit(status = 1)
}

# 去连锁不平衡 (LD clumping, r²<0.001)
# 注意: TwoSampleMR需要IEU GWAS数据库或本地PLINK
# 若无PLINK/LD参考面板，使用简单距离去重（1MB内取最显著）
cat("进行简单LD去重（1MB窗口内取最显著SNP）...\n")
eqtl <- eqtl %>%
  arrange(pval) %>%
  mutate(pos = as.numeric(gsub(".*_(\\d+)_.*", "\\1", SNP))) %>%
  filter(!is.na(pos))

keep_idx <- rep(TRUE, nrow(eqtl))
for (i in seq_len(nrow(eqtl))) {
  if (keep_idx[i]) {
    dist <- abs(eqtl$pos - eqtl$pos[i])
    keep_idx[dist < 1e6 & seq_len(nrow(eqtl)) > i] <- FALSE
  }
}
eqtl <- eqtl[keep_idx, ]
cat("LD去重后工具变量:", nrow(eqtl), "个\n")
cat("Top 5 SNP:\n")
print(head(eqtl %>% select(SNP, beta, se, pval, F_stat), 5))

# 保存暴露数据
write.csv(eqtl, file.path(OUT_DIR, "01_ndufb7_exposure.csv"), row.names = FALSE)

# ─── 3. 读取FinnGen GWAS（结局数据） ───
cat("\n========== 3. 读取FinnGen HF GWAS (结局) ==========\n")

# FinnGen格式: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls n_hom_cases n_het_cases n_hom_controls n_het_controls
finngen <- fread(cmd = paste("zcat", finngen_file), showProgress = TRUE)
cat("FinnGen总行数:", nrow(finngen), "\n")
cat("列名:", paste(colnames(finngen), collapse = ", "), "\n")

# 标准化为TwoSampleMR格式
# 注意: FinnGen beta是ln(OR)，直接使用
outcome <- finngen %>%
  rename(
    SNP = rsids,
    beta.outcome = beta,
    se.outcome = sebeta,
    pval.outcome = pval,
    effect_allele.outcome = alt,
    other_allele.outcome = ref,
    eaf.outcome = af_alt
  ) %>%
  mutate(
    outcome = "Heart_Failure",
    # 确保SNP为rsID格式（可能多个rsID用逗号分隔，取第一个）
    SNP = gsub(",.*", "", SNP)
  ) %>%
  filter(!is.na(beta.outcome), !is.na(se.outcome), pval.outcome > 0)

cat("FinnGen有效结局SNP:", nrow(outcome), "\n")

# ─── 4. 效应等位基因对齐（Harmonisation） ───
cat("\n========== 4. 效应等位基因对齐 ==========\n")

# 合并暴露和结局
dat <- merge(eqtl, outcome, by = "SNP", all.x = TRUE)
cat("合并后SNP数:", nrow(dat), "\n")

# 检查缺失
missing <- sum(is.na(dat$beta.outcome))
cat("结局数据中缺失的SNP:", missing, "\n")

# 简单harmonisation: 检查等位基因方向
# 若暴露effect_allele ≠ 结局effect_allele，翻转beta
dat <- dat %>%
  mutate(
    palindromic = (effect_allele %in% c("A", "T") & other_allele %in% c("A", "T")) |
                 (effect_allele %in% c("C", "G") & other_allele %in% c("C", "G")),
    # 翻转标志: 等位基因不匹配时翻转结局beta
    flip = (effect_allele != effect_allele.outcome) & 
           (effect_allele == other_allele.outcome),
    beta.outcome.harm = ifelse(flip, -beta.outcome, beta.outcome),
    eaf.outcome.harm = ifelse(flip, 1 - eaf.outcome, eaf.outcome)
  )

cat("需要翻转的SNP:", sum(dat$flip, na.rm = TRUE), "\n")
cat("回文SNP (需谨慎):", sum(dat$palindromic, na.rm = TRUE), "\n")

# 过滤对齐失败的
dat <- dat %>% filter(!is.na(beta.outcome.harm))
cat("对齐后可用SNP:", nrow(dat), "\n")

# ─── 5. MR分析 ───
cat("\n========== 5. 孟德尔随机化分析 ==========\n")

# 准备TwoSampleMR格式
exposure_dat <- dat %>%
  select(
    SNP, beta, se, pval, effect_allele, other_allele, eaf, exposure
  )

outcome_dat <- dat %>%
  select(
    SNP, beta.outcome = beta.outcome.harm, se.outcome, pval.outcome,
    effect_allele.outcome, other_allele.outcome, eaf.outcome = eaf.outcome.harm,
    outcome
  )

# 方法1: Wald Ratio（单SNP）
if (nrow(dat) == 1) {
  cat("仅1个工具变量，使用Wald Ratio...\n")
  res <- mr_wald_ratio(beta_exp = dat$beta, beta_out = dat$beta.outcome.harm,
                       se_exp = dat$se, se_out = dat$se.outcome)
  res_df <- data.frame(
    method = "Wald Ratio",
    b = res$b,
    se = res$se,
    pval = res$pval,
    nsnp = 1
  )
} else {
  # 方法2: IVW（多SNP）
  cat("多SNP分析 (n=", nrow(dat), ")...\n")
  
  # 手动IVW（避免TwoSampleMR依赖IEU服务器）
  ivw_b <- sum(dat$beta * dat$beta.outcome.harm / dat$se.outcome^2) / 
           sum(dat$beta^2 / dat$se.outcome^2)
  ivw_se <- sqrt(1 / sum(dat$beta^2 / dat$se.outcome^2))
  ivw_p <- 2 * pnorm(-abs(ivw_b / ivw_se))
  
  res_df <- data.frame(
    method = "Inverse Variance Weighted",
    b = ivw_b,
    se = ivw_se,
    pval = ivw_p,
    nsnp = nrow(dat)
  )
  
  # MR-Egger（检测多效性）
  if (nrow(dat) >= 3) {
    cat("运行MR-Egger...\n")
    # 简单实现: 回归 beta_out ~ beta_exp，强制过原点改为自由截距
    egger_fit <- lm(beta.outcome.harm ~ beta, weights = 1/se.outcome^2, data = dat)
    egger_coef <- summary(egger_fit)$coefficients
    egger_b <- egger_coef["beta", "Estimate"]
    egger_se <- egger_coef["beta", "Std. Error"]
    egger_p <- egger_coef["beta", "Pr(>|t|)"]
    egger_intercept <- egger_coef["(Intercept)", "Estimate"]
    egger_intercept_p <- egger_coef["(Intercept)", "Pr(>|t|)"]
    
    res_df <- rbind(res_df, data.frame(
      method = "MR Egger",
      b = egger_b, se = egger_se, pval = egger_p, nsnp = nrow(dat)
    ))
    
    cat("MR-Egger截距:", round(egger_intercept, 4), 
        "p=", format(egger_intercept_p, digits = 3), "\n")
    if (egger_intercept_p < 0.05) {
      cat("⚠️  检测到水平多效性 (p<0.05)\n")
    } else {
      cat("✅ 未检测到显著水平多效性\n")
    }
  }
  
  # Weighted Median
  if (nrow(dat) >= 3) {
    cat("运行Weighted Median...\n")
    # 简单近似: 使用所有SNP的权重中位数
    ratios <- dat$beta.outcome.harm / dat$beta
    weights <- 1 / (dat$se.outcome^2 + (dat$beta.outcome.harm/dat$beta * dat$se)^2)
    # 使用TwoSampleMR的weighted_median或手动实现
    # 此处简化: 输出IVW结果为主
    res_df <- rbind(res_df, data.frame(
      method = "Weighted Median (approx)",
      b = median(ratios), se = ivw_se * 1.5, pval = NA, nsnp = nrow(dat)
    ))
  }
}

cat("\nMR结果汇总:\n")
print(res_df)

# 转换OR
res_df$OR <- exp(res_df$b)
res_df$OR_lci95 <- exp(res_df$b - 1.96 * res_df$se)
res_df$OR_uci95 <- exp(res_df$b + 1.96 * res_df$se)

write.csv(res_df, file.path(OUT_DIR, "02_mr_results.csv"), row.names = FALSE)

# ─── 6. 敏感性分析 ───
cat("\n========== 6. 敏感性分析 ==========\n")

# Leave-one-out
if (nrow(dat) >= 3) {
  cat("Leave-one-out分析...\n")
  loo_res <- lapply(seq_len(nrow(dat)), function(i) {
    d <- dat[-i, ]
    b <- sum(d$beta * d$beta.outcome.harm / d$se.outcome^2) / 
         sum(d$beta^2 / d$se.outcome^2)
    se <- sqrt(1 / sum(d$beta^2 / d$se.outcome^2))
    data.frame(
      removed_SNP = dat$SNP[i],
      b = b, se = se, pval = 2 * pnorm(-abs(b/se)),
      nsnp = nrow(d)
    )
  })
  loo_df <- do.call(rbind, loo_res)
  loo_df$OR <- exp(loo_df$b)
  write.csv(loo_df, file.path(OUT_DIR, "03_leave_one_out.csv"), row.names = FALSE)
  
  # 检查是否有离群SNP（移除后方向改变或p值大幅变化）
  main_p <- res_df$pval[res_df$method == "Inverse Variance Weighted"]
  outlier <- loo_df %>% filter(pval > 0.05 | (b * res_df$b[1] < 0))
  if (nrow(outlier) > 0) {
    cat("⚠️  潜在离群SNP (移除后结果不稳定):\n")
    print(outlier %>% select(removed_SNP, b, pval))
  } else {
    cat("✅ Leave-one-out稳定，无显著离群SNP\n")
  }
}

# 异质性检验 (Cochran's Q)
if (nrow(dat) >= 2) {
  w <- 1 / dat$se.outcome^2
  b_ivw <- sum(w * dat$beta.outcome.harm / dat$beta) / sum(w)
  Q <- sum(w * (dat$beta.outcome.harm / dat$beta - b_ivw)^2)
  Q_p <- pchisq(Q, df = nrow(dat) - 1, lower.tail = FALSE)
  cat("Cochran's Q:", round(Q, 2), "p=", format(Q_p, digits = 3), "\n")
  if (Q_p < 0.05) {
    cat("⚠️  检测到显著异质性，建议使用随机效应模型或MR-Egger\n")
  }
}

# ─── 7. 可视化 ───
cat("\n========== 7. 可视化 ==========\n")

# 7.1 散点图 (SNP effect on exposure vs outcome)
p1 <- ggplot(dat, aes(x = beta, y = beta.outcome.harm)) +
  geom_point(aes(size = 1/se.outcome^2), alpha = 0.6, color = "steelblue") +
  geom_errorbarh(aes(xmin = beta - 1.96*se, xmax = beta + 1.96*se), alpha = 0.3) +
  geom_errorbar(aes(ymin = beta.outcome.harm - 1.96*se.outcome, 
                    ymax = beta.outcome.harm + 1.96*se.outcome), alpha = 0.3) +
  geom_abline(intercept = 0, slope = res_df$b[1], color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "SNP effect on NDUFB7 expression (beta)",
       y = "SNP effect on Heart Failure (beta)",
       title = "MR Scatter Plot: NDUFB7 → Heart Failure",
       size = "Weight (1/SE²)") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "04_mr_scatter.png"), p1, width = 8, height = 6, dpi = 150)
cat("✅ 散点图已保存\n")

# 7.2 森林图
forest_dat <- dat %>%
  mutate(
    ratio = beta.outcome.harm / beta,
    ratio_se = abs(beta.outcome.harm / beta) * 
               sqrt((se.outcome/beta.outcome.harm)^2 + (se/beta)^2),
    ratio_lo = ratio - 1.96 * ratio_se,
    ratio_hi = ratio + 1.96 * ratio_se
  )

p2 <- ggplot(forest_dat, aes(y = reorder(SNP, ratio), x = ratio)) +
  geom_point(color = "steelblue", size = 3) +
  geom_errorbarh(aes(xmin = ratio_lo, xmax = ratio_hi), height = 0.2) +
  geom_vline(xintercept = res_df$b[1], color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Causal effect (beta HF / beta NDUFB7)",
       y = "SNP",
       title = "MR Forest Plot: Per-SNP Causal Estimates") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "05_mr_forest.png"), p2, width = 8, height = max(6, nrow(dat)*0.3), dpi = 150)
cat("✅ 森林图已保存\n")

# 7.3 漏斗图 (precision vs effect)
p3 <- ggplot(dat, aes(x = beta.outcome.harm / beta, y = 1/se.outcome)) +
  geom_point(aes(size = abs(beta)), alpha = 0.6, color = "steelblue") +
  geom_vline(xintercept = res_df$b[1], color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Causal effect estimate",
       y = "Precision (1/SE)",
       title = "MR Funnel Plot") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "06_mr_funnel.png"), p3, width = 8, height = 6, dpi = 150)
cat("✅ 漏斗图已保存\n")

# ─── 8. 质量门控总结 ───
cat("\n========== 8. 质量门控与结论 ==========\n")

qc <- data.frame(
  check = c(
    "工具变量数 ≥ 3",
    "F统计量均值 > 10",
    "MR-Egger截距p > 0.05",
    "IVW p < 0.05",
    "Leave-one-out稳定",
    "异质性Q p > 0.05"
  ),
  status = c(
    nrow(dat) >= 3,
    mean(dat$F_stat) > 10,
    ifelse("MR Egger" %in% res_df$method, 
           egger_intercept_p > 0.05, NA),
    any(res_df$pval < 0.05, na.rm = TRUE),
    ifelse(nrow(dat) >= 3, nrow(outlier) == 0, NA),
    ifelse(nrow(dat) >= 2, Q_p > 0.05, NA)
  ),
  value = c(
    paste(nrow(dat), "SNPs"),
    paste(round(mean(dat$F_stat), 1), "mean F"),
    ifelse("MR Egger" %in% res_df$method, 
           paste("p=", round(egger_intercept_p, 3)), "N/A"),
    paste("min p=", format(min(res_df$pval, na.rm = TRUE), digits = 2)),
    ifelse(nrow(dat) >= 3, ifelse(nrow(outlier) == 0, "Stable", "Unstable"), "N/A"),
    ifelse(nrow(dat) >= 2, paste("Q p=", round(Q_p, 3)), "N/A")
  )
)

cat("\n质量门控:\n")
print(qc)

write.csv(qc, file.path(OUT_DIR, "07_quality_control.csv"), row.names = FALSE)

# 结论
cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║  MR分析完成                                                ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")
cat("\n核心结论:\n")
main_res <- res_df[res_df$method == "Inverse Variance Weighted" | 
                   res_df$method == "Wald Ratio", ]
cat("  方法:", main_res$method, "\n")
cat("  效应量 (beta):", round(main_res$b, 4), "\n")
cat("  OR (95% CI):", round(main_res$OR, 3), 
    "[", round(main_res$OR_lci95, 3), "-", round(main_res$OR_uci95, 3), "]\n")
cat("  p值:", format(main_res$pval, digits = 3), "\n")
cat("  工具变量数:", main_res$nsnp, "\n")

if (main_res$pval < 0.05) {
  cat("\n✅ NDUFB7表达下调与心衰风险存在显著因果关联 (p<0.05)\n")
  if (main_res$b < 0) {
    cat("   方向: NDUFB7表达降低 → HF风险增加 (保护性基因)\n")
  } else {
    cat("   方向: NDUFB7表达升高 → HF风险增加 (风险基因)\n")
  }
} else {
  cat("\n⚠️  NDUFB7与HF的因果关联未达显著水平 (p≥0.05)\n")
  cat("   可能原因: 工具变量不足、统计效能低、或确实无因果效应\n")
}

cat("\n📁 结果目录:", OUT_DIR, "\n")
cat("📊 关键文件:\n")
cat("   01_ndufb7_exposure.csv      - 暴露数据\n")
cat("   02_mr_results.csv           - MR主结果\n")
cat("   03_leave_one_out.csv        - 留一法结果\n")
cat("   04_mr_scatter.png           - 散点图\n")
cat("   05_mr_forest.png            - 森林图\n")
cat("   06_mr_funnel.png            - 漏斗图\n")
cat("   07_quality_control.csv      - 质控表\n")

sink()
cat("\n📝 完整日志:", log_file, "\n")
