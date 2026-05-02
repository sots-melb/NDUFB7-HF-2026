#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
outdir <- file.path(PROJECT, "03_results/V165A_Causal_Hierarchy")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# V164A原始数据手工整理
fer1_pathway <- data.frame(
  Gene = c("GPX4", "SLC7A11", "FTH1", "ACSL4", "LPCAT3", "PTGS2"),
  DMSO_Mean = c(2263, NA, NA, 2643, 1734, 207),
  Fer1_Mean = c(2815, NA, NA, 3100, 1139, 969),
  Direction = c("UP", "UNCLEAR", "UNCLEAR", "UP", "DOWN", "UNCLEAR"),
  Interpretation = c(
    "Fer-1 upregulates GPX4 defense (+24%), consistent with known mechanism",
    "Data format ambiguous, needs re-extraction",
    "Data format ambiguous, needs re-extraction", 
    "Unexpected UP (Fer-1 should suppress drivers); possible feedback activation",
    "Fer-1 downregulates LPCAT3 (-34%), consistent with reduced lipid peroxidation",
    "One outlier (2856) dominates; needs outlier-robust analysis"
  ),
  stringsAsFactors = FALSE
)

fwrite(fer1_pathway, file.path(outdir, "V165B_Fer1_Pathway_Evidence.csv"))
message("[DONE] Fer-1 pathway evidence saved")
