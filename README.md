# NDUFB7-HF-2026

**NDUFB7 as a mitochondrial energy crisis gatekeeper in heart failure: A pan-death suppressor independent of reactive oxygen species**

## Project Structure
NDUFB7_HF_2026_04_20/
├── 01_data/           # Raw and processed data
│   ├── 01_raw_geo/    # GEO downloads
│   └── 02_gpl_annotation/  # Platform annotations
├── 02_code/           # Analysis scripts (V001-V185)
├── 03_results/        # Figure panels and statistical outputs
│   ├── V178_Fixed/    # GSE57338/GSE59867 gene-level data
│   ├── V181_Comprehensive/  # Mechanism analysis
│   ├── V183_Fixed/    # Fig3 panels
│   ├── V184_Fixed/    # DOX validation
│   ├── V185_Fig4/     # Fig4 panels + Results P4
│   └── ...            # Intermediate analyses
├── 04_manuscript/     # Draft manuscript and supplementary
│   ├── V186_Draft/    # Full manuscript markdown
│   ├── V187_Supplementary/  # Supplementary tables
│   └── V188_Methods/  # Detailed methods
└── README.md
plain
复制

## Key Findings
1. NDUFB7 is consistently downregulated in HF across 5 independent cohorts
2. Genetic evidence supports causal NDUFB7→HF relationship (SMR)
3. NDUFB7 regulates pan-death via ROS-independent energy crisis mechanism
4. DOX directly suppresses NDUFB7; Fer-1 cannot rescue (upstream positioning)
5. Visium spatial transcriptomics reveals infarct-core NDUFB7 nadir

## Reproducibility
All analyses performed in R v4.3.1 and Python 3.10 on Ubuntu 22.04. 
See `04_manuscript/V188_Methods/Methods_Detail.md` for complete parameters.

## Data Availability
- GEO: GSE57338, GSE59867, GSE168742, GSE157282, GSE243655
- GTEx Portal: v8 Heart-LV eQTL
- eQTLGen Consortium: whole-blood eQTL meta-analysis
- HERMES GWAS: heart failure summary statistics

## Citation
[To be updated upon publication]
