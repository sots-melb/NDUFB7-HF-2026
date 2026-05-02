=== V60 文件结构最终验证报告 ===
时间: Fri Apr 24 09:30:49 PM CST 2026

## 一、一级目录完整性检查

| 目录 | 状态 | 大小 | 备注 |
|------|------|------|------|
| 01_data/ | ✅ | 93G | 标准目录 |
| 02_scripts/ | ✅ | 321K | 标准目录 |
| 03_results/ | ✅ | 1.2G | 标准目录 |
| 04_paper/ | ✅ | 3.0K | 标准目录 |
| 05_logs/ | ❌ | — | 缺失！ |
| 08_admin/ | ✅ | 47M | 标准目录 |
| 09_tmp/ | ✅ | 353K | 标准目录 |

## 二、关键文件存在性检查

| 文件 | 状态 | 大小 | 用途 |
|------|------|------|------|
| README.md | ❌ | — | 项目说明 |
| 00_PROJECT_STRUCTURE.txt | ✅ | 3.5K | 结构图 |
| 01_data/00_data_governance/00_data_sources_v60.md | ✅ | 6.0K | 数据溯源 |
| 01_data/00_data_governance/00_download_log_v60.csv | ✅ | 2.0K | 下载日志 |
| 04_paper/03_supplementary/DataAvailability_NatureCommunications.txt | ✅ | 1.5K | NC格式DAS |
| 03_results/01_figures/SuppFig3_ForestPlot_4platforms.pdf | ✅ | 30K | 森林图PDF |
| 03_results/01_figures/SuppFig3_ForestPlot_4platforms.tiff | ✅ | 226K | 森林图TIFF |

## 三、Tier S核心资产验证

| 资产 | 路径 | 大小 | 状态 |
|------|------|------|------|
| S01 | 01_data/04_spatial_geo/Kuppe_Nature_2022 | 999M | ✅ |
| S02 | 01_data/01_raw_geo/GSE57338 | 4.2G | ✅ |
| S03 | 01_data/01_raw_geo/GSE116250 | 20M | ✅ |
| S04 | 01_data/01_raw_geo/GSE55296 | 1.8M | ✅ |
| S05 | 01_data/05_proteomics_pride/PXD010154 | 3.1G | ✅ |
| S06 | 01_data/03_mr_data/eQTLGen | 309M | ✅ |
| S07 | 01_data/03_mr_data/gtex | 3.2G | ✅ |

## 四、幽灵目录状态

幽灵目录仍存在: /home/y411869/Projects/NDUFB7_HF_{2026_04_20} (2.7G)
备份位置: /home/y411869/Projects/.backups/ghost_merge_20260424_2114
**建议**: 确认合并无误后执行: `rm -rf "/home/y411869/Projects/NDUFB7_HF_{2026_04_20}"`

## 五、磁盘使用

Filesystem                                                                                                                                      Size  Used Avail Use% Mounted on
10.10.3.60:6789,10.10.3.61:6789,10.10.3.62:6789:/volumes/csi/csi-vol-2510602b-dd51-4751-8640-d5b637adc695/60b5114b-d88f-4961-9e02-4b7e6596b279  1.5T  203G  1.3T  14% /home

## 六、结构优化建议

1. **清理幽灵目录**: 确认备份后删除 `NDUFB7_HF_{2026_04_20}`
2. **清理Downloads**: 使用生成的清理计划脚本释放空间
3. **补充GDS4772**: 若需完整四平台验证，下载GDS4772
4. **解压GSE183852**: .Robj.gz需要解压或直接用Seurat读取
5. **GitHub上传**: 将09_tmp/git_archive_20260424/推送至GitHub获取DOI

=== 验证完成 ===
报告: /home/y411869/Projects/NDUFB7_HF_2026_04_20/03_results/05_audit_reports/v60_structure_final_verify_20260424_2130.md
