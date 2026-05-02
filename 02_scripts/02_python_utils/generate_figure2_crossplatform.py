import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: 跨平台NDUFB7表达比较
ax1 = axes[0, 0]
platforms = ['GSE57338\n(Affy 1.1 ST)\n313 samples', 'GSE116250\n(RNA-seq RPKM)\n64 samples', 
             'GDS4772\n(Affy 1.0 ST)\n17 samples', 'PXD010154\n(Proteomics)\n1 sample']
dcm_means = [7.95, 587.80, 8.50, 3.36e9]  # 最后一项是iBAQ
nf_means = [7.88, 512.76, 8.41, 0]
icm_means = [7.85, 607.81, 0, 0]

x = np.arange(len(platforms))
width = 0.25

# 标准化到各自NF的相对值（用于可视化）
# GSE57338: DCM/NF=1.009, ICM/NF=0.996
# GSE116250: DCM/NF=1.146, ICM/NF=1.185
# GDS4772: DCM/NF=1.011

ax1.bar(x - width, [1.009, 1.146, 1.011, 1], width, label='DCM/Condition', color='#E41A1C', alpha=0.8)
ax1.bar(x + width, [0.996, 1.185, 0, 0], width, label='ICM', color='#377EB8', alpha=0.8)
ax1.axhline(y=1.0, color='black', linestyle='--', linewidth=1)
ax1.set_xticks(x)
ax1.set_xticklabels(platforms, fontsize=9)
ax1.set_ylabel('Relative Expression (normalized to NF)', fontsize=11)
ax1.set_title('A. Cross-Platform NDUFB7 Comparison\n(Normalized to Non-Failing)', fontsize=12, fontweight='bold')
ax1.legend()
ax1.set_ylim(0.9, 1.25)

# Panel B: 统计显著性
ax2 = axes[0, 1]
comparisons = ['GSE57338\nDCM vs NF', 'GSE57338\nICM vs NF', 'GSE116250\nDCM vs NF', 
               'GSE116250\nICM vs NF', 'GDS4772\nDCM vs NF']
p_values = [0.3, 0.3, 0.0081, 0.0348, 0.4421]
colors_p = ['gray', 'gray', '#E41A1C', '#E41A1C', 'gray']
bars = ax2.barh(range(len(comparisons)), [-np.log10(p) for p in p_values], color=colors_p, edgecolor='black')
ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
ax2.set_yticks(range(len(comparisons)))
ax2.set_yticklabels(comparisons, fontsize=9)
ax2.set_xlabel('-log10(p-value)', fontsize=11)
ax2.set_title('B. Statistical Significance Across Platforms', fontsize=12, fontweight='bold')
ax2.legend()

# Panel C: RPKM偏倚示意图
ax3 = axes[1, 0]
genes = ['NDUFB7\n(137aa)', 'MT-CO1\n(513aa)', 'NDUFB8\n(156aa)', 'NDUFS1\n(727aa)']
lengths = [137, 513, 156, 727]
# 假设相同reads数，RPKM与长度成反比
rpkm_relative = [max(lengths)/l for l in lengths]
ax3.bar(range(len(genes)), rpkm_relative, color=['#E41A1C', '#4DAF4A', '#FF7F00', '#984EA3'], edgecolor='black')
ax3.set_xticks(range(len(genes)))
ax3.set_xticklabels(genes, fontsize=9)
ax3.set_ylabel('RPKM Relative Bias\n(same reads, shorter = higher RPKM)', fontsize=10)
ax3.set_title('C. RPKM Length Normalization Bias\n(Shorter genes artificially inflated)', fontsize=12, fontweight='bold')

# Panel D: 证据整合总结
ax4 = axes[1, 1]
ax4.axis('off')
ax4.set_xlim(0, 10)
ax4.set_ylim(0, 10)

ax4.text(5, 9, 'Evidence Integration (V48)', ha='center', fontsize=14, fontweight='bold')

ax4.text(5, 7.5, '✅ Consistent:\n• Protein detectable (HPA + PRIDE)\n• Spatial heterogeneity high (Visium)\n• Bulk stable (Affymetrix, 313 samples)', 
         ha='center', fontsize=10, color='darkgreen')

ax4.text(5, 5.0, '⚠️ Inconsistent:\n• RNA-seq RPKM: DCM/ICM > NF (significant)\n• Affymetrix: DCM ≈ NF ≈ ICM (not significant)\n• Likely RPKM bias for 137-aa NDUFB7', 
         ha='center', fontsize=10, color='darkorange')

ax4.text(5, 2.5, '🎯 Conclusion:\n• Spatial pattern (FZ loss) is robust\n• Bulk direction remains uncertain\n• Protein-level validation is critical', 
         ha='center', fontsize=10, color='darkblue')

plt.tight_layout()
plt.savefig('/home/y411869/Projects/NDUFB7_HF_2026_04_20/03_results/figures/Figure2_crossplatform_v48.png', dpi=300, bbox_inches='tight')
plt.close()

print("✅ Figure 2已保存:")
print("  ~/Projects/NDUFB7_HF_2026_04_20/03_results/figures/Figure2_crossplatform_v48.png")
