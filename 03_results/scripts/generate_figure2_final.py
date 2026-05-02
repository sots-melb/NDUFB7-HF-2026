import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: 四平台NDUFB7相对表达（标准化到Healthy=1）
ax1 = axes[0, 0]
platforms = ['GSE57338\nAffy 1.1 ST\nn=313', 'GDS4772\nAffy 1.0 ST\nn=17', 
             'GSE116250\nRNA-seq RPKM\nn=64', 'GSE55296\nRNA-seq count\nn=37']
x = np.arange(len(platforms))
width = 0.25

# 相对值（Healthy=1）
dcm_rel = [1.009, 1.011, 1.146, 1.254]  # GSE55296: 1347.8/1075.1=1.254
icm_rel = [0.996, 0, 1.185, 0.926]       # GSE55296: 995.4/1075.1=0.926

bars1 = ax1.bar(x - width/2, dcm_rel, width, label='DCM', color='#E41A1C', alpha=0.8, edgecolor='black')
bars2 = ax1.bar(x + width/2, [icm_rel[0], 0, icm_rel[2], icm_rel[3]], width, label='Ischemic', color='#377EB8', alpha=0.8, edgecolor='black')
ax1.axhline(y=1.0, color='black', linestyle='--', linewidth=1)
ax1.set_xticks(x)
ax1.set_xticklabels(platforms, fontsize=9)
ax1.set_ylabel('Relative Expression (Healthy=1)', fontsize=11)
ax1.set_title('A. Cross-Platform NDUFB7: DCM vs Ischemic vs Healthy', fontsize=12, fontweight='bold')
ax1.legend()
ax1.set_ylim(0.85, 1.35)

# 添加显著性标记
ax1.text(0, 1.02, 'ns', ha='center', fontsize=8, color='gray')
ax1.text(1, 1.02, 'ns', ha='center', fontsize=8, color='gray')
ax1.text(2, 1.15, '**', ha='center', fontsize=10, color='#E41A1C')
ax1.text(3, 1.28, '*', ha='center', fontsize=10, color='#E41A1C')
ax1.text(3, 0.95, 'ns', ha='center', fontsize=8, color='gray')

# Panel B: 统计显著性汇总
ax2 = axes[0, 1]
comparisons = ['GSE57338\nDCM vs NF', 'GSE57338\nICM vs NF', 'GSE116250\nDCM vs NF', 
               'GSE116250\nICM vs NF', 'GDS4772\nDCM vs NF', 'GSE55296\nDCM vs H', 
               'GSE55296\nICM vs H', 'GSE55296\nDCM vs ICM']
p_values = [0.3, 0.3, 0.0081, 0.0348, 0.4421, 0.2036, 0.9753, 0.0455]
colors_p = ['gray', 'gray', '#E41A1C', '#E41A1C', 'gray', 'orange', 'gray', '#E41A1C']

bars = ax2.barh(range(len(comparisons)), [-np.log10(max(p, 1e-10)) for p in p_values], 
                color=colors_p, edgecolor='black')
ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
ax2.set_yticks(range(len(comparisons)))
ax2.set_yticklabels(comparisons, fontsize=8)
ax2.set_xlabel('-log10(p-value)', fontsize=11)
ax2.set_title('B. Statistical Significance Across 4 Platforms', fontsize=12, fontweight='bold')
ax2.legend(loc='lower right')

# Panel C: RPKM偏倚示意图
ax3 = axes[1, 0]
genes = ['NDUFB7\n(137aa)', 'MT-CO1\n(513aa)', 'NDUFB8\n(156aa)', 'NDUFS1\n(727aa)']
lengths = [137, 513, 156, 727]
# 相同reads下，RPKM与长度成反比
rpkm_bias = [max(lengths)/l for l in lengths]
ax3.bar(range(len(genes)), rpkm_bias, color=['#E41A1C', '#4DAF4A', '#FF7F00', '#984EA3'], 
        edgecolor='black', linewidth=1)
ax3.set_xticks(range(len(genes)))
ax3.set_xticklabels(genes, fontsize=9)
ax3.set_ylabel('RPKM Inflation Factor\n(same reads)', fontsize=10)
ax3.set_title('C. RPKM Length Normalization Bias\n(137-aa NDUFB7 artificially inflated 5.3x)', 
              fontsize=12, fontweight='bold')
for i, v in enumerate(rpkm_bias):
    ax3.text(i, v + 0.1, f'{v:.1f}x', ha='center', fontsize=10, fontweight='bold')

# Panel D: 证据整合总结
ax4 = axes[1, 1]
ax4.axis('off')
ax4.set_xlim(0, 10)
ax4.set_ylim(0, 10)

ax4.text(5, 9, 'Evidence Integration (V49)', ha='center', fontsize=14, fontweight='bold')

ax4.text(5, 7.3, 'Robust Findings:\n'
         '• Spatial: FZ depleted vs IZ preserved (p<0.0001)\n'
         '• Protein: Detectable in human heart (12 peptides)\n'
         '• Bulk: No overall HF downregulation (4 platforms)', 
         ha='center', fontsize=10, color='darkgreen')

ax4.text(5, 5.0, 'Platform-Dependent:\n'
         '• Affymetrix: DCM trend > NF, no significance\n'
         '• RNA-seq count: DCM > Ischemic (p=0.046)\n'
         '• RNA-seq RPKM: Artificial upregulation (bias)', 
         ha='center', fontsize=10, color='darkorange')

ax4.text(5, 2.5, 'Conclusion:\n'
         'NDUFB7 pathophysiology = spatial redistribution\n'
         '(fibrotic zone loss) + etiology-specific gradient\n'
         '(DCM retention), not uniform transcriptional loss', 
         ha='center', fontsize=10, color='darkblue', fontweight='bold')

plt.tight_layout()
plt.savefig('/home/y411869/Projects/NDUFB7_HF_2026_04_20/03_results/figures/Figure2_crossplatform_v49.png', 
            dpi=300, bbox_inches='tight')
plt.close()

print("✅ Figure 2 (v49) 已保存:")
print("  ~/Projects/NDUFB7_HF_2026_04_20/03_results/figures/Figure2_crossplatform_v49.png")
