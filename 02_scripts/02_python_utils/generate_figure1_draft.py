import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 读取六样本数据
csv_path = "/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/04_spatial_geo/Kuppe_Nature_2022/Kuppe_six_sample_NDUFB7.csv"
df = pd.read_csv(csv_path)

# 设置顺序：Healthy → MI（从低到高或逻辑顺序）
order = ['GT_P13', 'GT_P9', 'GT_P15', 'FZ_P20', 'IZ_BZ_P2', 'IZ_P3']
df_ordered = df.set_index('Sample').reindex(order).reset_index()

# 颜色映射
colors = ['#2E86AB', '#2E86AB', '#2E86AB', '#A23B72', '#F18F01', '#C73E1D']  # 蓝→紫→橙→红

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: 阳性率柱状图
ax1 = axes[0, 0]
bars = ax1.bar(range(len(df_ordered)), df_ordered['Pos_rate_pct'], color=colors, edgecolor='black', linewidth=0.5)
ax1.set_xticks(range(len(df_ordered)))
ax1.set_xticklabels(df_ordered['Sample'], rotation=45, ha='right')
ax1.set_ylabel('NDUFB7 Positive Rate (%)', fontsize=11)
ax1.set_title('A. NDUFB7 Detection Rate Across Regions', fontsize=12, fontweight='bold')
ax1.axhline(y=50, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
# 添加数值标签
for i, (bar, val) in enumerate(zip(bars, df_ordered['Pos_rate_pct'])):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
             f'{val:.1f}%', ha='center', va='bottom', fontsize=9)
# 添加分组标签
ax1.text(1, 95, 'Healthy', ha='center', fontsize=10, color='#2E86AB', fontweight='bold')
ax1.text(4, 95, 'Post-MI', ha='center', fontsize=10, color='#C73E1D', fontweight='bold')

# Panel B: 均值+中位数双指标
ax2 = axes[0, 1]
x = np.arange(len(df_ordered))
width = 0.35
bars1 = ax2.bar(x - width/2, df_ordered['Mean'], width, label='Mean', color=colors, alpha=0.8, edgecolor='black')
bars2 = ax2.bar(x + width/2, df_ordered['Median'], width, label='Median', color=colors, alpha=0.4, edgecolor='black')
ax2.set_xticks(x)
ax2.set_xticklabels(df_ordered['Sample'], rotation=45, ha='right')
ax2.set_ylabel('Expression (log-normalized)', fontsize=11)
ax2.set_title('B. NDUFB7 Expression Level (Mean vs Median)', fontsize=12, fontweight='bold')
ax2.legend()
ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.3)

# Panel C: 样本量与阳性数
ax3 = axes[1, 0]
scatter = ax3.scatter(df_ordered['N_spots'], df_ordered['Pos_rate_pct'], 
                     s=df_ordered['Positive']/10, c=colors, alpha=0.7, edgecolors='black')
ax3.set_xlabel('Total Spots', fontsize=11)
ax3.set_ylabel('Positive Rate (%)', fontsize=11)
ax3.set_title('C. Detection Rate vs Tissue Coverage', fontsize=12, fontweight='bold')
for i, row in df_ordered.iterrows():
    ax3.annotate(row['Sample'], (row['N_spots'], row['Pos_rate_pct']), 
                textcoords="offset points", xytext=(5, 5), fontsize=8)

# Panel D: 关键发现示意图（文本+箭头）
ax4 = axes[1, 1]
ax4.axis('off')
ax4.set_xlim(0, 10)
ax4.set_ylim(0, 10)

# 绘制发现流程图
ax4.text(5, 9, 'Key Findings (V46)', ha='center', fontsize=14, fontweight='bold')
ax4.text(5, 7.5, '• Acute Infarct Zone (IZ): NDUFB7 preserved\n  (64.9% positive, median=1.26)', 
         ha='center', fontsize=10, color='#C73E1D')
ax4.text(5, 5.5, '• Chronic Fibrotic Zone (FZ): NDUFB7 depleted\n  (40.2% positive, median=0.00)', 
         ha='center', fontsize=10, color='#A23B72')
ax4.text(5, 3.5, '• Healthy controls: extreme heterogeneity\n  (35.9% - 90.1% across patients)', 
         ha='center', fontsize=10, color='#2E86AB')
ax4.text(5, 1.5, '• Bulk GSE57338: no HF difference (p=0.52)\n  Spatial heterogeneity > tissue average', 
         ha='center', fontsize=10, color='black')

# 添加箭头
ax4.annotate('', xy=(8, 6.5), xytext=(8, 7.5),
            arrowprops=dict(arrowstyle='->', color='gray', lw=2))
ax4.text(8.5, 7, 'Acute\n→Chronic', fontsize=9, color='gray')

plt.tight_layout()
plt.savefig('/home/y411869/Projects/NDUFB7_HF_2026_04_20/03_results/figures/Figure1_draft_v46.png', dpi=300, bbox_inches='tight')
plt.close()

print("✅ Figure 1草图已保存:")
print("  ~/Projects/NDUFB7_HF_2026_04_20/03_results/figures/Figure1_draft_v46.png")
