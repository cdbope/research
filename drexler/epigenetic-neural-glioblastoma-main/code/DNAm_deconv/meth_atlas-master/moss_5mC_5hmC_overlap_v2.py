#!/usr/bin/env python3
"""
Compare Moss deconvolution results between:
- Cortical neuron-specific atlas (500 probes) vs Full Moss atlas (~7890 probes)

Both use 5mC+5hmC merged data.

Input files:
- Cortical-only: GBM_Moss_cortical_neurone_only_signatures.csv
  (using cortical_neuron_reference_atlas.csv - 500 specific probes)
- Full atlas: combined_deconv_betas_merged_deconv_output_cortical_neuron.csv
  (using original reference_atlas.csv - ~7890 probes)

Output: Confusion matrix, scatter plot, and per-sample comparison.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, cohen_kappa_score, accuracy_score
from scipy.stats import pearsonr, spearmanr

# File paths
CORTICAL_ONLY_FILE = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/GBM_Moss_cortical_neurone_only_signatures.csv'
FULL_ATLAS_FILE = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/combined_deconv_betas_merged_deconv_output_cortical_neuron.csv'
OUTPUT_DIR = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/'

# Load data
print("Loading data...")
cortical_only = pd.read_csv(CORTICAL_ONLY_FILE)
full_atlas = pd.read_csv(FULL_ATLAS_FILE)

print(f"Cortical-only atlas (500 probes) samples: {len(cortical_only)}")
print(f"Full Moss atlas (~7890 probes) samples: {len(full_atlas)}")

# Rename columns for clarity
cortical_only = cortical_only.rename(columns={
    'Cortical_neurons': 'Score_Cortical_Only',
    'Neural_Classification': 'Class_Cortical_Only'
})

full_atlas = full_atlas.rename(columns={
    'Cortical_neurons': 'Score_Full_Atlas',
    'Neural_Classification': 'Class_Full_Atlas'
})

# Merge dataframes
merged_df = pd.merge(cortical_only, full_atlas, on='Sample_id', how='inner')
merged_df['Agreement'] = merged_df['Class_Cortical_Only'] == merged_df['Class_Full_Atlas']
merged_df = merged_df.sort_values('Sample_id').reset_index(drop=True)

total = len(merged_df)
agree = merged_df['Agreement'].sum()
disagree = total - agree

print(f"Matched samples: {total}")

# Calculate counts
class_cortical_high = (merged_df['Class_Cortical_Only'] == 'High-Neural').sum()
class_cortical_low = (merged_df['Class_Cortical_Only'] == 'Low-Neural').sum()
class_full_high = (merged_df['Class_Full_Atlas'] == 'High-Neural').sum()
class_full_low = (merged_df['Class_Full_Atlas'] == 'Low-Neural').sum()

# Create confusion matrix
labels = ['Low-Neural', 'High-Neural']
cm = confusion_matrix(merged_df['Class_Cortical_Only'], merged_df['Class_Full_Atlas'], labels=labels)
accuracy = accuracy_score(merged_df['Class_Cortical_Only'], merged_df['Class_Full_Atlas'])
kappa = cohen_kappa_score(merged_df['Class_Cortical_Only'], merged_df['Class_Full_Atlas'])

# Calculate correlations
pearson_r, pearson_p = pearsonr(merged_df['Score_Cortical_Only'], merged_df['Score_Full_Atlas'])
spearman_r, spearman_p = spearmanr(merged_df['Score_Cortical_Only'], merged_df['Score_Full_Atlas'])

# Print results
print("\n" + "=" * 60)
print("CORTICAL-ONLY ATLAS vs FULL MOSS ATLAS COMPARISON")
print("(Both using 5mC+5hmC merged data)")
print("=" * 60)

print(f"\nTotal samples: {total}")
print(f"Agreement:     {agree} ({100*agree/total:.1f}%)")
print(f"Disagreement:  {disagree} ({100*disagree/total:.1f}%)")

print("\n" + "-" * 60)
print("CLASSIFICATION COUNTS")
print("-" * 60)
print(f"Cortical-only (500 probes) - High-Neural: {class_cortical_high} ({100*class_cortical_high/total:.1f}%)")
print(f"Cortical-only (500 probes) - Low-Neural:  {class_cortical_low} ({100*class_cortical_low/total:.1f}%)")
print(f"Full atlas (~7890 probes) - High-Neural:  {class_full_high} ({100*class_full_high/total:.1f}%)")
print(f"Full atlas (~7890 probes) - Low-Neural:   {class_full_low} ({100*class_full_low/total:.1f}%)")

print("\n" + "-" * 60)
print("CONFUSION MATRIX")
print("-" * 60)
print(f"                           Full Atlas")
print(f"                        Low      High")
print(f"Cortical-only  Low      {cm[0,0]:3d}      {cm[0,1]:3d}")
print(f"               High     {cm[1,0]:3d}      {cm[1,1]:3d}")
print(f"\nAccuracy: {accuracy:.2f}, Cohen's kappa: {kappa:.2f}")

print("\n" + "-" * 60)
print("SCORE CORRELATIONS")
print("-" * 60)
print(f"Pearson r:  {pearson_r:.3f} (p={pearson_p:.2e})")
print(f"Spearman r: {spearman_r:.3f} (p={spearman_p:.2e})")

# Print score statistics
print("\n" + "-" * 60)
print("SCORE STATISTICS")
print("-" * 60)
print(f"Cortical-only: mean={merged_df['Score_Cortical_Only'].mean():.3f}, median={merged_df['Score_Cortical_Only'].median():.3f}, range=[{merged_df['Score_Cortical_Only'].min():.3f}, {merged_df['Score_Cortical_Only'].max():.3f}]")
print(f"Full atlas:    mean={merged_df['Score_Full_Atlas'].mean():.3f}, median={merged_df['Score_Full_Atlas'].median():.3f}, range=[{merged_df['Score_Full_Atlas'].min():.3f}, {merged_df['Score_Full_Atlas'].max():.3f}]")

# Create visualization
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 1. Confusion matrix
ax1 = axes[0]
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax1,
            xticklabels=['Low-Neural', 'High-Neural'],
            yticklabels=['Low-Neural', 'High-Neural'],
            annot_kws={'size': 18, 'weight': 'bold'})
ax1.set_xlabel('Full Atlas Classification\n(~7890 probes)', fontsize=11)
ax1.set_ylabel('Cortical-Only Classification\n(500 probes)', fontsize=11)
ax1.set_title(f'Classification Agreement\n(n={total}, Agreement={100*agree/total:.1f}%, kappa={kappa:.2f})',
              fontsize=12, fontweight='bold')

# 2. Scatter plot of scores
ax2 = axes[1]

# Color by agreement
colors = []
for _, row in merged_df.iterrows():
    if row['Agreement']:
        if row['Class_Cortical_Only'] == 'High-Neural':
            colors.append('#522404')  # Both High
        else:
            colors.append('#EC9D65')  # Both Low
    else:
        colors.append('#FF0000')  # Disagreement

ax2.scatter(merged_df['Score_Cortical_Only'], merged_df['Score_Full_Atlas'],
            c=colors, alpha=0.7, s=30, edgecolors='none')

# Add diagonal line
max_val = max(merged_df['Score_Cortical_Only'].max(), merged_df['Score_Full_Atlas'].max())
min_val = min(merged_df['Score_Cortical_Only'].min(), merged_df['Score_Full_Atlas'].min())
ax2.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, label='y=x')

ax2.set_xlabel('Cortical-Only Score (500 probes)', fontsize=11)
ax2.set_ylabel('Full Atlas Score (~7890 probes)', fontsize=11)
ax2.set_title(f'Score Comparison\n(Pearson r={pearson_r:.2f}, Red=Disagreement)',
              fontsize=12, fontweight='bold')
ax2.legend(loc='upper left', fontsize=9)

# 3. Bar chart of categories
ax3 = axes[2]
categories = ['Both Low', 'Both High', 'Cortical Low\nFull High', 'Cortical High\nFull Low']
counts = [cm[0,0], cm[1,1], cm[0,1], cm[1,0]]
colors_bar = ['#EC9D65', '#522404', '#FF6B6B', '#FF6B6B']

bars = ax3.bar(categories, counts, color=colors_bar, edgecolor='black', linewidth=1)

for bar, count in zip(bars, counts):
    ax3.annotate(f'{count}\n({100*count/total:.1f}%)',
                xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                xytext=(0, 5), textcoords="offset points",
                ha='center', va='bottom', fontsize=10, fontweight='bold')

ax3.set_ylabel('Number of Samples', fontsize=12)
ax3.set_title('Classification Agreement by Category', fontsize=12, fontweight='bold')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#EC9D65', label='Both Low-Neural'),
                   Patch(facecolor='#522404', label='Both High-Neural'),
                   Patch(facecolor='#FF6B6B', label='Disagreement')]
ax3.legend(handles=legend_elements, loc='upper right')

plt.suptitle('Cortical-Only Atlas (500 probes) vs Full Moss Atlas (~7890 probes)\n(Both using 5mC+5hmC merged data)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()

# Save figure
output_png = OUTPUT_DIR + 'moss_cortical_only_vs_full_atlas.png'
output_pdf = OUTPUT_DIR + 'moss_cortical_only_vs_full_atlas.pdf'
plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {output_png}")
print(f"Saved: {output_pdf}")

# Save per-sample results
output_csv = OUTPUT_DIR + 'moss_cortical_only_vs_full_atlas_per_sample.csv'
merged_df.to_csv(output_csv, index=False)
print(f"Saved: {output_csv}")

# Print disagreeing samples
print("\n" + "-" * 60)
print("SAMPLES WITH DISAGREEMENT")
print("-" * 60)
disagree_df = merged_df[~merged_df['Agreement']]
for _, row in disagree_df.iterrows():
    print(f"{row['Sample_id']}: Cortical-only={row['Class_Cortical_Only']} ({row['Score_Cortical_Only']:.3f}), Full={row['Class_Full_Atlas']} ({row['Score_Full_Atlas']:.3f})")

plt.close()
print("\nDone!")
