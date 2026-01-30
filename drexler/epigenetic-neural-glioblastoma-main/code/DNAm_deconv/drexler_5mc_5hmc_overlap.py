#!/usr/bin/env python3
"""
Compare Drexler classification results between 5mC only vs 5mC+5hmC merged.

Input files:
- 5mC only: Drexler_219_GBM_aggregated_predictions.csv
- 5mC+5hmC: drexler_5mc_5hmc_neural_classification_summary.tsv

Output: Confusion matrix, scatter plot, and per-sample comparison.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, cohen_kappa_score, accuracy_score
from scipy.stats import pearsonr, spearmanr

# File paths
DREXLER_5MC_FILE = '/home/chbope/extension/script/drexler/results/200gbm/prediction/individual/Drexler_219_GBM_aggregated_predictions.csv'
DREXLER_5MC_5HMC_FILE = '/home/chbope/extension/script/drexler/results/200gbm/mergemh/proba/drexler_5mc_5hmc_neural_classification_summary.tsv'
OUTPUT_DIR = '/home/chbope/extension/script/drexler/results/200gbm/mergemh/proba/'

# Load data
print("Loading data...")
drexler_5mc = pd.read_csv(DREXLER_5MC_FILE, index_col=0)
drexler_5mc_5hmc = pd.read_csv(DREXLER_5MC_5HMC_FILE, sep='\t')

# Reset index to make sample_id a column for 5mC data
drexler_5mc = drexler_5mc.reset_index()
drexler_5mc.columns = ['Sample_id', 'Prediction', 'Prediction_score', 'Neural_Classification']

# Rename columns for 5mC+5hmC data
drexler_5mc_5hmc = drexler_5mc_5hmc.rename(columns={'sample_id': 'Sample_id'})

print(f"5mC only samples: {len(drexler_5mc)}")
print(f"5mC+5hmC merged samples: {len(drexler_5mc_5hmc)}")

# Rename columns for clarity before merge
drexler_5mc = drexler_5mc.rename(columns={
    'Prediction_score': 'Score_5mC',
    'Neural_Classification': 'Class_5mC'
})

drexler_5mc_5hmc = drexler_5mc_5hmc.rename(columns={
    'Prediction_score': 'Score_5mC_5hmC',
    'Neural_Classification': 'Class_5mC_5hmC'
})

# Merge dataframes
merged_df = pd.merge(
    drexler_5mc[['Sample_id', 'Score_5mC', 'Class_5mC']],
    drexler_5mc_5hmc[['Sample_id', 'Score_5mC_5hmC', 'Class_5mC_5hmC']],
    on='Sample_id', how='inner'
)
merged_df['Agreement'] = merged_df['Class_5mC'] == merged_df['Class_5mC_5hmC']
merged_df = merged_df.sort_values('Sample_id').reset_index(drop=True)

total = len(merged_df)
agree = merged_df['Agreement'].sum()
disagree = total - agree

print(f"Matched samples: {total}")

# Calculate counts
class_5mc_high = (merged_df['Class_5mC'] == 'High-Neural').sum()
class_5mc_low = (merged_df['Class_5mC'] == 'Low-Neural').sum()
class_merged_high = (merged_df['Class_5mC_5hmC'] == 'High-Neural').sum()
class_merged_low = (merged_df['Class_5mC_5hmC'] == 'Low-Neural').sum()

# Create confusion matrix
labels = ['Low-Neural', 'High-Neural']
cm = confusion_matrix(merged_df['Class_5mC'], merged_df['Class_5mC_5hmC'], labels=labels)
accuracy = accuracy_score(merged_df['Class_5mC'], merged_df['Class_5mC_5hmC'])
kappa = cohen_kappa_score(merged_df['Class_5mC'], merged_df['Class_5mC_5hmC'])

# Calculate correlations
pearson_r, pearson_p = pearsonr(merged_df['Score_5mC'], merged_df['Score_5mC_5hmC'])
spearman_r, spearman_p = spearmanr(merged_df['Score_5mC'], merged_df['Score_5mC_5hmC'])

# Print results
print("\n" + "=" * 60)
print("DREXLER 5mC vs 5mC+5hmC COMPARISON")
print("=" * 60)

print(f"\nTotal samples: {total}")
print(f"Agreement:     {agree} ({100*agree/total:.1f}%)")
print(f"Disagreement:  {disagree} ({100*disagree/total:.1f}%)")

print("\n" + "-" * 60)
print("CLASSIFICATION COUNTS")
print("-" * 60)
print(f"5mC only    - High-Neural: {class_5mc_high} ({100*class_5mc_high/total:.1f}%)")
print(f"5mC only    - Low-Neural:  {class_5mc_low} ({100*class_5mc_low/total:.1f}%)")
print(f"5mC+5hmC    - High-Neural: {class_merged_high} ({100*class_merged_high/total:.1f}%)")
print(f"5mC+5hmC    - Low-Neural:  {class_merged_low} ({100*class_merged_low/total:.1f}%)")

print("\n" + "-" * 60)
print("CONFUSION MATRIX")
print("-" * 60)
print(f"                        5mC+5hmC")
print(f"                   Low      High")
print(f"5mC      Low       {cm[0,0]:3d}      {cm[0,1]:3d}")
print(f"         High      {cm[1,0]:3d}      {cm[1,1]:3d}")
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
print(f"5mC only:   mean={merged_df['Score_5mC'].mean():.3f}, median={merged_df['Score_5mC'].median():.3f}, range=[{merged_df['Score_5mC'].min():.3f}, {merged_df['Score_5mC'].max():.3f}]")
print(f"5mC+5hmC:   mean={merged_df['Score_5mC_5hmC'].mean():.3f}, median={merged_df['Score_5mC_5hmC'].median():.3f}, range=[{merged_df['Score_5mC_5hmC'].min():.3f}, {merged_df['Score_5mC_5hmC'].max():.3f}]")

# Create visualization
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 1. Confusion matrix
ax1 = axes[0]
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax1,
            xticklabels=['Low-Neural', 'High-Neural'],
            yticklabels=['Low-Neural', 'High-Neural'],
            annot_kws={'size': 18, 'weight': 'bold'})
ax1.set_xlabel('5mC+5hmC Classification', fontsize=12)
ax1.set_ylabel('5mC only Classification', fontsize=12)
ax1.set_title(f'Classification Agreement\n(n={total}, Agreement={100*agree/total:.1f}%, kappa={kappa:.2f})',
              fontsize=12, fontweight='bold')

# 2. Scatter plot of scores
ax2 = axes[1]

# Color by agreement
colors = []
for _, row in merged_df.iterrows():
    if row['Agreement']:
        if row['Class_5mC'] == 'High-Neural':
            colors.append('#522404')  # Both High
        else:
            colors.append('#EC9D65')  # Both Low
    else:
        colors.append('#FF0000')  # Disagreement

ax2.scatter(merged_df['Score_5mC'], merged_df['Score_5mC_5hmC'],
            c=colors, alpha=0.7, s=30, edgecolors='none')

# Add diagonal line
ax2.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='y=x')

# Add threshold line at 0.5 for both axes
ax2.axvline(x=0.5, color='blue', linestyle=':', linewidth=1, alpha=0.5)
ax2.axhline(y=0.5, color='green', linestyle=':', linewidth=1, alpha=0.5)

ax2.set_xlabel('5mC only Score (Prediction probability)', fontsize=11)
ax2.set_ylabel('5mC+5hmC Score (Prediction probability)', fontsize=11)
ax2.set_title(f'Score Comparison\n(Pearson r={pearson_r:.2f}, Red=Disagreement)',
              fontsize=12, fontweight='bold')
ax2.set_xlim(0.4, 1.0)
ax2.set_ylim(0.4, 1.0)
ax2.legend(loc='upper left', fontsize=9)

# 3. Bar chart of categories
ax3 = axes[2]
categories = ['Both Low', 'Both High', '5mC Low\n5mC+5hmC High', '5mC High\n5mC+5hmC Low']
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

plt.suptitle('Drexler Classification: 5mC only vs 5mC+5hmC merged',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()

# Save figure
output_png = OUTPUT_DIR + 'drexler_5mC_5hmC_overlap.png'
output_pdf = OUTPUT_DIR + 'drexler_5mC_5hmC_overlap.pdf'
plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {output_png}")
print(f"Saved: {output_pdf}")

# Save per-sample results
output_csv = OUTPUT_DIR + 'drexler_5mC_5hmC_overlap_per_sample.csv'
merged_df.to_csv(output_csv, index=False)
print(f"Saved: {output_csv}")

# Print disagreeing samples
print("\n" + "-" * 60)
print("SAMPLES WITH DISAGREEMENT")
print("-" * 60)
disagree_df = merged_df[~merged_df['Agreement']]
for _, row in disagree_df.iterrows():
    print(f"{row['Sample_id']}: 5mC={row['Class_5mC']} ({row['Score_5mC']:.3f}), 5mC+5hmC={row['Class_5mC_5hmC']} ({row['Score_5mC_5hmC']:.3f})")

plt.close()
print("\nDone!")
