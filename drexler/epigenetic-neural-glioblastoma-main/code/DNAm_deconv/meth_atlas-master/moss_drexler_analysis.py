#!/usr/bin/env python3
"""
Comparative analysis between Moss Deconvolution and Drexler Logistic Regression
neural classification results.

This script compares:
1. Moss atlas deconvolution (Cortical_neurons proportion with median threshold)
   - Cortical_neurons: actual neural proportion (0-1 scale)
   - Classification based on median threshold (0.229)

2. Drexler et al. logistic regression (1289 CpG model)
   - Prediction: 0=Low-Neural, 1=High-Neural
   - Prediction_score: probability/confidence of the predicted class (NOT a neural score)

Usage:
    python moss_drexler_analysis.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read both files
moss_file = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/GBM_Moss_signatures.csv'
drexler_file = '/home/chbope/extension/script/drexler/results/200gbm/prediction/individual/Drexler_219_GBM_aggregated_predictions.csv'

moss = pd.read_csv(moss_file, index_col=0)
drexler = pd.read_csv(drexler_file, index_col=0)

print("=" * 70)
print("COMPARATIVE ANALYSIS: Moss Deconvolution vs Drexler Logistic Regression")
print("=" * 70)

print(f"\n1. DATASET OVERVIEW")
print("-" * 40)
print(f"   Moss Deconvolution samples:    {len(moss)}")
print(f"   Drexler LR Prediction samples: {len(drexler)}")

# Find overlap
common_samples = moss.index.intersection(drexler.index)
moss_only = moss.index.difference(drexler.index)
drexler_only = drexler.index.difference(moss.index)

print(f"\n2. SAMPLE OVERLAP")
print("-" * 40)
print(f"   Common samples:     {len(common_samples)}")
print(f"   Only in Moss:       {len(moss_only)}")
print(f"   Only in Drexler:    {len(drexler_only)}")

if len(moss_only) > 0:
    print(f"\n   Samples only in Moss: {list(moss_only)}")
if len(drexler_only) > 0:
    print(f"\n   Samples only in Drexler: {list(drexler_only)}")

# Merge on common samples
merged = moss.loc[common_samples].join(drexler.loc[common_samples], lsuffix='_Moss', rsuffix='_Drexler')

print(f"\n3. CLASSIFICATION COMPARISON (n={len(merged)})")
print("-" * 40)

# Count classifications
moss_high = (merged['Neural_Classification_Moss'] == 'High-Neural').sum()
moss_low = (merged['Neural_Classification_Moss'] == 'Low-Neural').sum()
drex_high = (merged['Neural_Classification_Drexler'] == 'High-Neural').sum()
drex_low = (merged['Neural_Classification_Drexler'] == 'Low-Neural').sum()

print(f"\n   Moss Deconvolution:")
print(f"     High-Neural: {moss_high} ({100*moss_high/len(merged):.1f}%)")
print(f"     Low-Neural:  {moss_low} ({100*moss_low/len(merged):.1f}%)")

print(f"\n   Drexler LR Prediction:")
print(f"     High-Neural: {drex_high} ({100*drex_high/len(merged):.1f}%)")
print(f"     Low-Neural:  {drex_low} ({100*drex_low/len(merged):.1f}%)")

# Agreement analysis
agree = (merged['Neural_Classification_Moss'] == merged['Neural_Classification_Drexler']).sum()
disagree = len(merged) - agree

print(f"\n4. CLASSIFICATION AGREEMENT")
print("-" * 40)
print(f"   Agreement:    {agree} ({100*agree/len(merged):.1f}%)")
print(f"   Disagreement: {disagree} ({100*disagree/len(merged):.1f}%)")

# Confusion matrix style
both_high = ((merged['Neural_Classification_Moss'] == 'High-Neural') &
             (merged['Neural_Classification_Drexler'] == 'High-Neural')).sum()
both_low = ((merged['Neural_Classification_Moss'] == 'Low-Neural') &
            (merged['Neural_Classification_Drexler'] == 'Low-Neural')).sum()
moss_high_drex_low = ((merged['Neural_Classification_Moss'] == 'High-Neural') &
                       (merged['Neural_Classification_Drexler'] == 'Low-Neural')).sum()
moss_low_drex_high = ((merged['Neural_Classification_Moss'] == 'Low-Neural') &
                       (merged['Neural_Classification_Drexler'] == 'High-Neural')).sum()

print(f"\n5. CONFUSION MATRIX")
print("-" * 40)
print(f"                          Drexler Prediction")
print(f"                       High-Neural  Low-Neural")
print(f"   Moss      High-Neural    {both_high:3d}         {moss_high_drex_low:3d}")
print(f"   Deconv    Low-Neural     {moss_low_drex_high:3d}         {both_low:3d}")

# Detailed disagreements
print(f"\n6. DISAGREEMENT DETAILS")
print("-" * 40)

disagreements = merged[merged['Neural_Classification_Moss'] != merged['Neural_Classification_Drexler']]

print(f"\n   a) Moss=High, Drexler=Low ({moss_high_drex_low} samples):")
subset = disagreements[disagreements['Neural_Classification_Moss'] == 'High-Neural']
if len(subset) > 0:
    for idx, row in subset.iterrows():
        print(f"      {idx}: Moss_Cortical_neurons={row['Cortical_neurons']:.3f}, Drexler_prob={row['Prediction_score']:.3f}")

print(f"\n   b) Moss=Low, Drexler=High ({moss_low_drex_high} samples):")
subset = disagreements[disagreements['Neural_Classification_Moss'] == 'Low-Neural']
if len(subset) > 0:
    for idx, row in subset.iterrows():
        print(f"      {idx}: Moss_Cortical_neurons={row['Cortical_neurons']:.3f}, Drexler_prob={row['Prediction_score']:.3f}")

# Cohen's Kappa
po = agree / len(merged)  # observed agreement
moss_high_pct = moss_high / len(merged)
moss_low_pct = moss_low / len(merged)
drex_high_pct = drex_high / len(merged)
drex_low_pct = drex_low / len(merged)
pe = (moss_high_pct * drex_high_pct) + (moss_low_pct * drex_low_pct)  # expected agreement
kappa = (po - pe) / (1 - pe) if pe < 1 else 0

print(f"\n7. STATISTICAL METRICS")
print("-" * 40)
print(f"   Cohen's Kappa: {kappa:.3f}")
if kappa < 0:
    interpretation = "No agreement"
elif kappa < 0.20:
    interpretation = "Slight agreement"
elif kappa < 0.40:
    interpretation = "Fair agreement"
elif kappa < 0.60:
    interpretation = "Moderate agreement"
elif kappa < 0.80:
    interpretation = "Substantial agreement"
else:
    interpretation = "Almost perfect agreement"
print(f"   Interpretation: {interpretation}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
The two methods show {100*agree/len(merged):.1f}% agreement ({agree}/{len(merged)} samples).

Key differences:
- Moss deconvolution (median threshold): {moss_high} High-Neural, {moss_low} Low-Neural
- Drexler LR (1289 CpG model): {drex_high} High-Neural, {drex_low} Low-Neural

The {disagree} disagreements ({100*disagree/len(merged):.1f}%) may reflect:
1. Different underlying methodologies (deconvolution vs. logistic regression)
2. Different probe sets used (Moss atlas ~7,890 vs. Drexler 1,289 CpGs)
3. Threshold effects for samples near decision boundaries

Note: Moss Cortical_neurons is a neural proportion (0-1), while Drexler
Prediction_score is the probability of the predicted class (not comparable)
""")

# Save merged comparison file
output_file = 'moss_drexler_comparison.csv'
merged.to_csv(output_file)
print(f"\nSaved comparison file: {output_file}")

# ============================================================================
# VISUALIZATIONS
# ============================================================================
print("\n" + "=" * 70)
print("GENERATING VISUALIZATIONS")
print("=" * 70)

# Define colors
colors = {'High-Neural': '#522404', 'Low-Neural': '#EC9D65'}
agree_color = '#4CAF50'  # green
disagree_color = '#F44336'  # red

# Create figure with multiple subplots
fig = plt.figure(figsize=(16, 12))

# --- Plot 1: Confusion Matrix Heatmap ---
ax1 = fig.add_subplot(2, 3, 1)
conf_matrix = np.array([[both_high, moss_high_drex_low],
                        [moss_low_drex_high, both_low]])
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', ax=ax1,
            xticklabels=['High-Neural', 'Low-Neural'],
            yticklabels=['High-Neural', 'Low-Neural'],
            annot_kws={'size': 14, 'weight': 'bold'})
ax1.set_xlabel('Drexler LR Prediction', fontsize=12)
ax1.set_ylabel('Moss Deconvolution', fontsize=12)
ax1.set_title(f'Confusion Matrix\n(Agreement: {100*agree/len(merged):.1f}%)', fontsize=14, fontweight='bold')

# --- Plot 2: Agreement Pie Chart ---
ax2 = fig.add_subplot(2, 3, 2)
agreement_counts = [agree, disagree]
agreement_labels = [f'Agree\n({agree})', f'Disagree\n({disagree})']
wedges, texts, autotexts = ax2.pie(agreement_counts, labels=agreement_labels,
                                    colors=[agree_color, disagree_color],
                                    autopct='%1.1f%%', startangle=90,
                                    explode=(0, 0.05))
for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontweight('bold')
ax2.set_title(f'Classification Agreement\n(n={len(merged)})', fontsize=14, fontweight='bold')

# --- Plot 3: Side-by-side Bar Chart ---
ax3 = fig.add_subplot(2, 3, 3)
x = np.arange(2)
width = 0.35
bars1 = ax3.bar(x - width/2, [moss_high, moss_low], width, label='Moss Deconv',
                color=['#522404', '#EC9D65'], edgecolor='black', linewidth=1.5)
bars2 = ax3.bar(x + width/2, [drex_high, drex_low], width, label='Drexler LR',
                color=['#522404', '#EC9D65'], edgecolor='black', linewidth=1.5,
                hatch='///', alpha=0.7)
ax3.set_ylabel('Number of Samples', fontsize=12)
ax3.set_xticks(x)
ax3.set_xticklabels(['High-Neural', 'Low-Neural'])
ax3.legend()
ax3.set_title('Classification Counts by Method', fontsize=14, fontweight='bold')
# Add count labels on bars
for bar in bars1:
    height = bar.get_height()
    ax3.annotate(f'{int(height)}', xy=(bar.get_x() + bar.get_width()/2, height),
                 xytext=(0, 3), textcoords='offset points', ha='center', fontweight='bold')
for bar in bars2:
    height = bar.get_height()
    ax3.annotate(f'{int(height)}', xy=(bar.get_x() + bar.get_width()/2, height),
                 xytext=(0, 3), textcoords='offset points', ha='center', fontweight='bold')

# --- Plot 4: Moss Cortical_neurons Distribution by Drexler Classification ---
ax4 = fig.add_subplot(2, 3, 4)
drex_high_samples = merged[merged['Neural_Classification_Drexler'] == 'High-Neural']['Cortical_neurons']
drex_low_samples = merged[merged['Neural_Classification_Drexler'] == 'Low-Neural']['Cortical_neurons']

# Box plot
bp = ax4.boxplot([drex_high_samples, drex_low_samples],
                  labels=['Drexler\nHigh-Neural', 'Drexler\nLow-Neural'],
                  patch_artist=True)
bp['boxes'][0].set_facecolor('#522404')
bp['boxes'][0].set_alpha(0.7)
bp['boxes'][1].set_facecolor('#EC9D65')
bp['boxes'][1].set_alpha(0.7)

# Add individual points
np.random.seed(42)
for i, data in enumerate([drex_high_samples, drex_low_samples], 1):
    x_jitter = np.random.normal(i, 0.04, len(data))
    ax4.scatter(x_jitter, data, alpha=0.5, s=20, c='black')

ax4.axhline(y=0.229, color='red', linestyle='--', linewidth=1.5, label='Moss threshold (0.229)')
ax4.set_ylabel('Moss Cortical_neurons', fontsize=12)
ax4.set_title('Moss Neural Proportion\nby Drexler Classification', fontsize=14, fontweight='bold')
ax4.legend(loc='upper right', fontsize=9)

# --- Plot 5: Scatter plot with agreement highlighting ---
ax5 = fig.add_subplot(2, 3, 5)
# Create agreement column
merged['Agreement'] = merged['Neural_Classification_Moss'] == merged['Neural_Classification_Drexler']

# Plot disagreements first (so they appear behind)
disagree_data = merged[~merged['Agreement']]
agree_data = merged[merged['Agreement']]

ax5.scatter(disagree_data['Cortical_neurons'], disagree_data.index,
            c=disagree_color, s=50, alpha=0.7, label=f'Disagree (n={len(disagree_data)})', marker='x')
ax5.scatter(agree_data['Cortical_neurons'], agree_data.index,
            c=agree_color, s=30, alpha=0.5, label=f'Agree (n={len(agree_data)})', marker='o')

ax5.axvline(x=0.229, color='black', linestyle='--', linewidth=1.5, label='Moss threshold')
ax5.set_xlabel('Moss Cortical_neurons', fontsize=12)
ax5.set_ylabel('Sample', fontsize=12)
ax5.set_title('Sample Distribution\n(Agreement vs Disagreement)', fontsize=14, fontweight='bold')
ax5.legend(loc='upper right', fontsize=9)
ax5.set_yticks([])

# --- Plot 6: Stacked bar showing agreement breakdown ---
ax6 = fig.add_subplot(2, 3, 6)
categories = ['Both\nHigh-Neural', 'Both\nLow-Neural', 'Moss High\nDrexler Low', 'Moss Low\nDrexler High']
counts = [both_high, both_low, moss_high_drex_low, moss_low_drex_high]
bar_colors = [agree_color, agree_color, disagree_color, disagree_color]
bars = ax6.bar(categories, counts, color=bar_colors, edgecolor='black', linewidth=1.5)

# Add count labels
for bar, count in zip(bars, counts):
    height = bar.get_height()
    ax6.annotate(f'{count}\n({100*count/len(merged):.1f}%)',
                 xy=(bar.get_x() + bar.get_width()/2, height),
                 xytext=(0, 3), textcoords='offset points', ha='center',
                 fontweight='bold', fontsize=10)

ax6.set_ylabel('Number of Samples', fontsize=12)
ax6.set_title('Detailed Agreement Breakdown', fontsize=14, fontweight='bold')

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=agree_color, edgecolor='black', label='Agreement'),
                   Patch(facecolor=disagree_color, edgecolor='black', label='Disagreement')]
ax6.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()

# Save figure
plt.savefig('moss_drexler_comparison_plots.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('moss_drexler_comparison_plots.pdf', bbox_inches='tight', facecolor='white')
print("\nSaved plots:")
print("  moss_drexler_comparison_plots.png")
print("  moss_drexler_comparison_plots.pdf")

plt.close()

print("\n" + "=" * 70)
print("DONE!")
print("=" * 70)
