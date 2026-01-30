#!/usr/bin/env python3
"""
Compare Drexler vs Moss neural classification per sample.

Shows per-sample agreement/disagreement between:
- Drexler: Neural signature probability-based classification
- Moss: Cortical neurons deconvolution-based classification
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, cohen_kappa_score, accuracy_score

# File paths
DREXLER_FILE = '/home/chbope/extension/script/drexler/results/200gbm/mergemh/proba//neural_classification_summary.tsv'
MOSS_FILE = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/combined_deconv_betas_merged_deconv_output_cortical_neuron.csv'
OUTPUT_DIR = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/'

# Load data
print("Loading data...")
drexler_df = pd.read_csv(DREXLER_FILE, sep='\t')
moss_df = pd.read_csv(MOSS_FILE)

# Rename columns for clarity
drexler_df = drexler_df.rename(columns={
    'sample_id': 'Sample_id',
    'Neural_Classification': 'Drexler_Class'
})

moss_df = moss_df.rename(columns={
    'Neural_Classification': 'Moss_Class'
})

# Merge dataframes
merged_df = pd.merge(drexler_df[['Sample_id', 'Drexler_Class']],
                      moss_df[['Sample_id', 'Moss_Class']],
                      on='Sample_id', how='inner')

# Add agreement column
merged_df['Agreement'] = merged_df['Drexler_Class'] == merged_df['Moss_Class']

# Sort by sample ID
merged_df = merged_df.sort_values('Sample_id').reset_index(drop=True)

total = len(merged_df)
agree = merged_df['Agreement'].sum()
disagree = total - agree

# Calculate counts
drexler_high = (merged_df['Drexler_Class'] == 'High-Neural').sum()
drexler_low = (merged_df['Drexler_Class'] == 'Low-Neural').sum()
moss_high = (merged_df['Moss_Class'] == 'High-Neural').sum()
moss_low = (merged_df['Moss_Class'] == 'Low-Neural').sum()

# Create confusion matrix
labels = ['Low-Neural', 'High-Neural']
cm = confusion_matrix(merged_df['Drexler_Class'], merged_df['Moss_Class'], labels=labels)
accuracy = accuracy_score(merged_df['Drexler_Class'], merged_df['Moss_Class'])
kappa = cohen_kappa_score(merged_df['Drexler_Class'], merged_df['Moss_Class'])

# Print results
print("\n" + "=" * 60)
print("PER-SAMPLE CLASSIFICATION COMPARISON")
print("=" * 60)
print(f"\nTotal samples: {total}")
print(f"Agreement:     {agree} ({100*agree/total:.1f}%)")
print(f"Disagreement:  {disagree} ({100*disagree/total:.1f}%)")

print("\n" + "-" * 60)
print("CLASSIFICATION COUNTS")
print("-" * 60)
print(f"Drexler - High-Neural: {drexler_high} ({100*drexler_high/total:.1f}%)")
print(f"Drexler - Low-Neural:  {drexler_low} ({100*drexler_low/total:.1f}%)")
print(f"Moss    - High-Neural: {moss_high} ({100*moss_high/total:.1f}%)")
print(f"Moss    - Low-Neural:  {moss_low} ({100*moss_low/total:.1f}%)")

print("\n" + "-" * 60)
print("CONFUSION MATRIX")
print("-" * 60)
print(f"                        Moss")
print(f"                   Low      High")
print(f"Drexler  Low       {cm[0,0]:3d}      {cm[0,1]:3d}")
print(f"         High      {cm[1,0]:3d}      {cm[1,1]:3d}")
print(f"\nAccuracy: {accuracy:.2f}, Cohen's κ: {kappa:.2f}")

# Print disagreeing samples
print("\n" + "-" * 60)
print("SAMPLES WITH DISAGREEMENT")
print("-" * 60)
disagree_df = merged_df[~merged_df['Agreement']]
for _, row in disagree_df.iterrows():
    print(f"{row['Sample_id']}: Drexler={row['Drexler_Class']}, Moss={row['Moss_Class']}")

# Create visualization
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 1. Confusion matrix
ax1 = axes[0]
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax1,
            xticklabels=['Low-Neural', 'High-Neural'],
            yticklabels=['Low-Neural', 'High-Neural'],
            annot_kws={'size': 20, 'weight': 'bold'})
ax1.set_xlabel('Moss Classification', fontsize=12)
ax1.set_ylabel('Drexler Classification', fontsize=12)
ax1.set_title(f'Per-Sample Agreement\n(n={total}, Agreement={100*agree/total:.1f}%, κ={kappa:.2f})',
              fontsize=12, fontweight='bold')

# 2. Per-sample comparison strip plot
ax2 = axes[1]

# Create data for plotting - each sample gets two points (Drexler and Moss)
plot_data = []
for idx, row in merged_df.iterrows():
    # Drexler
    plot_data.append({
        'Sample': row['Sample_id'],
        'Method': 'Drexler',
        'Classification': 1 if row['Drexler_Class'] == 'High-Neural' else 0,
        'Agreement': row['Agreement']
    })
    # Moss
    plot_data.append({
        'Sample': row['Sample_id'],
        'Method': 'Moss',
        'Classification': 1 if row['Moss_Class'] == 'High-Neural' else 0,
        'Agreement': row['Agreement']
    })

plot_df = pd.DataFrame(plot_data)

# Create summary bar showing agreement/disagreement
categories = ['Both Low', 'Both High', 'Drexler Low\nMoss High', 'Drexler High\nMoss Low']
counts = [cm[0,0], cm[1,1], cm[0,1], cm[1,0]]
colors = ['#EC9D65', '#522404', '#FF6B6B', '#FF6B6B']

bars = ax2.bar(categories, counts, color=colors, edgecolor='black', linewidth=1)

# Add count labels
for bar, count in zip(bars, counts):
    ax2.annotate(f'{count}\n({100*count/total:.1f}%)',
                xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                xytext=(0, 5), textcoords="offset points",
                ha='center', va='bottom', fontsize=11, fontweight='bold')

ax2.set_ylabel('Number of Samples', fontsize=12)
ax2.set_title('Classification Agreement by Category', fontsize=12, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#EC9D65', label='Both Low-Neural'),
                   Patch(facecolor='#522404', label='Both High-Neural'),
                   Patch(facecolor='#FF6B6B', label='Disagreement')]
ax2.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()

# Save figure
output_png = OUTPUT_DIR + 'confusion_drexler_moss_proba.png'
output_pdf = OUTPUT_DIR + 'confusion_drexler_moss_proba.pdf'
plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {output_png}")
print(f"Saved: {output_pdf}")

# Save detailed per-sample results
output_csv = OUTPUT_DIR + 'confusion_drexler_moss_per_sample_proba.csv'
merged_df.to_csv(output_csv, index=False)
print(f"Saved: {output_csv}")

plt.close()
print("\nDone!")
