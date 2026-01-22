#!/usr/bin/env python3
"""
Threshold extraction and classification plot for GBM Moss signatures.

This script:
1. Calculates median threshold from Cortical_neurons values (following Drexler et al. approach)
2. Classifies samples as High-Neural (>= median) or Low-Neural (< median)
3. Creates a beeswarm-like scatter plot showing the distribution

Usage:
    python threshold_mean_extraction_and_plot.py

Based on Drexler et al. (Nature Medicine 2024):
"The combined dataset was dichotomized into low- and high-neural tumors
using the median neural proportion of 0.41."
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the data
input_file = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/GBM_Moss_signatures.csv'
df = pd.read_csv(input_file, index_col=0)

# Calculate statistics
n_total = len(df)
scores = df['Cortical_neurons']

# Calculate median threshold (following Drexler et al. approach)
median_threshold = scores.median()
mean_value = scores.mean()
std_value = scores.std()
min_value = scores.min()
max_value = scores.max()

print("=" * 60)
print("GBM Moss Signatures - Threshold Extraction")
print("=" * 60)
print(f"Total samples: {n_total}")
print(f"\nCortical_neurons statistics:")
print(f"  Mean: {mean_value:.4f}")
print(f"  Median: {median_threshold:.4f}")
print(f"  Std: {std_value:.4f}")
print(f"  Min: {min_value:.4f}")
print(f"  Max: {max_value:.4f}")

# Save threshold to file
with open('gbm_threshold.txt', 'w') as f:
    f.write(f"GBM Neural Signature Threshold\n")
    f.write(f"==============================\n\n")
    f.write(f"Following Drexler et al. (Nature Medicine 2024) approach:\n")
    f.write(f"'The combined dataset was dichotomized into low- and high-neural\n")
    f.write(f"tumors using the median neural proportion.'\n\n")
    f.write(f"Dataset: {n_total} GBM samples\n\n")
    f.write(f"Statistics:\n")
    f.write(f"  Mean: {mean_value:.4f}\n")
    f.write(f"  Median (threshold): {median_threshold:.4f}\n")
    f.write(f"  Std: {std_value:.4f}\n")
    f.write(f"  Min: {min_value:.4f}\n")
    f.write(f"  Max: {max_value:.4f}\n\n")
    f.write(f"Classification rule:\n")
    f.write(f"  High-Neural: Cortical_neurons >= {median_threshold:.4f}\n")
    f.write(f"  Low-Neural:  Cortical_neurons < {median_threshold:.4f}\n")

print(f"\nSaved threshold to: gbm_threshold.txt")

# Classify samples based on median threshold
df['Neural_Classification'] = df['Cortical_neurons'].apply(
    lambda x: 'High-Neural' if x >= median_threshold else 'Low-Neural'
)

# Count classifications
high_neural = (df['Neural_Classification'] == 'High-Neural').sum()
low_neural = (df['Neural_Classification'] == 'Low-Neural').sum()
pct_high = 100 * high_neural / n_total
pct_low = 100 * low_neural / n_total

print(f"\nClassification results:")
print(f"  High-Neural: {high_neural} ({pct_high:.1f}%)")
print(f"  Low-Neural: {low_neural} ({pct_low:.1f}%)")

# Save updated file with classification
df.to_csv(input_file)
print(f"\nUpdated file with Neural_Classification: {input_file}")

# Create the plot
print("\nGenerating plot...")

fig, ax = plt.subplots(figsize=(8, 15))

# Create jittered x positions for beeswarm-like effect
np.random.seed(42)
x_jitter = np.random.normal(0, 0.08, len(df))

# Color points by classification
colors = ["#522404" if c == 'High-Neural' else "#EC9D65" for c in df['Neural_Classification']]

# Plot scatter points with labels for legend
high_mask = (df['Neural_Classification'] == 'High-Neural').values
low_mask = (df['Neural_Classification'] == 'Low-Neural').values

ax.scatter(x_jitter[high_mask], df.loc[high_mask, 'Cortical_neurons'].values,
           c='#522404', alpha=0.7, s=50, edgecolors='none',
           label=f'High-Neural ({pct_high:.1f}%, n={high_neural})')
ax.scatter(x_jitter[low_mask], df.loc[low_mask, 'Cortical_neurons'].values,
           c='#EC9D65', alpha=0.7, s=50, edgecolors='none',
           label=f'Low-Neural ({pct_low:.1f}%, n={low_neural})')

# Add threshold line
ax.axhline(y=median_threshold, color='black', linestyle='-', linewidth=1.5,
           label=f'Threshold ({median_threshold:.2f})')

# Add legend outside the plot
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10, frameon=False)

# Set axis labels and title
ax.set_ylabel('Neural signature', fontsize=14)
ax.set_xlim(-0.5, 0.5)
ax.set_xticks([])

# Set title
ax.set_title(f'200GBMs\n(n={n_total}, threshold={median_threshold:.2f})', fontsize=14, fontweight='bold')

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.tight_layout()

# Save figure
plt.savefig('GBM_Moss_neural_signature_plot.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('GBM_Moss_neural_signature_plot.pdf', bbox_inches='tight', facecolor='white')

print(f"\nSaved plots:")
print(f"  GBM_Moss_neural_signature_plot.png")
print(f"  GBM_Moss_neural_signature_plot.pdf")

plt.close()

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
