#!/usr/bin/env python3
"""
Neural classification based on Cortical_neurons threshold of 0.5.

Reads Moss deconvolution results and classifies samples as:
- High-Neural: Cortical_neurons >= 0.5
- Low-Neural: Cortical_neurons < 0.5

Usage:
    python threshold_05_classification.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Configuration
INPUT_FILE = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/GBM_Moss_signatures.csv'
THRESHOLD = 0.23
OUTPUT_FILE = 'GBM_Moss_threshold_023_classification.csv'

# Read data
print("=" * 60)
print("Neural Classification (Threshold = 0.5)")
print("=" * 60)

df = pd.read_csv(INPUT_FILE)
print(f"Loaded {len(df)} samples from {INPUT_FILE}")

# Apply classification with threshold 0.5
df['Neural_Classification_05'] = df['Cortical_neurons'].apply(
    lambda x: 'High-Neural' if x >= THRESHOLD else 'Low-Neural'
)

# Create output dataframe
output_df = df[['Sample_id', 'Cortical_neurons', 'Neural_Classification_05']].copy()
output_df.columns = ['Sample_id', 'Cortical_neurons', 'Neural_Classification']

# Save results
output_df.to_csv(OUTPUT_FILE, index=False)
print(f"\nSaved classification results to: {OUTPUT_FILE}")

# Count classifications
high_count = (output_df['Neural_Classification'] == 'High-Neural').sum()
low_count = (output_df['Neural_Classification'] == 'Low-Neural').sum()
total = len(output_df)

print(f"\nClassification Summary (threshold = {THRESHOLD}):")
print(f"  High-Neural: {high_count} ({100*high_count/total:.1f}%)")
print(f"  Low-Neural:  {low_count} ({100*low_count/total:.1f}%)")

# Compare with original classification (threshold 0.229)
if 'Neural_Classification' in df.columns:
    original_high = (df['Neural_Classification'] == 'High-Neural').sum()
    original_low = (df['Neural_Classification'] == 'Low-Neural').sum()
    print(f"\nOriginal Classification (threshold = 0.229):")
    print(f"  High-Neural: {original_high} ({100*original_high/total:.1f}%)")
    print(f"  Low-Neural:  {original_low} ({100*original_low/total:.1f}%)")

# ============================================================================
# VISUALIZATION
# ============================================================================
print("\n" + "=" * 60)
print("Generating Plots")
print("=" * 60)

# Define colors
high_color = '#522404'  # Dark brown for High-Neural
low_color = '#EC9D65'   # Light orange for Low-Neural

# Create figure with 2 subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# --- Plot 1: Pie Chart ---
ax1 = axes[0]
counts = [high_count, low_count]
labels = [f'High-Neural\n(n={high_count})', f'Low-Neural\n(n={low_count})']
colors = [high_color, low_color]

wedges, texts, autotexts = ax1.pie(
    counts,
    labels=labels,
    colors=colors,
    autopct='%1.1f%%',
    startangle=90,
    explode=(0.02, 0.02),
    textprops={'fontsize': 11}
)
for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontweight('bold')
    autotext.set_fontsize(12)

ax1.set_title(f'Neural Classification Distribution\n(Threshold = {THRESHOLD})',
              fontsize=14, fontweight='bold')

# --- Plot 2: Violin/Swarm Plot (like the figure shown) ---
ax2 = axes[1]

# Sort by Cortical_neurons for better visualization
df_sorted = output_df.sort_values('Cortical_neurons')

# Create jittered x positions for swarm-like effect
np.random.seed(42)
n_samples = len(df_sorted)

# Calculate kernel density for width adjustment (violin-like shape)
from scipy import stats

values = df_sorted['Cortical_neurons'].values
kde = stats.gaussian_kde(values)

# Create y positions and x jitter based on density
y_positions = values
density_at_points = kde(values)
# Normalize density to control spread
max_density = density_at_points.max()
normalized_density = density_at_points / max_density

# Jitter x based on density (wider where more points)
x_jitter = np.random.uniform(-1, 1, n_samples) * normalized_density * 0.4

# Color points based on classification
colors_points = [high_color if c == 'High-Neural' else low_color
                 for c in df_sorted['Neural_Classification']]

# Plot points
ax2.scatter(x_jitter, y_positions, c=colors_points, s=25, alpha=0.7, edgecolors='none')

# Add threshold line
ax2.axhline(y=THRESHOLD, color='red', linestyle='-', linewidth=2, label=f'Threshold = {THRESHOLD}')

# Add text labels for regions
ax2.text(0.9, 0.75, 'High-Neural', fontsize=14, fontweight='bold',
         color=high_color, ha='left', va='center', transform=ax2.transAxes,
         rotation=0)
ax2.text(0.9, 0.10, 'Low-Neural', fontsize=14, fontweight='bold',
         color=low_color, ha='left', va='center', transform=ax2.transAxes,
         rotation= 0)

# Formatting
ax2.set_ylabel('Neural Signature\n(Cortical_neurons)', fontsize=12)
ax2.set_xlim(-0.6, 0.8)
ax2.set_ylim(-0.1, max(values) + 0.05)
ax2.set_xticks([])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)

# Add title
ax2.set_title('Neural Signature Distribution', fontsize=14, fontweight='bold')

plt.tight_layout()

# Save figures
plt.savefig('GBM_threshold_023_lassification_plots.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('GBM_threshold_023_lassification_plots.pdf', bbox_inches='tight', facecolor='white')
print("\nSaved plots:")
print("  GBM_threshold_023_lassification_plots.png")
print("  GBM_threshold_023_lassification_plots.pdf")

plt.close()

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
