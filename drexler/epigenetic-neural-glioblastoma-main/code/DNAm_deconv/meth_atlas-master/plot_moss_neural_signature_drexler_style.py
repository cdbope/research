#!/usr/bin/env python3
"""
Generate Moss Neural Signature Score distribution plot in Drexler et al. style
(Nature Medicine 2024)

This creates a beeswarm-like plot with:
- Vertical text labels "High-Neural" and "Low-Neural" on the right side
- Percentages displayed on the right
- Threshold line at 0.37
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read data
INPUT_FILE = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/combined_deconv_betas_merged_deconv_output_cortical_neuron.csv'
OUTPUT_DIR = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/'

df = pd.read_csv(INPUT_FILE)

# Threshold for classification
THRESHOLD = 0.37

# Calculate statistics
total = len(df)
high_neural = (df['Neural_Classification'] == 'High-Neural').sum()
low_neural = (df['Neural_Classification'] == 'Low-Neural').sum()
pct_high = 100 * high_neural / total
pct_low = 100 * low_neural / total

print(f"Total samples: {total}")
print(f"High-Neural: {high_neural} ({pct_high:.1f}%)")
print(f"Low-Neural: {low_neural} ({pct_low:.1f}%)")
print(f"Threshold: {THRESHOLD}")

# Create figure
fig, ax = plt.subplots(figsize=(2.5, 5))

# Create jittered x positions for beeswarm-like effect
np.random.seed(42)
x_jitter = np.random.normal(0, 0.08, len(df))

# Color points by classification
colors = ["#522404" if c == 'High-Neural' else "#EC9D65" for c in df['Neural_Classification']]

# Plot points
ax.scatter(x_jitter, df['Cortical_neurons'], c=colors, s=25, alpha=0.7, edgecolors='none')

# Add threshold line (black)
ax.axhline(y=THRESHOLD, color='black', linewidth=2, linestyle='-', zorder=0)

# Determine y-axis limits based on data
y_min = df['Cortical_neurons'].min()
y_max = df['Cortical_neurons'].max()
y_range = y_max - y_min
y_pad = y_range * 0.1
y_lim_min = max(0, y_min - y_pad)
y_lim_max = y_max + y_pad

# Calculate positions for labels
high_center = (THRESHOLD + y_lim_max) / 2
low_center = (y_lim_min + THRESHOLD) / 2

# Add labels for High-Neural and Low-Neural regions (vertical, on right side)
ax.text(0.35, high_center, 'High-Neural', fontsize=11, color='#8B0000', fontweight='bold',
        rotation=90, va='center')
ax.text(0.35, low_center, 'Low-Neural', fontsize=11, color='#D2691E', fontweight='bold',
        rotation=90, va='center')

# Add percentage labels
ax.text(0.28, high_center + 0.08, f'{pct_high:.1f}%', fontsize=10, color='#8B0000', fontweight='bold')
ax.text(0.28, low_center + 0.08, f'{pct_low:.1f}%', fontsize=10, color='#D2691E', fontweight='bold')

# Formatting
ax.set_ylabel('Neural signature', fontsize=12)
ax.set_xlim(-0.3, 0.4)
ax.set_ylim(y_lim_min, y_lim_max)

# Set y-ticks at regular intervals
ax.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.set_xticks([])

# Add title
ax.set_title(f'200GBMs (Moss 5mC+5hmC)\n(n={total})', fontsize=11, fontweight='bold')

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.tight_layout()

# Save figure
output_png = OUTPUT_DIR + 'GBM_Moss_merged_neural_signature_plot_drexler_style.png'
output_pdf = OUTPUT_DIR + 'GBM_Moss_merged_neural_signature_plot_drexler_style.pdf'
plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {output_png}")
print(f"Saved: {output_pdf}")

plt.close()
print("\nDone!")
