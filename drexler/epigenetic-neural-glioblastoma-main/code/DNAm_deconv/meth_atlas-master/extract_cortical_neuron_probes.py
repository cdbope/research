#!/usr/bin/env python3
"""
Extract CpG probes that are most informative for Cortical_neurons classification.

This script identifies probes where cortical neurons have distinctly different
methylation levels compared to other cell types, making them most useful for
sparse/low-coverage data.

Three types of informative probes:
1. Neuron-HIGH: High methylation in neurons, low in other tissues
2. Neuron-LOW: Low methylation in neurons, high in other tissues
3. Neuron-VARIABLE: High variance across samples (most discriminative)

Output:
- cortical_neuron_specific_probes.csv: Filtered probe list with specificity scores
- cortical_neuron_reference_atlas.csv: Reduced atlas with only informative probes
"""

import pandas as pd
import numpy as np
from scipy import stats

# File paths
ATLAS_FILE = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/reference_atlas.csv'
OUTPUT_DIR = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/'

# Parameters
MIN_FOLD_CHANGE = 2.0  # Minimum fold change vs other cell types
MIN_DIFF = 0.15        # Minimum absolute difference in methylation
TOP_N_PROBES = 500     # Number of top probes to select

print("Loading reference atlas...")
atlas = pd.read_csv(ATLAS_FILE, index_col=0)

# Remove duplicate CpG IDs (keep first occurrence)
atlas = atlas[~atlas.index.duplicated(keep='first')]

print(f"Atlas shape: {atlas.shape}")
print(f"Cell types: {list(atlas.columns)}")

# Get cortical neuron values
neuron_col = 'Cortical_neurons'
neuron_values = atlas[neuron_col]

# Get other cell types (excluding neurons)
other_cols = [c for c in atlas.columns if c != neuron_col]
other_values = atlas[other_cols]

# Calculate statistics for each probe
print("\nCalculating probe specificity scores...")

results = []
for probe_id in atlas.index:
    neuron_val = neuron_values[probe_id]
    other_vals = other_values.loc[probe_id].values

    # Statistics
    other_mean = np.mean(other_vals)
    other_std = np.std(other_vals)
    other_min = np.min(other_vals)
    other_max = np.max(other_vals)

    # Difference metrics
    diff_vs_mean = neuron_val - other_mean
    abs_diff = abs(diff_vs_mean)

    # Fold change (avoid division by zero)
    if other_mean > 0.01:
        fold_change = neuron_val / other_mean
    else:
        fold_change = neuron_val / 0.01

    # Z-score: how many SDs is neuron value from other cell types
    if other_std > 0.001:
        z_score = (neuron_val - other_mean) / other_std
    else:
        z_score = 0

    # Specificity score: combination of metrics
    # High score = more specific to neurons
    specificity = abs(z_score) * abs_diff

    # Determine probe type
    if neuron_val > other_mean + MIN_DIFF and neuron_val > other_max:
        probe_type = 'Neuron-HIGH'
    elif neuron_val < other_mean - MIN_DIFF and neuron_val < other_min:
        probe_type = 'Neuron-LOW'
    elif abs_diff > MIN_DIFF:
        probe_type = 'Neuron-DIFF'
    else:
        probe_type = 'Non-specific'

    results.append({
        'CpG': probe_id,
        'Neuron_beta': neuron_val,
        'Other_mean': other_mean,
        'Other_std': other_std,
        'Other_min': other_min,
        'Other_max': other_max,
        'Diff_vs_mean': diff_vs_mean,
        'Abs_diff': abs_diff,
        'Fold_change': fold_change,
        'Z_score': z_score,
        'Specificity_score': specificity,
        'Probe_type': probe_type
    })

results_df = pd.DataFrame(results)

# Summary statistics
print("\n" + "=" * 60)
print("PROBE TYPE SUMMARY")
print("=" * 60)
for ptype in results_df['Probe_type'].unique():
    count = (results_df['Probe_type'] == ptype).sum()
    print(f"  {ptype}: {count}")

# Filter to specific probes
specific_probes = results_df[results_df['Probe_type'] != 'Non-specific'].copy()
print(f"\nProbes with neuron-specific signal: {len(specific_probes)}")

# Sort by specificity score and select top N
specific_probes = specific_probes.sort_values('Specificity_score', ascending=False)
top_probes = specific_probes.head(TOP_N_PROBES)

print(f"Selected top {len(top_probes)} most informative probes")

# Print top 20 probes
print("\n" + "=" * 60)
print("TOP 20 MOST INFORMATIVE PROBES FOR CORTICAL NEURONS")
print("=" * 60)
print(f"{'CpG':<12} {'Neuron':<8} {'Others':<8} {'Diff':<8} {'Z-score':<8} {'Type'}")
print("-" * 60)
for _, row in top_probes.head(20).iterrows():
    print(f"{row['CpG']:<12} {row['Neuron_beta']:.3f}    {row['Other_mean']:.3f}    {row['Diff_vs_mean']:+.3f}   {row['Z_score']:+.2f}     {row['Probe_type']}")

# Save probe list
output_probes = OUTPUT_DIR + 'cortical_neuron_specific_probes.csv'
top_probes.to_csv(output_probes, index=False)
print(f"\nSaved: {output_probes}")

# Create reduced reference atlas with only top probes
top_cpgs = top_probes['CpG'].tolist()
reduced_atlas = atlas.loc[atlas.index.isin(top_cpgs)]
reduced_atlas = reduced_atlas.loc[top_cpgs]  # Maintain order

output_atlas = OUTPUT_DIR + 'cortical_neuron_reference_atlas.csv'
reduced_atlas.to_csv(output_atlas)
print(f"Saved: {output_atlas}")

# Also save just the CpG IDs for use with generate_beta scripts
output_cpg_list = OUTPUT_DIR + 'cortical_neuron_cpg_list.txt'
with open(output_cpg_list, 'w') as f:
    for cpg in top_cpgs:
        f.write(cpg + '\n')
print(f"Saved: {output_cpg_list}")

# Statistics on selected probes
print("\n" + "=" * 60)
print("SELECTED PROBES STATISTICS")
print("=" * 60)
print(f"Total probes: {len(top_probes)}")
print(f"  Neuron-HIGH: {(top_probes['Probe_type'] == 'Neuron-HIGH').sum()}")
print(f"  Neuron-LOW: {(top_probes['Probe_type'] == 'Neuron-LOW').sum()}")
print(f"  Neuron-DIFF: {(top_probes['Probe_type'] == 'Neuron-DIFF').sum()}")

print(f"\nNeuron beta value distribution:")
print(f"  Mean: {top_probes['Neuron_beta'].mean():.3f}")
print(f"  Std:  {top_probes['Neuron_beta'].std():.3f}")
print(f"  Min:  {top_probes['Neuron_beta'].min():.3f}")
print(f"  Max:  {top_probes['Neuron_beta'].max():.3f}")

print(f"\nSpecificity score distribution:")
print(f"  Mean: {top_probes['Specificity_score'].mean():.3f}")
print(f"  Min:  {top_probes['Specificity_score'].min():.3f}")
print(f"  Max:  {top_probes['Specificity_score'].max():.3f}")

print("\nDone!")
