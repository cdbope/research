#!/usr/bin/env python3
"""
Extract the top 50 Cortical_neurons probes sorted by highest methylation value.

Steps:
1. Load reference atlas
2. Sort by Cortical_neurons value (highest to lowest)
3. Select top 50 probes

Input: reference_atlas.csv
Output: top_50_Cortical_neurons_probes_from_reference.csv (50 probes with CpGs and Cortical_neurons values)
"""

import pandas as pd

# File paths
INPUT_FILE = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/cortical_neuron_reference_atlas.csv'
OUTPUT_FILE = '/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/top_50_Cortical_neurons_probes_from_500.csv'

# Load reference atlas
print("Loading reference atlas...")
atlas = pd.read_csv(INPUT_FILE, index_col=0)

print(f"Atlas shape: {atlas.shape}")

# Extract only Cortical_neurons column
cortical_neurons = atlas[['Cortical_neurons']].copy()

# Sort by Cortical_neurons value from highest to lowest
cortical_neurons_sorted = cortical_neurons.sort_values('Cortical_neurons', ascending=False)

# Take the top 50 (highest Cortical_neurons values)
top_50 = cortical_neurons_sorted.head(50)

# Reset index to make CpGs a column
top_50 = top_50.reset_index()
top_50.columns = ['CpGs', 'Cortical_neurons']

# Save to CSV
top_50.to_csv(OUTPUT_FILE, index=False)

print(f"\nExtracted top 50 probes:")
print(f"  Output file: {OUTPUT_FILE}")
print(f"  Probes: {len(top_50)}")

# Print summary
print(f"\nCortical_neurons value statistics:")
print(f"  Mean:   {top_50['Cortical_neurons'].mean():.3f}")
print(f"  Median: {top_50['Cortical_neurons'].median():.3f}")
print(f"  Min:    {top_50['Cortical_neurons'].min():.3f}")
print(f"  Max:    {top_50['Cortical_neurons'].max():.3f}")

# Show first 10 probes (highest Cortical_neurons values)
print(f"\nTop 10 probes (sorted by highest Cortical_neurons value):")
print(top_50.head(10).to_string(index=False))

print("\nDone!")
