#!/usr/bin/env python3
"""
Quick fix script to recalculate gene sample counts correctly.
Uses the existing matrix to count unique samples per gene.
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict

# Paths
MATRIX_FILE = "/home/chbope/extension/script/svmeta/results/matrices/sample_by_sv_matrix.csv"
EVENTS_FILE = "/home/chbope/extension/script/svmeta/results/matrices/sv_events.csv"
GENE_ANNOTATION = "/home/chbope/extension/script/svmeta/external_datasets/refseq_genes_hg38.bed"
GBM_DRIVERS = "/home/chbope/extension/script/svmeta/gbm_driver_genes.txt"
OUTPUT_FILE = "/home/chbope/extension/script/svmeta/results/genes/gene_sv_counts.csv"

print("Loading data...")
print(f"  Matrix: {MATRIX_FILE}")
matrix_df = pd.read_csv(MATRIX_FILE, index_col=0)
print(f"  Shape: {matrix_df.shape}")

print(f"  Events: {EVENTS_FILE}")
events_df = pd.read_csv(EVENTS_FILE)
print(f"  Events: {len(events_df)}")

print(f"  Genes: {GENE_ANNOTATION}")
genes_df = pd.read_csv(
    GENE_ANNOTATION,
    sep="\t",
    names=["chr", "start", "end", "gene_symbol"],
    comment="#",
    dtype={"chr": str}
)
genes_df["start"] = pd.to_numeric(genes_df["start"], errors="coerce")
genes_df["end"] = pd.to_numeric(genes_df["end"], errors="coerce")
genes_df = genes_df.dropna(subset=["start", "end"])

# Normalize chr names
genes_df["chr_norm"] = genes_df["chr"].str.replace("chr", "").str.upper()
events_df["chr_norm"] = events_df["chr"].astype(str).str.replace("chr", "").str.upper()

print(f"  Genes loaded: {len(genes_df)}")

print(f"  Driver genes: {GBM_DRIVERS}")
driver_genes = set()
if os.path.exists(GBM_DRIVERS):
    with open(GBM_DRIVERS) as f:
        driver_genes = {line.strip() for line in f if line.strip()}
print(f"  Driver genes: {len(driver_genes)}")

print("\nAnalyzing gene-SV overlaps...")
gene_sv_counts = defaultdict(lambda: {
    "n_svs": 0,
    "sv_ids": set(),
    "DEL": 0,
    "DUP": 0,
    "INV": 0,
    "BND": 0,
    "INS": 0,
    "gene_id": None,
})

for idx, sv in events_df.iterrows():
    if idx % 10000 == 0:
        print(f"  Processed {idx}/{len(events_df)} SVs...")

    # Find overlapping genes
    overlapping = genes_df[
        (genes_df["chr_norm"] == sv["chr_norm"]) &
        (genes_df["start"] <= sv["end"]) &
        (genes_df["end"] >= sv["pos"])
    ]

    for _, row in overlapping.iterrows():
        symbol = row["gene_symbol"]

        gene_sv_counts[symbol]["n_svs"] += 1
        gene_sv_counts[symbol]["sv_ids"].add(sv["sv_id"])

        # Map SVTYPE
        svtype = sv["svtype"]
        if svtype == "TRA":
            svtype = "BND"

        if svtype in ["DEL", "DUP", "INV", "BND", "INS"]:
            gene_sv_counts[symbol][svtype] += 1

print(f"\nGenes with SVs: {len(gene_sv_counts)}")

print("\nCalculating unique sample counts per gene...")
gene_data = []
for idx, (symbol, data) in enumerate(gene_sv_counts.items()):
    if idx % 1000 == 0:
        print(f"  Processed {idx}/{len(gene_sv_counts)} genes...")

    # Calculate number of unique SAMPLES with SVs in this gene
    sv_ids = list(data["sv_ids"])
    if len(sv_ids) > 0:
        # Get columns for these SVs (if they exist in matrix)
        existing_sv_ids = [sv_id for sv_id in sv_ids if sv_id in matrix_df.columns]
        if existing_sv_ids:
            # Count samples where at least one of these SVs is present
            gene_sample_matrix = matrix_df[existing_sv_ids]
            n_unique_samples = (gene_sample_matrix.sum(axis=1) > 0).sum()
        else:
            n_unique_samples = 0
    else:
        n_unique_samples = 0

    gene_data.append({
        "gene": symbol,
        "gene_id": "",  # Not available in BED format
        "n_svs": data["n_svs"],
        "n_samples": n_unique_samples,  # CORRECTED: unique samples
        "DEL": data["DEL"],
        "DUP": data["DUP"],
        "INV": data["INV"],
        "BND": data["BND"],
        "INS": data["INS"],
        "is_driver": symbol in driver_genes,
    })

print("\nSaving results...")
gene_df = pd.DataFrame(gene_data).sort_values("n_svs", ascending=False)
gene_df.to_csv(OUTPUT_FILE, index=False)

print(f"✓ Saved: {OUTPUT_FILE}")
print(f"  Total genes: {len(gene_df)}")
print(f"  Max samples per gene: {gene_df['n_samples'].max()}")
print(f"  Genes with >200 samples (ERROR): {(gene_df['n_samples'] > 200).sum()}")

print("\nTop 10 genes by n_samples:")
print(gene_df[['gene', 'n_svs', 'n_samples', 'is_driver']].head(10).to_string(index=False))

print("\nDone!")
