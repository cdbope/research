#!/usr/bin/env python3
"""
Merge Fusion Events Script
Combines all fusion event files from multiple samples into a single file
Handles files that contain only headers (no data rows)
Creates publication-ready visualizations
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from glob import glob
from collections import Counter

# Set style for publication-quality figures
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

def merge_fusion_events(input_dir, output_file):
    """
    Merge all fusion event TSV files from input directory

    Args:
        input_dir: Directory containing fusion event TSV files
        output_file: Path to output merged file
    """
    # Find all fusion event files
    file_pattern = f"{input_dir}/*_filter_fusion_event.tsv"
    files = sorted(glob(file_pattern))

    print(f"Found {len(files)} fusion event files")

    if len(files) == 0:
        print("No fusion event files found!")
        return

    # Read first file to get header
    first_file = files[0]
    try:
        header_df = pd.read_csv(first_file, sep='\t', nrows=0)
        header = list(header_df.columns)
        print(f"Header columns: {header}")
    except Exception as e:
        print(f"Error reading header from {first_file}: {e}")
        return

    # Collect all data
    all_data = []
    files_with_data = 0
    files_header_only = 0

    for file_path in files:
        # Extract sample ID from filename
        sample_id = Path(file_path).name.replace('_filter_fusion_event.tsv', '')

        try:
            # Read file (try different encodings)
            try:
                df = pd.read_csv(file_path, sep='\t')
            except UnicodeDecodeError:
                # Try latin-1 encoding
                df = pd.read_csv(file_path, sep='\t', encoding='latin-1')

            # Check if file has data (more than just header)
            if len(df) > 0:
                # Add sample ID column
                df.insert(0, 'Sample_ID', sample_id)
                all_data.append(df)
                files_with_data += 1
                print(f"  {sample_id}: {len(df)} fusion events")
            else:
                files_header_only += 1
                print(f"  {sample_id}: header only (no data)")

        except Exception as e:
            print(f"  Error reading {sample_id}: {e}")

    # Summary
    print(f"\n{'='*60}")
    print(f"Files processed: {len(files)}")
    print(f"Files with data: {files_with_data}")
    print(f"Files with header only: {files_header_only}")
    print(f"{'='*60}")

    # Combine all data
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        print(f"\nTotal fusion events: {len(combined_df)}")
        print(f"Samples with fusion events: {combined_df['Sample_ID'].nunique()}")

        # Save combined data
        combined_df.to_csv(output_file, sep='\t', index=False)
        print(f"\nMerged data saved to: {output_file}")

        # Show summary statistics
        print(f"\n{'='*60}")
        print("FUSION EVENT SUMMARY")
        print(f"{'='*60}")

        # Count by sample
        sample_counts = combined_df['Sample_ID'].value_counts()
        print(f"\nTop 10 samples with most fusion events:")
        for sample, count in sample_counts.head(10).items():
            print(f"  {sample}: {count}")

        # Count by SV type if available
        if 'svtype' in combined_df.columns:
            print(f"\nFusion events by SV type:")
            svtype_counts = combined_df['svtype'].value_counts()
            for svtype, count in svtype_counts.items():
                print(f"  {svtype}: {count}")

        # Most common genes if available
        if 'Genes' in combined_df.columns:
            print(f"\nTop 10 genes involved in fusion events:")
            # Split genes and count occurrences
            all_genes = []
            for genes in combined_df['Genes'].dropna():
                if isinstance(genes, str):
                    # Handle different gene separators
                    gene_list = genes.replace('--', ',').replace('-', ',').split(',')
                    all_genes.extend([g.strip() for g in gene_list if g.strip()])

            if all_genes:
                from collections import Counter
                gene_counts = Counter(all_genes)
                for gene, count in gene_counts.most_common(10):
                    print(f"  {gene}: {count}")

        print(f"{'='*60}")

        return combined_df

    else:
        print("\nNo fusion events found in any file (all files contain only headers)")
        return None

def create_fusion_visualizations(df, output_dir):
    """Create publication-ready visualizations for fusion events"""

    print(f"\n{'='*60}")
    print("GENERATING VISUALIZATIONS")
    print(f"{'='*60}\n")

    # 1. Fusion prevalence across samples
    print("[1/5] Creating fusion prevalence plot...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Pie chart
    total_samples = 200  # Total samples in cohort
    samples_with_fusions = df['Sample_ID'].nunique()
    samples_without_fusions = total_samples - samples_with_fusions

    colors = ['#e74c3c', '#95a5a6']
    labels = [f'With Fusions\n({samples_with_fusions} samples)',
              f'Without Fusions\n({samples_without_fusions} samples)']

    wedges, texts, autotexts = ax1.pie([samples_with_fusions, samples_without_fusions],
                                         labels=labels, colors=colors, autopct='%1.1f%%',
                                         startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax1.set_title('Fusion Event Prevalence in Cohort\n(n=200 samples)',
                  fontsize=13, fontweight='bold', pad=15)

    # Bar chart - fusion events per sample
    sample_counts = df['Sample_ID'].value_counts().sort_values(ascending=True)
    bars = ax2.barh(range(len(sample_counts)), sample_counts.values,
                    color='#3498db', edgecolor='black', alpha=0.8)
    ax2.set_yticks(range(len(sample_counts)))
    ax2.set_yticklabels(sample_counts.index, fontsize=10)
    ax2.set_xlabel('Number of Fusion Events', fontsize=11, fontweight='bold')
    ax2.set_title('Fusion Events per Sample', fontsize=13, fontweight='bold', pad=15)
    ax2.grid(axis='x', alpha=0.3)

    # Add count labels
    for i, count in enumerate(sample_counts.values):
        ax2.text(count, i, f' {count}', va='center', fontsize=10, fontweight='bold')

    plt.tight_layout()
    output_path = Path(output_dir) / 'fusion_prevalence.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"   Saved: {output_path}")
    plt.close()

    # 2. SV type distribution
    print("[2/5] Creating SV type distribution plot...")
    if 'svtype' in df.columns:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        svtype_counts = df['svtype'].value_counts()

        # Pie chart
        colors_sv = plt.cm.Set3(range(len(svtype_counts)))
        wedges, texts, autotexts = ax1.pie(svtype_counts.values, labels=svtype_counts.index,
                                             colors=colors_sv, autopct='%1.1f%%',
                                             startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
        ax1.set_title('Distribution by SV Type', fontsize=13, fontweight='bold', pad=15)

        # Bar chart
        bars = ax2.bar(range(len(svtype_counts)), svtype_counts.values,
                       color=colors_sv, edgecolor='black', alpha=0.8)
        ax2.set_xticks(range(len(svtype_counts)))
        ax2.set_xticklabels(svtype_counts.index, fontsize=11, fontweight='bold')
        ax2.set_ylabel('Number of Events', fontsize=11, fontweight='bold')
        ax2.set_title('Fusion Events by SV Type', fontsize=13, fontweight='bold', pad=15)
        ax2.grid(axis='y', alpha=0.3)

        # Add count labels on bars
        for i, (bar, count) in enumerate(zip(bars, svtype_counts.values)):
            percentage = (count / svtype_counts.sum()) * 100
            ax2.text(i, count, f'{count}\n({percentage:.1f}%)',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

        plt.tight_layout()
        output_path = Path(output_dir) / 'fusion_svtype_distribution.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"   Saved: {output_path}")
        plt.close()

    # 3. Gene fusion network/partners
    print("[3/5] Creating gene fusion partners plot...")
    if 'Genes' in df.columns:
        # Extract all genes involved
        all_genes = []
        gene_pairs = []

        for idx, row in df.iterrows():
            genes = row['Genes']
            if isinstance(genes, str) and genes.strip():
                gene = genes.strip()
                all_genes.append(gene)

        # Group by fusion ID to get gene pairs
        fusion_ids = df['ID'].unique()
        for fusion_id in fusion_ids:
            fusion_df = df[df['ID'] == fusion_id]
            if len(fusion_df) == 2:
                genes = fusion_df['Genes'].tolist()
                if all(isinstance(g, str) for g in genes):
                    gene_pairs.append(tuple(sorted(genes)))

        # Count gene occurrences
        gene_counts = Counter(all_genes)

        # Plot top genes
        if gene_counts:
            top_genes = dict(gene_counts.most_common(10))

            fig, ax = plt.subplots(figsize=(10, 7))
            genes = list(top_genes.keys())
            counts = list(top_genes.values())

            bars = ax.barh(range(len(genes)), counts, color='#e67e22',
                          edgecolor='black', alpha=0.8)
            ax.set_yticks(range(len(genes)))
            ax.set_yticklabels(genes, fontsize=11, fontweight='bold')
            ax.set_xlabel('Number of Fusion Events', fontsize=11, fontweight='bold')
            ax.set_title('Genes Involved in Fusion Events', fontsize=13, fontweight='bold', pad=15)
            ax.grid(axis='x', alpha=0.3)

            # Add count labels
            for i, count in enumerate(counts):
                ax.text(count, i, f' {count}', va='center', fontsize=10, fontweight='bold')

            plt.tight_layout()
            output_path = Path(output_dir) / 'fusion_gene_frequency.png'
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"   Saved: {output_path}")
            plt.close()

    # 4. Fusion gene pairs visualization
    print("[4/5] Creating fusion gene pairs plot...")
    if gene_pairs:
        pair_counts = Counter(gene_pairs)

        if pair_counts:
            fig, ax = plt.subplots(figsize=(12, 6))

            pair_labels = [f"{g1}--{g2}" for g1, g2 in pair_counts.keys()]
            pair_values = list(pair_counts.values())

            bars = ax.bar(range(len(pair_labels)), pair_values,
                         color='#9b59b6', edgecolor='black', alpha=0.8)
            ax.set_xticks(range(len(pair_labels)))
            ax.set_xticklabels(pair_labels, rotation=45, ha='right', fontsize=11, fontweight='bold')
            ax.set_ylabel('Number of Occurrences', fontsize=11, fontweight='bold')
            ax.set_title('Recurrent Gene Fusion Pairs', fontsize=13, fontweight='bold', pad=15)
            ax.grid(axis='y', alpha=0.3)

            # Add count labels
            for i, (bar, count) in enumerate(zip(bars, pair_values)):
                ax.text(i, count, f'{count}', ha='center', va='bottom',
                       fontsize=10, fontweight='bold')

            plt.tight_layout()
            output_path = Path(output_dir) / 'fusion_gene_pairs.png'
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"   Saved: {output_path}")
            plt.close()

    # 5. Chromosome distribution
    print("[5/5] Creating chromosome distribution plot...")
    if 'chr' in df.columns:
        chr_counts = df['chr'].value_counts()

        # Sort chromosomes properly
        chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        chr_counts_sorted = chr_counts.reindex([c for c in chr_order if c in chr_counts.index])

        fig, ax = plt.subplots(figsize=(14, 6))
        bars = ax.bar(range(len(chr_counts_sorted)), chr_counts_sorted.values,
                     color='#16a085', edgecolor='black', alpha=0.8)
        ax.set_xticks(range(len(chr_counts_sorted)))
        ax.set_xticklabels([c.replace('chr', '') for c in chr_counts_sorted.index],
                          fontsize=10, fontweight='bold')
        ax.set_xlabel('Chromosome', fontsize=11, fontweight='bold')
        ax.set_ylabel('Number of Fusion Breakpoints', fontsize=11, fontweight='bold')
        ax.set_title('Fusion Event Distribution Across Chromosomes',
                    fontsize=13, fontweight='bold', pad=15)
        ax.grid(axis='y', alpha=0.3)

        # Add count labels
        for i, (bar, count) in enumerate(zip(bars, chr_counts_sorted.values)):
            ax.text(i, count, f'{count}', ha='center', va='bottom',
                   fontsize=9, fontweight='bold')

        plt.tight_layout()
        output_path = Path(output_dir) / 'fusion_chromosome_distribution.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"   Saved: {output_path}")
        plt.close()

    print(f"\n{'='*60}")
    print("All visualizations generated successfully!")
    print(f"{'='*60}")

def main():
    # Input directory
    input_dir = "/home/chbope/extension/data/200GMBs/results/fusion_event"

    # Output file
    output_dir = Path("/home/chbope/extension/script/mut_landscape/results")
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / "merged_fusion_events.tsv"

    print("="*60)
    print("FUSION EVENTS MERGER & VISUALIZATION")
    print("="*60)
    print(f"\nInput directory: {input_dir}")
    print(f"Output file: {output_file}\n")

    # Merge files
    combined_df = merge_fusion_events(input_dir, output_file)

    # Create visualizations if there's data
    if combined_df is not None and len(combined_df) > 0:
        create_fusion_visualizations(combined_df, output_dir)
    else:
        print("\nNo fusion events found - skipping visualization generation")

if __name__ == "__main__":
    main()
