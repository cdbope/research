#!/usr/bin/env python3
"""
Mutation Landscape Analysis Script
Analyzes variant calling data from Clair3 and ClairS-TO
Processes multiple samples and aggregates results
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from glob import glob

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

def load_all_variant_data(data_dir):
    """Load and aggregate variant calling data from all samples"""
    # Find all sample files
    file_pattern = f"{data_dir}/*_merge_annotation_filter_snvs_allcall.csv"
    files = sorted(glob(file_pattern))

    print(f"Found {len(files)} sample files")

    all_data = []
    sample_counts = {}

    for file_path in files:
        # Extract sample ID from filename
        sample_id = Path(file_path).name.replace('_merge_annotation_filter_snvs_allcall.csv', '')

        try:
            df = pd.read_csv(file_path, sep='\t')
            df['Sample_ID'] = sample_id
            all_data.append(df)
            sample_counts[sample_id] = len(df)
            print(f"  Loaded {sample_id}: {len(df)} variants")
        except Exception as e:
            print(f"  Error loading {sample_id}: {e}")

    # Combine all data
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        print(f"\nTotal variants across all samples: {len(combined_df)}")

        # Clean data - convert AF and Depth to numeric, handling any non-numeric values
        combined_df['AF'] = pd.to_numeric(combined_df['AF'], errors='coerce')
        combined_df['Depth'] = pd.to_numeric(combined_df['Depth'], errors='coerce')

        # Remove rows with invalid AF or Depth
        initial_len = len(combined_df)
        combined_df = combined_df.dropna(subset=['AF', 'Depth'])
        if len(combined_df) < initial_len:
            print(f"Removed {initial_len - len(combined_df)} variants with invalid AF/Depth values")

        return combined_df, sample_counts
    else:
        raise ValueError("No data files could be loaded!")

def load_variant_data(file_path):
    """Load and preprocess variant calling data from a single file"""
    print(f"Loading data from: {file_path}")
    df = pd.read_csv(file_path, sep='\t')
    print(f"Loaded {len(df)} variants")
    return df

def plot_gene_mutation_frequency(df, top_n=20, output_dir='.'):
    """Plot the frequency of mutations per gene"""
    gene_counts = df['Gene.refGene'].value_counts().head(top_n)
    total_variants = gene_counts.sum()

    plt.figure(figsize=(12, 8))
    bars = plt.barh(range(len(gene_counts)), gene_counts.values, color='steelblue', edgecolor='black')
    plt.yticks(range(len(gene_counts)), gene_counts.index)
    plt.xlabel('Number of Variants', fontsize=12)
    plt.ylabel('Gene', fontsize=12)
    plt.title('Top Mutated Genes', fontsize=14, fontweight='bold')

    # Add percentage labels at the end of each bar
    for i, (gene, count) in enumerate(gene_counts.items()):
        percentage = (count / total_variants) * 100
        plt.text(count, i, f' {percentage:.1f}%',
                va='center', fontsize=9, fontweight='bold')

    plt.tight_layout()

    output_path = Path(output_dir) / 'gene_mutation_frequency.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

    return gene_counts

def plot_mutation_types(df, output_dir='.'):
    """Plot distribution of mutation functional types"""
    fig = plt.figure(figsize=(16, 10))

    # Functional region - top plot
    ax1 = plt.subplot(2, 2, (1, 2))
    func_counts = df['Func.refGene'].value_counts()
    colors1 = plt.cm.Set3(range(len(func_counts)))
    wedges1, texts1, autotexts1 = ax1.pie(func_counts.values, autopct='%1.1f%%',
                                            startangle=90, colors=colors1,
                                            textprops={'fontsize': 10, 'fontweight': 'bold'})
    ax1.set_title('Functional Region Distribution', fontsize=14, fontweight='bold', pad=20)

    # Add legend for functional region
    ax1.legend(wedges1, [f'{label} (n={count})' for label, count in func_counts.items()],
              title="Functional Regions",
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize=10)

    # Exonic function - use bar chart instead of pie for better readability
    exonic_df = df[df['Func.refGene'] == 'exonic']
    if len(exonic_df) > 0:
        ax2 = plt.subplot(2, 1, 2)
        exonic_counts = exonic_df['ExonicFunc.refGene'].value_counts()
        total_exonic = exonic_counts.sum()

        # Create horizontal bar chart
        bars = ax2.barh(range(len(exonic_counts)), exonic_counts.values,
                       color=plt.cm.Set2(range(len(exonic_counts))),
                       edgecolor='black', alpha=0.8)
        ax2.set_yticks(range(len(exonic_counts)))
        ax2.set_yticklabels(exonic_counts.index, fontsize=10)
        ax2.set_xlabel('Number of Variants', fontsize=12, fontweight='bold')
        ax2.set_title('Exonic Function Distribution', fontsize=14, fontweight='bold', pad=15)
        ax2.grid(axis='x', alpha=0.3)

        # Add count and percentage labels
        for i, (bar, count) in enumerate(zip(bars, exonic_counts.values)):
            percentage = (count / total_exonic) * 100
            ax2.text(count, i, f'  {count} ({percentage:.1f}%)',
                    va='center', fontsize=10, fontweight='bold')

    plt.tight_layout()
    output_path = Path(output_dir) / 'mutation_types.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_allele_frequency_distribution(df, output_dir='.'):
    """Plot variant allele frequency distribution"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Histogram
    ax1.hist(df['AF'], bins=30, color='coral', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Allele Frequency (AF)', fontsize=12)
    ax1.set_ylabel('Number of Variants', fontsize=12)
    ax1.set_title('Allele Frequency Distribution', fontsize=12, fontweight='bold')
    ax1.axvline(df['AF'].median(), color='red', linestyle='--',
                label=f'Median: {df["AF"].median():.3f}')
    ax1.legend()

    # Box plot by gene (top genes)
    top_genes = df['Gene.refGene'].value_counts().head(10).index
    df_top = df[df['Gene.refGene'].isin(top_genes)]
    df_top.boxplot(column='AF', by='Gene.refGene', ax=ax2)
    ax2.set_xlabel('Gene', fontsize=12)
    ax2.set_ylabel('Allele Frequency (AF)', fontsize=12)
    ax2.set_title('AF Distribution by Top Genes', fontsize=12, fontweight='bold')
    plt.suptitle('')  # Remove automatic title
    ax2.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    output_path = Path(output_dir) / 'allele_frequency_distribution.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_clinical_significance(df, output_dir='.'):
    """Plot clinical significance distribution"""
    # Filter out 0 values
    df_clnsig = df[df['CLNSIG'] != '0'].copy()

    if len(df_clnsig) > 0:
        plt.figure(figsize=(10, 6))
        clnsig_counts = df_clnsig['CLNSIG'].value_counts()
        total_clnsig = clnsig_counts.sum()

        ax = clnsig_counts.plot(kind='bar', color='teal', edgecolor='black')
        plt.xlabel('Clinical Significance', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.title('Clinical Significance Distribution', fontsize=14, fontweight='bold')
        plt.xticks(rotation=45, ha='right')

        # Add percentage labels on top of each bar
        for i, (category, count) in enumerate(clnsig_counts.items()):
            percentage = (count / total_clnsig) * 100
            plt.text(i, count, f'{percentage:.1f}%',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

        plt.tight_layout()

        output_path = Path(output_dir) / 'clinical_significance.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_path}")
        plt.close()
    else:
        print("No clinical significance data available")

def plot_variant_caller_comparison(df, output_dir='.'):
    """Compare variants called by different callers"""
    plt.figure(figsize=(10, 6))

    # Count variants by caller
    caller_counts = df['Variant_caller'].value_counts()
    ax = caller_counts.plot(kind='bar', color='mediumpurple', edgecolor='black')
    plt.xlabel('Variant Caller', fontsize=12)
    plt.ylabel('Number of Variants', fontsize=12)
    plt.title('Variants by Caller', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    output_path = Path(output_dir) / 'variant_caller_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_chromosome_distribution(df, output_dir='.'):
    """Plot variant distribution across chromosomes"""
    plt.figure(figsize=(14, 6))

    chr_counts = df['Chr'].value_counts().sort_index()
    # Sort chromosomes properly (chr1, chr2, ..., chr22, chrX, chrY)
    chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    chr_counts = chr_counts.reindex([c for c in chr_order if c in chr_counts.index])

    ax = chr_counts.plot(kind='bar', color='skyblue', edgecolor='black')
    plt.xlabel('Chromosome', fontsize=12)
    plt.ylabel('Number of Variants', fontsize=12)
    plt.title('Variant Distribution Across Chromosomes', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45)
    plt.tight_layout()

    output_path = Path(output_dir) / 'chromosome_distribution.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_depth_distribution(df, output_dir='.'):
    """Plot sequencing depth distribution"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Overall depth distribution
    ax1.hist(df['Depth'], bins=30, color='lightgreen', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Sequencing Depth', fontsize=12)
    ax1.set_ylabel('Number of Variants', fontsize=12)
    ax1.set_title('Sequencing Depth Distribution', fontsize=12, fontweight='bold')
    ax1.axvline(df['Depth'].median(), color='red', linestyle='--',
                label=f'Median: {df["Depth"].median():.0f}')
    ax1.legend()

    # Depth vs AF scatter
    ax2.scatter(df['Depth'], df['AF'], alpha=0.5, color='darkgreen', s=50)
    ax2.set_xlabel('Sequencing Depth', fontsize=12)
    ax2.set_ylabel('Allele Frequency', fontsize=12)
    ax2.set_title('Depth vs Allele Frequency', fontsize=12, fontweight='bold')

    plt.tight_layout()
    output_path = Path(output_dir) / 'depth_distribution.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_sample_mutation_burden(sample_counts, output_dir='.'):
    """Plot mutation burden across samples"""
    plt.figure(figsize=(16, 6))

    samples = list(sample_counts.keys())
    counts = list(sample_counts.values())

    # Sort by count
    sorted_data = sorted(zip(samples, counts), key=lambda x: x[1], reverse=True)
    samples, counts = zip(*sorted_data)

    plt.bar(range(len(samples)), counts, color='indianred', edgecolor='black', alpha=0.7)
    plt.xlabel('Sample', fontsize=12)
    plt.ylabel('Number of Variants', fontsize=12)
    plt.title('Mutation Burden Across Samples', fontsize=14, fontweight='bold')
    plt.axhline(np.mean(counts), color='blue', linestyle='--',
                label=f'Mean: {np.mean(counts):.1f}')
    plt.axhline(np.median(counts), color='green', linestyle='--',
                label=f'Median: {np.median(counts):.1f}')
    plt.legend()
    plt.tight_layout()

    output_path = Path(output_dir) / 'sample_mutation_burden.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_gene_recurrence_heatmap(df, top_genes=20, output_dir='.'):
    """Create a heatmap showing gene mutation recurrence across samples"""
    # Get top mutated genes
    top_gene_list = df['Gene.refGene'].value_counts().head(top_genes).index

    # Create a matrix: samples x genes
    samples = sorted(df['Sample_ID'].unique())
    matrix = pd.DataFrame(0, index=samples, columns=top_gene_list)

    for sample in samples:
        sample_df = df[df['Sample_ID'] == sample]
        gene_counts = sample_df['Gene.refGene'].value_counts()
        for gene in top_gene_list:
            if gene in gene_counts:
                matrix.loc[sample, gene] = gene_counts[gene]

    # Create heatmap
    plt.figure(figsize=(20, 12))
    sns.heatmap(matrix.T, cmap='YlOrRd', cbar_kws={'label': 'Mutation Count'},
                linewidths=0.5, linecolor='gray')
    plt.xlabel('Sample', fontsize=12)
    plt.ylabel('Gene', fontsize=12)
    plt.title(f'Mutation Recurrence Heatmap (Top {top_genes} Genes)', fontsize=14, fontweight='bold')
    plt.xticks(rotation=90, fontsize=6)
    plt.yticks(fontsize=8)
    plt.tight_layout()

    output_path = Path(output_dir) / 'gene_recurrence_heatmap.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

    # Save matrix to CSV
    csv_path = Path(output_dir) / 'gene_recurrence_matrix.csv'
    matrix.to_csv(csv_path)
    print(f"Saved: {csv_path}")

def plot_gene_frequency_across_samples(df, top_n=20, output_dir='.'):
    """Plot how frequently each gene is mutated across different samples"""
    top_genes = df['Gene.refGene'].value_counts().head(top_n).index

    gene_sample_freq = {}
    total_samples = df['Sample_ID'].nunique()

    for gene in top_genes:
        # Count how many samples have this gene mutated
        samples_with_gene = df[df['Gene.refGene'] == gene]['Sample_ID'].nunique()
        gene_sample_freq[gene] = samples_with_gene

    # Sort by frequency
    sorted_genes = sorted(gene_sample_freq.items(), key=lambda x: x[1], reverse=True)
    genes, frequencies = zip(*sorted_genes)

    plt.figure(figsize=(12, 8))
    bars = plt.barh(range(len(genes)), frequencies, color='steelblue', edgecolor='black')
    plt.yticks(range(len(genes)), genes)
    plt.xlabel('Number of Samples with Mutation', fontsize=12)
    plt.ylabel('Gene', fontsize=12)
    plt.title('Gene Mutation Frequency Across Samples', fontsize=14, fontweight='bold')

    # Add percentage labels at the end of each bar
    for i, (bar, freq) in enumerate(zip(bars, frequencies)):
        percentage = (freq / total_samples) * 100
        plt.text(freq, bar.get_y() + bar.get_height()/2,
                f' {percentage:.1f}%',
                va='center', fontsize=9, fontweight='bold')

    plt.tight_layout()

    output_path = Path(output_dir) / 'gene_sample_frequency.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

    return gene_sample_freq

def plot_target_gene_coverage(df, bed_file, output_dir='.'):
    """Plot coverage of target genes from the OCC protein coding bed file"""
    # Parse bed file to extract gene names
    target_genes = []
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            # Extract gene_name from the attributes field (10th column, index 9)
            if len(fields) >= 10:
                attributes = fields[9]
                for attr in attributes.split(';'):
                    if 'gene_name=' in attr:
                        gene_name = attr.split('gene_name=')[1]
                        target_genes.append(gene_name)
                        break

    target_genes = sorted(set(target_genes))  # Remove duplicates
    print(f"  Found {len(target_genes)} unique target genes in bed file")

    # Get genes found in variant data
    genes_in_data = set(df['Gene.refGene'].unique())
    print(f"  Found {len(genes_in_data)} unique genes in variant data")

    # Calculate coverage
    covered_genes = [g for g in target_genes if g in genes_in_data]
    uncovered_genes = [g for g in target_genes if g not in genes_in_data]

    coverage_pct = (len(covered_genes) / len(target_genes)) * 100

    # Create visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # 1. Coverage pie chart
    coverage_data = [len(covered_genes), len(uncovered_genes)]
    coverage_labels = [f'Covered\n({len(covered_genes)} genes)',
                      f'Not Covered\n({len(uncovered_genes)} genes)']
    colors = ['#2ecc71', '#e74c3c']

    wedges, texts, autotexts = ax1.pie(coverage_data, labels=coverage_labels,
                                         autopct='%1.1f%%', colors=colors,
                                         startangle=90, textprops={'fontsize': 12, 'fontweight': 'bold'})
    ax1.set_title(f'Target Gene Coverage\n({len(target_genes)} target genes total)',
                  fontsize=14, fontweight='bold', pad=20)

    # 2. Top covered target genes by mutation count
    covered_gene_counts = df[df['Gene.refGene'].isin(covered_genes)]['Gene.refGene'].value_counts().head(30)

    if len(covered_gene_counts) > 0:
        bars = ax2.barh(range(len(covered_gene_counts)), covered_gene_counts.values,
                       color='#3498db', edgecolor='black', alpha=0.8)
        ax2.set_yticks(range(len(covered_gene_counts)))
        ax2.set_yticklabels(covered_gene_counts.index, fontsize=10)
        ax2.set_xlabel('Number of Mutations', fontsize=12, fontweight='bold')
        ax2.set_title(f'Top 20 Mutated Target Genes\n({coverage_pct:.1f}% of target genes have mutations)',
                     fontsize=14, fontweight='bold', pad=15)
        ax2.grid(axis='x', alpha=0.3)

        # Add count labels
        for i, count in enumerate(covered_gene_counts.values):
            ax2.text(count, i, f' {count}', va='center', fontsize=9, fontweight='bold')

    plt.tight_layout()
    output_path = Path(output_dir) / 'target_gene_coverage.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

    # Save detailed coverage report
    report_path = Path(output_dir) / 'target_gene_coverage_report.txt'
    with open(report_path, 'w') as f:
        f.write("=== TARGET GENE COVERAGE REPORT ===\n\n")
        f.write(f"Total target genes: {len(target_genes)}\n")
        f.write(f"Covered genes (with mutations): {len(covered_genes)} ({coverage_pct:.1f}%)\n")
        f.write(f"Uncovered genes (no mutations): {len(uncovered_genes)} ({100-coverage_pct:.1f}%)\n\n")

        f.write("=== TOP 20 MUTATED TARGET GENES ===\n")
        for gene, count in covered_gene_counts.items():
            samples_with_gene = df[df['Gene.refGene'] == gene]['Sample_ID'].nunique()
            f.write(f"{gene}: {count} mutations in {samples_with_gene} samples\n")

        f.write(f"\n=== UNCOVERED TARGET GENES ({len(uncovered_genes)} genes) ===\n")
        for gene in sorted(uncovered_genes):
            f.write(f"{gene}\n")

    print(f"Saved: {report_path}")

    return {
        'total_target_genes': len(target_genes),
        'covered_genes': len(covered_genes),
        'uncovered_genes': len(uncovered_genes),
        'coverage_percentage': coverage_pct
    }

def generate_summary_stats(df, sample_counts=None, output_dir='.'):
    """Generate summary statistics table"""
    summary = {
        'Total Variants': len(df),
        'Total Samples': df['Sample_ID'].nunique() if 'Sample_ID' in df.columns else 1,
        'Unique Genes': df['Gene.refGene'].nunique(),
        'Chromosomes': df['Chr'].nunique(),
        'Mean Depth': df['Depth'].mean(),
        'Median Depth': df['Depth'].median(),
        'Mean AF': df['AF'].mean(),
        'Median AF': df['AF'].median(),
        'Exonic Variants': (df['Func.refGene'] == 'exonic').sum(),
        'Upstream Variants': (df['Func.refGene'] == 'upstream').sum(),
        'Variants with COSMIC ID': (df['COSMIC100'] != '0').sum(),
        'Clinically Significant': (df['CLNSIG'] != '0').sum(),
    }

    # Save to file
    output_path = Path(output_dir) / 'summary_statistics.txt'
    with open(output_path, 'w') as f:
        f.write("=== MUTATION LANDSCAPE SUMMARY STATISTICS ===\n\n")
        for key, value in summary.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.2f}\n")
            else:
                f.write(f"{key}: {value}\n")

        if sample_counts:
            f.write(f"\n=== SAMPLE STATISTICS ===\n")
            f.write(f"Mean variants per sample: {np.mean(list(sample_counts.values())):.2f}\n")
            f.write(f"Median variants per sample: {np.median(list(sample_counts.values())):.2f}\n")
            f.write(f"Min variants per sample: {min(sample_counts.values())}\n")
            f.write(f"Max variants per sample: {max(sample_counts.values())}\n")

        f.write("\n=== TOP 20 MUTATED GENES ===\n")
        top_genes = df['Gene.refGene'].value_counts().head(20)
        for gene, count in top_genes.items():
            if 'Sample_ID' in df.columns:
                samples_with_gene = df[df['Gene.refGene'] == gene]['Sample_ID'].nunique()
                f.write(f"{gene}: {count} mutations in {samples_with_gene} samples\n")
            else:
                f.write(f"{gene}: {count}\n")

    print(f"Saved: {output_path}")

    # Print to console
    print("\n=== SUMMARY STATISTICS ===")
    for key, value in summary.items():
        if isinstance(value, float):
            print(f"{key}: {value:.2f}")
        else:
            print(f"{key}: {value}")

def main():
    # Input directory containing all sample files
    data_dir = "/home/chbope/extension/data/200GMBs/results/merge_annot_clair3andclairsto_update"

    # Output directory
    output_dir = Path("/home/chbope/extension/script/mut_landscape/results")
    output_dir.mkdir(exist_ok=True)

    print("="*80)
    print("MUTATION LANDSCAPE ANALYSIS - AGGREGATED ACROSS ALL SAMPLES")
    print("="*80)

    # Load all data
    print("\nLoading variant data from all samples...")
    df, sample_counts = load_all_variant_data(data_dir)

    # Save combined data to CSV for future reference
    combined_csv = output_dir / 'combined_variants_all_samples.csv'
    df.to_csv(combined_csv, index=False)
    print(f"\nSaved combined data to: {combined_csv}")

    # Generate all plots
    print("\n" + "="*80)
    print("GENERATING VISUALIZATIONS")
    print("="*80)

    print("\n[1/12] Sample-level analysis...")
    plot_sample_mutation_burden(sample_counts, output_dir=output_dir)

    print("[2/12] Gene recurrence heatmap...")
    plot_gene_recurrence_heatmap(df, top_genes=30, output_dir=output_dir)

    print("[3/12] Gene frequency across samples...")
    plot_gene_frequency_across_samples(df, top_n=20, output_dir=output_dir)

    print("[4/12] Overall gene mutation frequency...")
    plot_gene_mutation_frequency(df, top_n=30, output_dir=output_dir)

    print("[5/12] Mutation types...")
    plot_mutation_types(df, output_dir=output_dir)

    print("[6/12] Allele frequency distribution...")
    plot_allele_frequency_distribution(df, output_dir=output_dir)

    print("[7/12] Clinical significance...")
    plot_clinical_significance(df, output_dir=output_dir)

    print("[8/12] Variant caller comparison...")
    plot_variant_caller_comparison(df, output_dir=output_dir)

    print("[9/12] Chromosome distribution...")
    plot_chromosome_distribution(df, output_dir=output_dir)

    print("[10/12] Depth distribution...")
    plot_depth_distribution(df, output_dir=output_dir)

    print("[11/12] Target gene coverage...")
    bed_file = "/home/chbope/extension/nWGS_manuscript_data/data/reference/OCC.protein_coding.bed"
    target_coverage = plot_target_gene_coverage(df, bed_file, output_dir=output_dir)

    # Generate summary statistics
    print("[12/12] Summary statistics...")
    generate_summary_stats(df, sample_counts=sample_counts, output_dir=output_dir)

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nResults saved to: {output_dir}")
    print(f"Total samples analyzed: {len(sample_counts)}")
    print(f"Total variants: {len(df)}")
    print(f"Unique genes: {df['Gene.refGene'].nunique()}")
    print(f"\nTarget Gene Coverage:")
    print(f"  - Total target genes: {target_coverage['total_target_genes']}")
    print(f"  - Covered (with mutations): {target_coverage['covered_genes']} ({target_coverage['coverage_percentage']:.1f}%)")
    print(f"  - Uncovered (no mutations): {target_coverage['uncovered_genes']}")
    print("\nGenerated files:")
    print("  - Combined variant data (CSV)")
    print("  - Gene recurrence matrix (CSV)")
    print("  - Target gene coverage report (TXT)")
    print("  - 11 visualization plots (PNG)")
    print("  - Summary statistics (TXT)")
    print("="*80)

if __name__ == "__main__":
    main()
