#!/usr/bin/env python3
"""
Step 3: Build sample × gene matrix with CONSERVATIVE COUNTING

Conservative approach (following npae082.pdf):
"For each gene, we reported only one instance of a specific alteration type,
regardless of how many identical alterations occurred within that same gene."

This means: Max 1 SV per gene per sample (binary: 0 or 1)

Input: Merged VCF from step 02
Output: Gene-level matrix with binary counts
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input (read-only from svmeta - we don't modify these)
MERGED_VCF = "/home/chbope/extension/script/svmeta/results/merged/merged_SV.vcf.gz"
GENE_ANNOTATION_BED = "/home/chbope/extension/script/svmeta/external_datasets/refseq_genes_hg38.bed"
GBM_DRIVERS = "/home/chbope/extension/script/svmeta/gbm_driver_genes.txt"

# Output (all outputs go to svmetaunique)
OUTPUT_DIR = "/home/chbope/extension/script/svmetaunique/results"

# ============================================================================

def load_gene_annotations(gene_bed_file):
    """
    Load gene annotations from BED file
    Expected format: chr start end gene_name
    """
    print(f"Loading gene annotations from: {gene_bed_file}")

    if not os.path.exists(gene_bed_file):
        print(f"⚠️  Gene BED file not found: {gene_bed_file}")
        print("Creating from GFF3...")
        return create_gene_bed_from_gff3()

    genes = []
    with open(gene_bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                genes.append({
                    'chr': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'gene': fields[3]
                })

    genes_df = pd.DataFrame(genes)
    print(f"✓ Loaded {len(genes_df)} genes")
    return genes_df

def create_gene_bed_from_gff3():
    """Create gene BED from GFF3 if BED doesn't exist"""
    gff3_file = "/home/chbope/extension/nWGS_manuscript_data/data/reference/gencode.v48.annotation.gff3"

    if not os.path.exists(gff3_file):
        print(f"⚠️  GFF3 file not found: {gff3_file}")
        return pd.DataFrame()

    print(f"Creating gene BED from GFF3: {gff3_file}")

    genes = []
    with gzip.open(gff3_file, 'rt') if gff3_file.endswith('.gz') else open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]

            # Extract gene name
            gene_name = None
            for attr in attributes.split(';'):
                if attr.startswith('gene_name='):
                    gene_name = attr.split('=')[1]
                    break

            if gene_name:
                genes.append({
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'gene': gene_name
                })

    genes_df = pd.DataFrame(genes)

    # Save as BED in svmetaunique folder
    bed_dir = "/home/chbope/extension/script/svmetaunique/external_datasets"
    os.makedirs(bed_dir, exist_ok=True)
    bed_file = f"{bed_dir}/refseq_genes_hg38.bed"

    with open(bed_file, 'w') as f:
        for _, row in genes_df.iterrows():
            f.write(f"{row['chr']}\t{row['start']}\t{row['end']}\t{row['gene']}\n")

    print(f"✓ Created gene BED: {bed_file}")
    print(f"✓ Extracted {len(genes_df)} genes")

    return genes_df

def overlap_sv_with_genes(sv, genes_df):
    """
    Find genes that overlap with an SV
    Returns list of gene names
    """
    sv_chr = sv['chr']
    sv_start = sv['pos']
    sv_end = sv['end']

    # Filter genes on same chromosome
    chr_genes = genes_df[genes_df['chr'] == sv_chr]

    # Find overlapping genes
    overlapping = chr_genes[
        (chr_genes['start'] <= sv_end) &
        (chr_genes['end'] >= sv_start)
    ]

    return overlapping['gene'].unique().tolist()

def parse_merged_vcf_and_annotate_genes(vcf_file, genes_df):
    """
    Parse merged VCF and create gene-level matrix with CONSERVATIVE counting
    Conservative: Max 1 SV per gene per sample (binary)
    """
    print(f"\nParsing merged VCF: {vcf_file}")

    # Data structure: {gene: {sample: has_sv (0 or 1)}}
    gene_sample_matrix = defaultdict(lambda: defaultdict(int))

    # Track SV counts per gene per sample (for statistics)
    gene_sample_sv_counts = defaultdict(lambda: defaultdict(int))

    # Track SV types per gene per sample
    gene_sample_sv_types = defaultdict(lambda: defaultdict(lambda: {'DEL': 0, 'DUP': 0, 'INV': 0, 'BND': 0}))

    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'

    with open_func(vcf_file, mode) as f:
        samples = []

        for line in f:
            if line.startswith('##'):
                continue

            if line.startswith('#CHROM'):
                # Header line - extract sample names
                fields = line.strip().split('\t')
                samples = fields[9:]
                print(f"✓ Found {len(samples)} samples in VCF")
                continue

            # Parse variant line
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            svid = fields[2]
            info = fields[7]

            # Extract SV info
            svtype = None
            end = pos + 1000  # Default

            for item in info.split(';'):
                if item.startswith('SVTYPE='):
                    svtype = item.split('=')[1]
                elif item.startswith('END='):
                    end = int(item.split('=')[1])

            # Create SV object
            sv = {
                'chr': chrom,
                'pos': pos,
                'end': end,
                'svtype': svtype
            }

            # Find overlapping genes
            overlapping_genes = overlap_sv_with_genes(sv, genes_df)

            if not overlapping_genes:
                continue  # Skip SVs that don't overlap genes

            # Parse genotypes
            genotypes = fields[9:]

            for sample_idx, gt in enumerate(genotypes):
                if sample_idx >= len(samples):
                    break

                sample_name = samples[sample_idx]

                # Check if this sample has the SV
                gt_field = gt.split(':')[0]

                if gt_field in ['0/0', './.', '0']:
                    continue  # Reference genotype

                # This sample has an SV - mark all affected genes
                for gene in overlapping_genes:
                    # Conservative counting: mark as 1 (present) regardless of how many SVs
                    gene_sample_matrix[gene][sample_name] = 1

                    # Track actual SV counts (for statistics only)
                    gene_sample_sv_counts[gene][sample_name] += 1

                    # Track SV types
                    if svtype in ['DEL', 'DUP', 'INV', 'BND']:
                        gene_sample_sv_types[gene][sample_name][svtype] += 1

    print(f"✓ Parsed {len(gene_sample_matrix)} genes with SVs")

    # Convert to DataFrame (binary matrix)
    matrix_data = []
    for gene in sorted(gene_sample_matrix.keys()):
        row = {'gene': gene}
        for sample in samples:
            row[sample] = gene_sample_matrix[gene].get(sample, 0)
        matrix_data.append(row)

    matrix_df = pd.DataFrame(matrix_data)
    matrix_df = matrix_df.set_index('gene')

    return matrix_df, samples, gene_sample_sv_counts, gene_sample_sv_types

def calculate_gene_frequencies_conservative(matrix_df, total_samples=200):
    """
    Calculate gene-level frequencies using CONSERVATIVE counting
    Each gene is 0 (no SV) or 1 (has SV) per sample
    """
    print("\nCalculating gene frequencies (conservative counting)...")

    gene_stats = []

    for gene in matrix_df.index:
        # Count how many samples have this gene affected
        n_samples_affected = (matrix_df.loc[gene] == 1).sum()

        # Conservative frequency: binary presence
        conservative_freq = n_samples_affected / total_samples

        gene_stats.append({
            'gene': gene,
            'n_samples_affected': n_samples_affected,
            'frequency': conservative_freq,
            'frequency_pct': f"{conservative_freq*100:.1f}%"
        })

    gene_freq_df = pd.DataFrame(gene_stats)
    gene_freq_df = gene_freq_df.sort_values('frequency', ascending=False)

    print(f"✓ Calculated frequencies for {len(gene_freq_df)} genes")

    return gene_freq_df

def add_sv_type_breakdown(gene_freq_df, gene_sample_sv_types, samples):
    """Add SV type breakdown to gene frequency table"""

    gene_type_data = []

    for gene in gene_freq_df['gene']:
        del_count = sum(gene_sample_sv_types[gene][s]['DEL'] for s in samples)
        dup_count = sum(gene_sample_sv_types[gene][s]['DUP'] for s in samples)
        inv_count = sum(gene_sample_sv_types[gene][s]['INV'] for s in samples)
        bnd_count = sum(gene_sample_sv_types[gene][s]['BND'] for s in samples)

        gene_type_data.append({
            'gene': gene,
            'n_DEL': del_count,
            'n_DUP': dup_count,
            'n_INV': inv_count,
            'n_BND': bnd_count
        })

    type_df = pd.DataFrame(gene_type_data)

    # Merge with frequency table
    result_df = gene_freq_df.merge(type_df, on='gene', how='left')

    return result_df

def load_driver_genes(driver_file):
    """Load list of known GBM driver genes"""
    if not os.path.exists(driver_file):
        print(f"⚠️  Driver gene file not found: {driver_file}")
        return []

    with open(driver_file, 'r') as f:
        drivers = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    print(f"✓ Loaded {len(drivers)} driver genes")
    return drivers

def main():
    """Main analysis pipeline"""

    print("="*80)
    print("STEP 03: BUILD GENE MATRIX WITH CONSERVATIVE COUNTING")
    print("="*80)
    print("\nApproach: Max 1 SV per gene per sample (binary: 0 or 1)")
    print("Following npae082.pdf methodology")
    print("="*80)
    print()

    # Create output directories
    os.makedirs(f"{OUTPUT_DIR}/matrices", exist_ok=True)
    os.makedirs(f"{OUTPUT_DIR}/genes", exist_ok=True)

    # Load gene annotations
    genes_df = load_gene_annotations(GENE_ANNOTATION_BED)

    if genes_df.empty:
        print("⚠️  No gene annotations loaded. Exiting.")
        return

    # Parse merged VCF and create gene matrix
    matrix_df, samples, gene_sv_counts, gene_sv_types = parse_merged_vcf_and_annotate_genes(
        MERGED_VCF,
        genes_df
    )

    # Save binary matrix
    matrix_file = f"{OUTPUT_DIR}/matrices/gene_sample_matrix_conservative.csv"
    matrix_df.to_csv(matrix_file)
    print(f"\n✓ Saved gene×sample matrix: {matrix_file}")
    print(f"  Dimensions: {matrix_df.shape[0]} genes × {matrix_df.shape[1]} samples")

    # Calculate gene frequencies
    gene_freq_df = calculate_gene_frequencies_conservative(matrix_df, len(samples))

    # Add SV type breakdown
    gene_freq_df = add_sv_type_breakdown(gene_freq_df, gene_sv_types, samples)

    # Load driver genes
    driver_genes = load_driver_genes(GBM_DRIVERS)
    gene_freq_df['is_driver'] = gene_freq_df['gene'].isin(driver_genes)

    # Save gene frequencies
    freq_file = f"{OUTPUT_DIR}/genes/gene_frequencies_conservative.csv"
    gene_freq_df.to_csv(freq_file, index=False)
    print(f"✓ Saved gene frequencies: {freq_file}")

    # Print top genes
    print("\n" + "="*80)
    print("TOP 20 GENES (CONSERVATIVE COUNTING)")
    print("="*80)
    print("\n{:<15} {:>12} {:>10} {:>8} {:>8} {:>8} {:>8}".format(
        "Gene", "Frequency", "N Samples", "DEL", "DUP", "INV", "BND"
    ))
    print("-"*80)

    for _, row in gene_freq_df.head(20).iterrows():
        driver_mark = "★" if row['is_driver'] else " "
        print("{:<15} {:>12} {:>10} {:>8} {:>8} {:>8} {:>8} {}".format(
            row['gene'],
            row['frequency_pct'],
            row['n_samples_affected'],
            row['n_DEL'],
            row['n_DUP'],
            row['n_INV'],
            row['n_BND'],
            driver_mark
        ))

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total genes with SVs: {len(gene_freq_df)}")
    print(f"Total samples: {len(samples)}")
    print(f"Driver genes in top 20: {gene_freq_df.head(20)['is_driver'].sum()}")
    print(f"\nResults saved to: {OUTPUT_DIR}/")
    print("="*80)
    print("\nNext step: Run 04_external_dataset_comparison_conservative.py")

if __name__ == "__main__":
    main()
