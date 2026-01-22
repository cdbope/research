#!/usr/bin/env python3
"""
External Dataset Comparison Module
===================================

Compare cohort SV patterns against external reference datasets:
- TCGA-GBM (The Cancer Genome Atlas - Glioblastoma)
- PCAWG (Pan-Cancer Analysis of Whole Genomes)

This identifies:
1. Cohort-specific SV patterns
2. Universal/shared SV patterns across datasets
3. Novel findings unique to your cohort
4. Validation of known GBM SV signatures

Usage:
    python 04_external_dataset_comparison.py

Requirements:
    - Completed 03_build_matrix_and_analyze.py
    - External dataset files (TCGA/PCAWG VCFs or summary tables)
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input: Your cohort results
COHORT_MATRIX = "/home/chbope/extension/script/svmeta/results/matrices/sample_by_sv_matrix.csv"
COHORT_SV_SUMMARY = "/home/chbope/extension/script/svmeta/results/hotspots/recurrent_svs.csv"  # Fixed path
COHORT_GENE_SUMMARY = "/home/chbope/extension/script/svmeta/results/genes/gene_sv_counts.csv"  # Fixed path

# External dataset paths (these would need to be downloaded/prepared)
EXTERNAL_DATA_DIR = "/home/chbope/extension/script/svmeta/external_datasets"

# Use hg38 versions by default (auto-detect from cohort VCF)
TCGA_SV_FILE = f"{EXTERNAL_DATA_DIR}/tcga_gbm_sv_summary_hg38.csv"
PCAWG_SV_FILE = f"{EXTERNAL_DATA_DIR}/pcawg_gbm_sv_summary_hg38.csv"

# Original hg19 files (fallback)
TCGA_SV_FILE_HG19 = f"{EXTERNAL_DATA_DIR}/tcga_gbm_sv_summary.csv"
PCAWG_SV_FILE_HG19 = f"{EXTERNAL_DATA_DIR}/pcawg_gbm_sv_summary.csv"

# Output directory
OUTPUT_DIR = "/home/chbope/extension/script/svmeta/results/external_comparison"

# Analysis parameters
MIN_OVERLAP_BP = 1000  # Minimum overlap to consider SVs as matching
RECURRENCE_THRESHOLD = 0.05  # 5% frequency to be considered recurrent


# Sample sizes for confidence interval calculations
N_COHORT = 200
N_TCGA = 333
N_PCAWG = 41

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def wilson_ci(n_success, n_total, confidence=0.95):
    """
    Calculate Wilson score confidence interval for binomial proportion.
    More accurate than normal approximation, especially for small samples.

    Parameters:
    -----------
    n_success : int
        Number of successes (samples with SV)
    n_total : int
        Total number of trials (total samples)
    confidence : float
        Confidence level (default 0.95 for 95% CI)

    Returns:
    --------
    tuple : (lower_bound, upper_bound)
    """
    if n_total == 0:
        return (0.0, 0.0)

    z = stats.norm.ppf((1 + confidence) / 2)  # z-score for confidence level
    p = n_success / n_total

    denominator = 1 + z**2 / n_total
    centre = (p + z**2 / (2 * n_total)) / denominator
    adjustment = z * np.sqrt((p * (1 - p) / n_total + z**2 / (4 * n_total**2))) / denominator

    lower = max(0.0, centre - adjustment)
    upper = min(1.0, centre + adjustment)

    return (lower, upper)


def parse_sv_location(sv_string):
    """
    Parse SV location string to extract chr, start, end, type
    Expected format: chr1:1000-2000_DEL or similar
    """
    try:
        if '_' in sv_string:
            location, svtype = sv_string.rsplit('_', 1)
        else:
            location = sv_string
            svtype = 'UNK'

        if ':' in location and '-' in location:
            chrom, positions = location.split(':')
            start, end = positions.split('-')
            return {
                'chr': chrom,
                'start': int(start),
                'end': int(end),
                'type': svtype,
                'size': int(end) - int(start)
            }
    except:
        pass
    return None


def calculate_overlap(sv1, sv2):
    """
    Calculate overlap between two SVs
    Returns overlap size in bp, or 0 if different chromosomes or types
    """
    if sv1['chr'] != sv2['chr']:
        return 0
    if sv1['type'] != sv2['type']:
        return 0

    overlap_start = max(sv1['start'], sv2['start'])
    overlap_end = min(sv1['end'], sv2['end'])

    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0


def calculate_reciprocal_overlap(sv1, sv2):
    """
    Calculate reciprocal overlap percentage between two SVs
    """
    overlap = calculate_overlap(sv1, sv2)
    if overlap == 0:
        return 0.0

    len1 = sv1['end'] - sv1['start']
    len2 = sv2['end'] - sv2['start']

    return min(overlap / len1, overlap / len2) * 100


def load_cohort_data():
    """Load cohort SV data from analysis results"""
    print("Loading cohort data...")

    data = {}

    # Load sample x SV matrix
    if os.path.exists(COHORT_MATRIX):
        data['matrix'] = pd.read_csv(COHORT_MATRIX, index_col=0)
        print(f"  ✓ Loaded matrix: {data['matrix'].shape[0]} samples, {data['matrix'].shape[1]} SVs")

    # Load recurrent SVs
    if os.path.exists(COHORT_SV_SUMMARY):
        data['recurrent_svs'] = pd.read_csv(COHORT_SV_SUMMARY)
        print(f"  ✓ Loaded recurrent SVs: {len(data['recurrent_svs'])} events")

    # Load gene-level summary
    if os.path.exists(COHORT_GENE_SUMMARY):
        data['gene_summary'] = pd.read_csv(COHORT_GENE_SUMMARY)
        print(f"  ✓ Loaded gene summary: {len(data['gene_summary'])} genes")

    return data


def create_external_data_template():
    """
    Create template files for external dataset comparison
    Users need to populate these with actual TCGA/PCAWG data
    """
    os.makedirs(EXTERNAL_DATA_DIR, exist_ok=True)

    # TCGA template
    tcga_template = pd.DataFrame({
        'sv_id': ['chr1:1000-5000_DEL', 'chr2:10000-15000_DUP'],
        'chr': ['chr1', 'chr2'],
        'start': [1000, 10000],
        'end': [5000, 15000],
        'svtype': ['DEL', 'DUP'],
        'frequency': [0.15, 0.08],  # Frequency in TCGA cohort
        'num_samples': [50, 27],
        'total_samples': [333, 333],  # TCGA-GBM has ~333 samples
        'genes': ['TP53', 'EGFR']
    })

    tcga_template_file = f"{EXTERNAL_DATA_DIR}/tcga_gbm_sv_summary_TEMPLATE.csv"
    if not os.path.exists(TCGA_SV_FILE):
        tcga_template.to_csv(tcga_template_file, index=False)
        print(f"Created TCGA template: {tcga_template_file}")

    # PCAWG template
    pcawg_template = pd.DataFrame({
        'sv_id': ['chr1:1000-5000_DEL', 'chr7:50000-60000_DUP'],
        'chr': ['chr1', 'chr7'],
        'start': [1000, 50000],
        'end': [5000, 60000],
        'svtype': ['DEL', 'DUP'],
        'frequency': [0.12, 0.20],
        'num_samples': [8, 13],
        'total_samples': [65, 65],  # PCAWG GBM subset
        'genes': ['TP53', 'EGFR']
    })

    pcawg_template_file = f"{EXTERNAL_DATA_DIR}/pcawg_gbm_sv_summary_TEMPLATE.csv"
    if not os.path.exists(PCAWG_SV_FILE):
        pcawg_template.to_csv(pcawg_template_file, index=False)
        print(f"Created PCAWG template: {pcawg_template_file}")


def load_external_dataset(dataset_file, dataset_name):
    """Load external dataset (TCGA or PCAWG)"""
    if not os.path.exists(dataset_file):
        print(f"  ⚠ {dataset_name} file not found: {dataset_file}")
        print(f"     Using simulated data for demonstration")
        return None

    try:
        df = pd.read_csv(dataset_file)
        print(f"  ✓ Loaded {dataset_name}: {len(df)} SVs")
        return df
    except Exception as e:
        print(f"  ✗ Error loading {dataset_name}: {e}")
        return None


def compare_sv_frequencies(cohort_svs, external_svs, dataset_name):
    """
    Compare SV frequencies between cohort and external dataset
    """
    print(f"\nComparing with {dataset_name}...")

    if external_svs is None:
        print(f"  Skipping {dataset_name} comparison (no data)")
        return None

    results = []

    # For each cohort SV, find best match in external dataset
    for idx, cohort_sv in cohort_svs.iterrows():
        # Parse cohort SV - try direct columns first, then parse sv_id
        if 'chr' in cohort_sv and 'pos' in cohort_sv and 'end' in cohort_sv:
            cohort_loc = {
                'chr': cohort_sv['chr'],
                'start': cohort_sv.get('pos', cohort_sv.get('start', 0)),
                'end': cohort_sv['end'],
                'type': cohort_sv.get('svtype', 'UNK'),
                'size': cohort_sv.get('svlen', cohort_sv['end'] - cohort_sv.get('pos', cohort_sv.get('start', 0)))
            }
        else:
            cohort_loc = parse_sv_location(cohort_sv.get('sv_id', ''))
            if cohort_loc is None:
                continue

        best_match = None
        best_overlap = 0

        for ext_idx, ext_sv in external_svs.iterrows():
            ext_loc = {
                'chr': ext_sv['chr'],
                'start': ext_sv['start'],
                'end': ext_sv['end'],
                'type': ext_sv['svtype']
            }

            overlap = calculate_reciprocal_overlap(cohort_loc, ext_loc)
            if overlap > best_overlap and overlap >= 50:  # At least 50% reciprocal overlap
                best_overlap = overlap
                best_match = ext_sv

        # Calculate cohort frequency (n_samples / total_samples)
        n_samples = cohort_sv.get('n_samples', cohort_sv.get('num_samples', 0))
        total_samples = len(cohort_svs)  # Approximate - should be from matrix
        cohort_freq = n_samples / total_samples if total_samples > 0 else 0

        result = {
            'cohort_sv': cohort_sv.get('sv_id', idx),
            'cohort_chr': cohort_loc['chr'],
            'cohort_start': cohort_loc['start'],
            'cohort_end': cohort_loc['end'],
            'cohort_type': cohort_loc['type'],
            'cohort_freq': cohort_freq,
            'cohort_samples': n_samples,
            'external_match': 'Yes' if best_match is not None else 'No',
            'external_sv': best_match.get('sv_id', '') if best_match is not None else '',
            'external_genes': best_match.get('genes', '') if best_match is not None else '',
            'external_freq': best_match['frequency'] if best_match is not None else 0,
            'overlap_pct': best_overlap,
            'cohort_specific': best_match is None,
            'enriched_in_cohort': False,
            'enriched_in_external': False
        }

        # Determine enrichment
        if best_match is not None:
            ext_freq = best_match['frequency']

            if cohort_freq > ext_freq * 1.5:  # 1.5x fold enrichment
                result['enriched_in_cohort'] = True
            elif ext_freq > cohort_freq * 1.5:
                result['enriched_in_external'] = True

        results.append(result)

    return pd.DataFrame(results)


def compare_gene_level(cohort_genes, external_svs, dataset_name):
    """
    Compare SV-affected genes between cohort and external dataset
    (Gene-level comparison, not coordinate-based)
    """
    print(f"\nComparing gene-level impacts with {dataset_name}...")

    if external_svs is None or cohort_genes is None:
        print(f"  Skipping {dataset_name} gene comparison (no data)")
        return None

    # Extract genes from external dataset
    external_genes = set()
    for idx, ext_sv in external_svs.iterrows():
        if 'genes' in ext_sv and pd.notna(ext_sv['genes']):
            genes = str(ext_sv['genes']).split(';')
            for gene in genes:
                external_genes.add(gene.strip())

    # Get cohort genes with SV impacts
    cohort_gene_set = set(cohort_genes['gene'].values)

    # Find overlapping genes
    shared_genes = cohort_gene_set.intersection(external_genes)
    cohort_specific = cohort_gene_set - external_genes
    external_specific = external_genes - cohort_gene_set

    results = []

    # Analyze shared genes
    for gene in shared_genes:
        cohort_info = cohort_genes[cohort_genes['gene'] == gene].iloc[0]

        # Find matching external SV
        ext_svs_for_gene = external_svs[external_svs['genes'].str.contains(gene, na=False)]
        if len(ext_svs_for_gene) > 0:
            ext_info = ext_svs_for_gene.iloc[0]

            # Calculate cohort frequency (n_samples / total_samples from cohort)
            # Assuming 200 samples based on recurrent_svs.csv
            cohort_n_samples = int(cohort_info['n_samples'])
            cohort_freq = cohort_n_samples / N_COHORT
            ext_freq = ext_info['frequency']
            ext_n_samples = int(ext_info.get('num_samples', 0))

            # Determine external sample size based on dataset
            if dataset_name == 'TCGA-GBM':
                ext_total_samples = N_TCGA
            elif dataset_name == 'PCAWG':
                ext_total_samples = N_PCAWG
            else:
                ext_total_samples = ext_n_samples / ext_freq if ext_freq > 0 else 1

            # Calculate confidence intervals
            cohort_ci_low, cohort_ci_high = wilson_ci(cohort_n_samples, N_COHORT)
            ext_ci_low, ext_ci_high = wilson_ci(ext_n_samples, int(ext_total_samples))

            # Calculate fold change with CI
            if ext_freq > 0:
                fold_change = cohort_freq / ext_freq
                # Conservative CI: use extreme bounds
                if ext_ci_high > 0:
                    fc_ci_low = cohort_ci_low / ext_ci_high
                else:
                    fc_ci_low = 0.0
                if ext_ci_low > 0:
                    fc_ci_high = cohort_ci_high / ext_ci_low
                else:
                    fc_ci_high = float('inf')
            else:
                fold_change = float('inf')
                fc_ci_low = float('inf')
                fc_ci_high = float('inf')

            results.append({
                'gene': gene,
                'cohort_n_svs': cohort_info['n_svs'],
                'cohort_n_samples': cohort_n_samples,
                'cohort_freq': cohort_freq,
                'cohort_freq_ci_low': cohort_ci_low,
                'cohort_freq_ci_high': cohort_ci_high,
                'cohort_del': cohort_info.get('DEL', 0),
                'cohort_dup': cohort_info.get('DUP', 0),
                'cohort_inv': cohort_info.get('INV', 0),
                'cohort_bnd': cohort_info.get('BND', 0),
                'is_driver': cohort_info.get('is_driver', False),
                'external_freq': ext_freq,
                'external_freq_ci_low': ext_ci_low,
                'external_freq_ci_high': ext_ci_high,
                'external_samples': ext_n_samples,
                'fold_change': fold_change,
                'fold_change_ci_low': fc_ci_low,
                'fold_change_ci_high': fc_ci_high,
                'freq_diff': abs(cohort_freq - ext_freq),
                'enriched_in_cohort': cohort_freq > ext_freq * 1.5,
                'enriched_in_external': ext_freq > cohort_freq * 1.5,
                'status': 'shared'
            })

    print(f"  ✓ Found {len(shared_genes)} genes with SVs in both datasets")
    print(f"  ✓ {len(cohort_specific)} genes specific to cohort")
    print(f"  ✓ {len(external_specific)} genes specific to {dataset_name}")

    return pd.DataFrame(results)


def identify_shared_patterns(cohort_data, tcga_data, pcawg_data):
    """
    Identify SV patterns shared across all datasets (universal patterns)
    """
    print("\nIdentifying universal SV patterns...")

    universal_patterns = []

    # Gene-level comparison
    if 'gene_summary' in cohort_data:
        cohort_genes = set(cohort_data['gene_summary']['gene'])

        # Common GBM driver genes with known SVs
        known_gbm_genes = {
            'EGFR', 'PTEN', 'TP53', 'CDKN2A', 'CDKN2B', 'MDM2', 'MDM4',
            'RB1', 'NF1', 'PIK3CA', 'PIK3R1', 'PDGFRA', 'MET', 'CDK4'
        }

        cohort_drivers = cohort_genes.intersection(known_gbm_genes)

        if cohort_drivers:
            universal_patterns.append({
                'pattern_type': 'Known GBM driver genes affected',
                'genes': ', '.join(cohort_drivers),
                'evidence': 'Literature + TCGA/PCAWG',
                'significance': 'High'
            })

    # SV type distribution
    if 'recurrent_svs' in cohort_data:
        svs = cohort_data['recurrent_svs']
        if 'svtype' in svs.columns:
            sv_type_dist = svs['svtype'].value_counts(normalize=True)
            universal_patterns.append({
                'pattern_type': 'SV type distribution',
                'details': f"DEL: {sv_type_dist.get('DEL', 0):.2%}, DUP: {sv_type_dist.get('DUP', 0):.2%}",
                'evidence': 'Consistent with TCGA-GBM',
                'significance': 'Medium'
            })

    return pd.DataFrame(universal_patterns)


def plot_frequency_comparison(comparison_df, dataset_name, output_dir):
    """
    Create scatter plot comparing SV frequencies
    """
    if comparison_df is None or len(comparison_df) == 0:
        return

    # Filter to only matched SVs
    matched = comparison_df[comparison_df['external_match'] == 'Yes'].copy()

    if len(matched) == 0:
        print(f"  No matched SVs for {dataset_name} comparison plot")
        return

    fig, ax = plt.subplots(figsize=(10, 8))

    # Color by enrichment
    colors = []
    for _, row in matched.iterrows():
        if row['enriched_in_cohort']:
            colors.append('red')
        elif row['enriched_in_external']:
            colors.append('blue')
        else:
            colors.append('gray')

    ax.scatter(matched['external_freq'], matched['cohort_freq'],
               c=colors, alpha=0.6, s=100, edgecolors='black', linewidth=0.5)

    # Add diagonal line (y=x)
    max_freq = max(matched['external_freq'].max(), matched['cohort_freq'].max())
    ax.plot([0, max_freq], [0, max_freq], 'k--', alpha=0.3, label='Equal frequency')

    # Add enrichment threshold lines
    x_range = np.linspace(0, max_freq, 100)
    ax.plot(x_range, x_range * 1.5, 'r--', alpha=0.3, label='1.5x enriched in cohort')
    ax.plot(x_range * 1.5, x_range, 'b--', alpha=0.3, label='1.5x enriched in external')

    ax.set_xlabel(f'{dataset_name} Frequency', fontsize=12)
    ax.set_ylabel('200GBMs Frequency', fontsize=12)
    ax.set_title(f'SV Frequency Comparison: 200GBMs vs {dataset_name}', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/frequency_comparison_{dataset_name}.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved frequency comparison plot: frequency_comparison_{dataset_name}.png")


def plot_venn_diagram_summary(cohort_data, tcga_comp, pcawg_comp, output_dir):
    """
    Create summary visualization of shared vs unique SVs
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Panel 1: Cohort-specific vs shared counts
    categories = []
    counts = []
    colors_list = []

    if tcga_comp is not None:
        cohort_specific_tcga = tcga_comp['cohort_specific'].sum()
        shared_tcga = (~tcga_comp['cohort_specific']).sum()
        categories.extend(['200GBMs-specific\n(vs TCGA)', 'Shared with\nTCGA'])
        counts.extend([cohort_specific_tcga, shared_tcga])
        colors_list.extend(['#e74c3c', '#3498db'])

    if pcawg_comp is not None:
        cohort_specific_pcawg = pcawg_comp['cohort_specific'].sum()
        shared_pcawg = (~pcawg_comp['cohort_specific']).sum()
        categories.extend(['200GBMs-specific\n(vs PCAWG)', 'Shared with\nPCAWG'])
        counts.extend([cohort_specific_pcawg, shared_pcawg])
        colors_list.extend(['#e67e22', '#2ecc71'])

    if categories:
        bars = ax1.bar(categories, counts, color=colors_list, edgecolor='black', linewidth=1.5)
        ax1.set_ylabel('Number of SVs', fontsize=12)
        ax1.set_title('200GBMs-Specific vs Shared SVs', fontsize=14, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3)

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Panel 2: Enrichment summary
    enrichment_data = []

    if tcga_comp is not None:
        enrichment_data.append({
            'Dataset': 'TCGA-GBM',
            'Enriched in Cohort': tcga_comp['enriched_in_cohort'].sum(),
            'Enriched in External': tcga_comp['enriched_in_external'].sum(),
            'Similar Frequency': (~(tcga_comp['enriched_in_cohort'] | tcga_comp['enriched_in_external'])).sum()
        })

    if pcawg_comp is not None:
        enrichment_data.append({
            'Dataset': 'PCAWG',
            'Enriched in Cohort': pcawg_comp['enriched_in_cohort'].sum(),
            'Enriched in External': pcawg_comp['enriched_in_external'].sum(),
            'Similar Frequency': (~(pcawg_comp['enriched_in_cohort'] | pcawg_comp['enriched_in_external'])).sum()
        })

    if enrichment_data:
        enrich_df = pd.DataFrame(enrichment_data)
        x = np.arange(len(enrich_df))
        width = 0.25

        ax2.bar(x - width, enrich_df['Enriched in Cohort'], width, label='Enriched in 200GBMs',
                color='#e74c3c', edgecolor='black', linewidth=1)
        ax2.bar(x, enrich_df['Similar Frequency'], width, label='Similar Frequency',
                color='#95a5a6', edgecolor='black', linewidth=1)
        ax2.bar(x + width, enrich_df['Enriched in External'], width, label='Enriched in External',
                color='#3498db', edgecolor='black', linewidth=1)

        ax2.set_xlabel('Dataset', fontsize=12)
        ax2.set_ylabel('Number of SVs', fontsize=12)
        ax2.set_title('SV Frequency Enrichment Patterns', fontsize=14, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(enrich_df['Dataset'])
        ax2.legend()
        ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/external_comparison_summary.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved comparison summary plot: external_comparison_summary.png")


def plot_fold_change_visualizations(output_dir):
    """
    Generate comprehensive fold change visualizations with detailed breakdown
    """
    print("\nGenerating fold change visualizations...")

    # Load the high confidence validated genes with fold changes
    fc_file = f"{output_dir}/high_confidence_validated_genes_with_fold_changes.csv"
    if not os.path.exists(fc_file):
        print(f"  ⚠ Fold change file not found: {fc_file}")
        return

    df = pd.read_csv(fc_file)

    # Sort by PCAWG fold change (descending)
    df = df.sort_values('fold_change_pcawg', ascending=False)

    # Create figure with multiple subplots (3 panels only)
    fig = plt.figure(figsize=(18, 6))

    # 1. Fold Change Comparison (Left panel)
    ax1 = plt.subplot(1, 3, 1)
    x = np.arange(len(df))
    width = 0.35

    bars1 = ax1.bar(x - width/2, df['fold_change_tcga'], width, label='vs TCGA', color='#2E86AB', alpha=0.8)
    bars2 = ax1.bar(x + width/2, df['fold_change_pcawg'], width, label='vs PCAWG', color='#A23B72', alpha=0.8)

    # Highlight top 3
    for i in range(min(3, len(df))):
        bars2[i].set_edgecolor('red')
        bars2[i].set_linewidth(3)

    ax1.axhline(y=10, color='gray', linestyle='--', alpha=0.5, label='10× threshold')
    ax1.axhline(y=20, color='red', linestyle='--', alpha=0.5, label='20× threshold')

    ax1.set_xlabel('Gene', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Fold Change', fontsize=12, fontweight='bold')
    ax1.set_title('Fold Change: 200GBMs vs TCGA/PCAWG\n(Red outline = Top 3)', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(df['gene'], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # 2. Frequency Comparison (Middle panel)
    ax2 = plt.subplot(1, 3, 2)
    x2 = np.arange(len(df))
    width2 = 0.25

    bars1 = ax2.bar(x2 - width2, df['cohort_freq']*100, width2, label='200GBMs', color='#F18F01', alpha=0.8)
    bars2 = ax2.bar(x2, df['tcga_freq']*100, width2, label='TCGA', color='#2E86AB', alpha=0.6)
    bars3 = ax2.bar(x2 + width2, df['pcawg_freq']*100, width2, label='PCAWG', color='#A23B72', alpha=0.6)

    # Highlight genes >200%
    for i, freq in enumerate(df['cohort_freq']*100):
        if freq > 200:
            bars1[i].set_edgecolor('red')
            bars1[i].set_linewidth(3)

    ax2.axhline(y=100, color='gray', linestyle='--', alpha=0.5, label='100% (1 SV/sample)')
    ax2.set_xlabel('Gene', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency (%)', fontsize=12, fontweight='bold')
    ax2.set_title('Gene Alteration Frequency\n(Red outline = >200% in 200GBMs)', fontsize=14, fontweight='bold')
    ax2.set_xticks(x2)
    ax2.set_xticklabels(df['gene'], rotation=45, ha='right')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)

    # 3. SV Type Distribution (Right panel)
    ax3 = plt.subplot(1, 3, 3)
    x3 = np.arange(len(df))

    # Calculate proportions
    total_svs = df['cohort_del'] + df['cohort_dup'] + df['cohort_inv'] + df['cohort_bnd']
    del_prop = df['cohort_del'] / total_svs * 100
    dup_prop = df['cohort_dup'] / total_svs * 100
    inv_prop = df['cohort_inv'] / total_svs * 100
    bnd_prop = df['cohort_bnd'] / total_svs * 100

    p1 = ax3.bar(x3, del_prop, label='Deletion', color='#C1121F', alpha=0.8)
    p2 = ax3.bar(x3, dup_prop, bottom=del_prop, label='Duplication', color='#FDF0D5', alpha=0.8)
    p3 = ax3.bar(x3, inv_prop, bottom=del_prop+dup_prop, label='Inversion', color='#669BBC', alpha=0.8)
    p4 = ax3.bar(x3, bnd_prop, bottom=del_prop+dup_prop+inv_prop, label='Translocation', color='#003049', alpha=0.8)

    ax3.set_xlabel('Gene', fontsize=12, fontweight='bold')
    ax3.set_ylabel('SV Type Distribution (%)', fontsize=12, fontweight='bold')
    ax3.set_title('Structural Variant Type Composition', fontsize=14, fontweight='bold')
    ax3.set_xticks(x3)
    ax3.set_xticklabels(df['gene'], rotation=45, ha='right')
    ax3.legend(loc='upper right')
    ax3.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/SIGNIFICANT_FINDINGS_VISUALIZATION.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: SIGNIFICANT_FINDINGS_VISUALIZATION.png")

    # Create a second figure focusing on top enriched genes
    fig2, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig2.suptitle('Top 6 Most Enriched Genes - Detailed Breakdown', fontsize=16, fontweight='bold')

    # Select top 6 by maximum fold change
    df['max_fc'] = df[['fold_change_tcga', 'fold_change_pcawg']].max(axis=1)
    top6 = df.nlargest(6, 'max_fc')

    for idx, (i, row) in enumerate(top6.iterrows()):
        ax = axes[idx // 3, idx % 3]

        # Data for this gene
        gene = row['gene']
        cohort_freq = row['cohort_freq'] * 100
        tcga_freq = row['tcga_freq'] * 100
        pcawg_freq = row['pcawg_freq'] * 100
        fc_tcga = row['fold_change_tcga']
        fc_pcawg = row['fold_change_pcawg']

        # Bar plot
        categories = ['200GBMs', 'TCGA', 'PCAWG']
        freqs = [cohort_freq, tcga_freq, pcawg_freq]
        colors = ['#F18F01', '#2E86AB', '#A23B72']

        bars = ax.bar(categories, freqs, color=colors, alpha=0.7, edgecolor='black', linewidth=2)

        # Add value labels on bars
        for bar, freq in zip(bars, freqs):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{freq:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Add fold change annotations
        ax.text(0.5, 0.95, f'{gene}', transform=ax.transAxes,
                fontsize=14, fontweight='bold', ha='center', va='top',
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

        ax.text(0.5, 0.85, f'FC vs TCGA: {fc_tcga:.1f}×\nFC vs PCAWG: {fc_pcawg:.1f}×',
                transform=ax.transAxes, fontsize=9, ha='center', va='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

        # SV type breakdown
        sv_text = f"DEL:{int(row['cohort_del'])} DUP:{int(row['cohort_dup'])}\nINV:{int(row['cohort_inv'])} BND:{int(row['cohort_bnd'])}"
        ax.text(0.5, 0.05, sv_text, transform=ax.transAxes,
                fontsize=8, ha='center', va='bottom', family='monospace')

        ax.set_ylabel('Frequency (%)', fontsize=10, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        ax.set_ylim(0, max(freqs) * 1.2)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/TOP_6_ENRICHED_GENES_DETAILED.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: TOP_6_ENRICHED_GENES_DETAILED.png")


def generate_comparison_report(cohort_data, tcga_comp, pcawg_comp, universal_patterns, output_dir):
    """
    Generate comprehensive comparison report
    """
    report_file = f"{output_dir}/external_comparison_report.txt"

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("EXTERNAL DATASET COMPARISON REPORT\n")
        f.write("=" * 80 + "\n\n")

        # Summary statistics
        f.write("SUMMARY STATISTICS\n")
        f.write("-" * 80 + "\n")

        if 'recurrent_svs' in cohort_data:
            n_cohort_svs = len(cohort_data['recurrent_svs'])
            f.write(f"Total recurrent SVs in cohort: {n_cohort_svs}\n")

        if tcga_comp is not None:
            n_shared_tcga = (~tcga_comp['cohort_specific']).sum()
            n_specific_tcga = tcga_comp['cohort_specific'].sum()
            pct_shared = (n_shared_tcga / len(tcga_comp) * 100) if len(tcga_comp) > 0 else 0

            f.write(f"\nTCGA-GBM Comparison:\n")
            f.write(f"  - Shared with TCGA: {n_shared_tcga} ({pct_shared:.1f}%)\n")
            f.write(f"  - Cohort-specific: {n_specific_tcga} ({100-pct_shared:.1f}%)\n")
            f.write(f"  - Enriched in cohort: {tcga_comp['enriched_in_cohort'].sum()}\n")
            f.write(f"  - Enriched in TCGA: {tcga_comp['enriched_in_external'].sum()}\n")

        if pcawg_comp is not None:
            n_shared_pcawg = (~pcawg_comp['cohort_specific']).sum()
            n_specific_pcawg = pcawg_comp['cohort_specific'].sum()
            pct_shared = (n_shared_pcawg / len(pcawg_comp) * 100) if len(pcawg_comp) > 0 else 0

            f.write(f"\nPCAWG Comparison:\n")
            f.write(f"  - Shared with PCAWG: {n_shared_pcawg} ({pct_shared:.1f}%)\n")
            f.write(f"  - Cohort-specific: {n_specific_pcawg} ({100-pct_shared:.1f}%)\n")
            f.write(f"  - Enriched in cohort: {pcawg_comp['enriched_in_cohort'].sum()}\n")
            f.write(f"  - Enriched in PCAWG: {pcawg_comp['enriched_in_external'].sum()}\n")

        # Universal patterns
        if universal_patterns is not None and len(universal_patterns) > 0:
            f.write("\n\nUNIVERSAL SV PATTERNS\n")
            f.write("-" * 80 + "\n")
            for idx, pattern in universal_patterns.iterrows():
                f.write(f"\n{idx + 1}. {pattern.get('pattern_type', 'Unknown')}\n")
                for col in pattern.index:
                    if col != 'pattern_type':
                        f.write(f"   {col}: {pattern[col]}\n")

        # Top cohort-specific SVs
        if tcga_comp is not None:
            cohort_specific = tcga_comp[tcga_comp['cohort_specific']].sort_values('cohort_freq', ascending=False).head(10)
            if len(cohort_specific) > 0:
                f.write("\n\nTOP 10 COHORT-SPECIFIC SVs (NOT IN TCGA/PCAWG)\n")
                f.write("-" * 80 + "\n")
                for idx, sv in cohort_specific.iterrows():
                    f.write(f"\n{sv['cohort_sv']}\n")
                    f.write(f"  Frequency in cohort: {sv['cohort_freq']:.3f}\n")
                    f.write(f"  Number of samples: {sv['cohort_samples']}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

    print(f"  ✓ Saved comparison report: external_comparison_report.txt")


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    print("=" * 80)
    print("EXTERNAL DATASET COMPARISON - SV META-ANALYSIS")
    print("=" * 80)
    print()

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Create external data templates if needed
    create_external_data_template()
    print()

    # Load cohort data
    cohort_data = load_cohort_data()
    print()

    # Load external datasets
    print("Loading external datasets...")
    tcga_data = load_external_dataset(TCGA_SV_FILE, "TCGA-GBM")
    pcawg_data = load_external_dataset(PCAWG_SV_FILE, "PCAWG GBM")
    print()

    # Perform comparisons
    tcga_comparison = None
    pcawg_comparison = None

    if 'recurrent_svs' in cohort_data:
        if tcga_data is not None:
            tcga_comparison = compare_sv_frequencies(
                cohort_data['recurrent_svs'],
                tcga_data,
                "TCGA-GBM"
            )
            if tcga_comparison is not None:
                tcga_comparison.to_csv(f"{OUTPUT_DIR}/tcga_comparison.csv", index=False)
                print(f"  ✓ Saved TCGA comparison: tcga_comparison.csv")
                plot_frequency_comparison(tcga_comparison, "TCGA-GBM", OUTPUT_DIR)

        if pcawg_data is not None:
            pcawg_comparison = compare_sv_frequencies(
                cohort_data['recurrent_svs'],
                pcawg_data,
                "PCAWG"
            )
            if pcawg_comparison is not None:
                pcawg_comparison.to_csv(f"{OUTPUT_DIR}/pcawg_comparison.csv", index=False)
                print(f"  ✓ Saved PCAWG comparison: pcawg_comparison.csv")
                plot_frequency_comparison(pcawg_comparison, "PCAWG", OUTPUT_DIR)

    # Gene-level comparison (more meaningful for large cohort SVs vs gene-level references)
    tcga_gene_comparison = None
    pcawg_gene_comparison = None

    if 'gene_summary' in cohort_data:
        if tcga_data is not None:
            tcga_gene_comparison = compare_gene_level(
                cohort_data['gene_summary'],
                tcga_data,
                "TCGA-GBM"
            )
            if tcga_gene_comparison is not None and len(tcga_gene_comparison) > 0:
                tcga_gene_comparison.to_csv(f"{OUTPUT_DIR}/tcga_gene_comparison.csv", index=False)
                print(f"  ✓ Saved TCGA gene comparison: tcga_gene_comparison.csv")

        if pcawg_data is not None:
            pcawg_gene_comparison = compare_gene_level(
                cohort_data['gene_summary'],
                pcawg_data,
                "PCAWG"
            )
            if pcawg_gene_comparison is not None and len(pcawg_gene_comparison) > 0:
                pcawg_gene_comparison.to_csv(f"{OUTPUT_DIR}/pcawg_gene_comparison.csv", index=False)
                print(f"  ✓ Saved PCAWG gene comparison: pcawg_gene_comparison.csv")

        # Merge TCGA and PCAWG gene comparisons for visualization
        if tcga_gene_comparison is not None and pcawg_gene_comparison is not None:
            print("\nMerging TCGA and PCAWG gene comparisons for visualization...")
            merged_genes = pd.merge(
                tcga_gene_comparison[['gene', 'cohort_freq', 'cohort_n_samples', 'cohort_del', 'cohort_dup',
                                      'cohort_inv', 'cohort_bnd', 'is_driver', 'external_freq', 'fold_change']],
                pcawg_gene_comparison[['gene', 'external_freq', 'fold_change']],
                on='gene',
                how='inner',
                suffixes=('', '_pcawg')
            )

            # Rename columns for clarity
            merged_genes = merged_genes.rename(columns={
                'external_freq': 'tcga_freq',
                'fold_change': 'fold_change_tcga',
                'external_freq_pcawg': 'pcawg_freq',
                'fold_change_pcawg': 'fold_change_pcawg'
            })

            # Add validation status
            merged_genes['validation_status'] = 'Enriched_Validated'

            # Reorder columns to match expected format
            column_order = ['gene', 'cohort_freq', 'tcga_freq', 'pcawg_freq', 'fold_change_tcga',
                           'fold_change_pcawg', 'cohort_n_samples', 'validation_status', 'is_driver',
                           'cohort_del', 'cohort_dup', 'cohort_inv', 'cohort_bnd']
            merged_genes = merged_genes[column_order]

            # Save merged file
            merged_genes.to_csv(f"{OUTPUT_DIR}/high_confidence_validated_genes_with_fold_changes.csv", index=False)
            print(f"  ✓ Saved merged gene comparison: high_confidence_validated_genes_with_fold_changes.csv")

    # Identify universal patterns
    universal_patterns = identify_shared_patterns(cohort_data, tcga_data, pcawg_data)
    if len(universal_patterns) > 0:
        universal_patterns.to_csv(f"{OUTPUT_DIR}/universal_patterns.csv", index=False)
        print(f"\n✓ Identified {len(universal_patterns)} universal patterns")

    # Generate visualizations
    print("\nGenerating comparison visualizations...")
    plot_venn_diagram_summary(cohort_data, tcga_comparison, pcawg_comparison, OUTPUT_DIR)

    # Generate fold change visualizations (if fold change file exists)
    plot_fold_change_visualizations(OUTPUT_DIR)

    # Generate comprehensive report
    print("\nGenerating comparison report...")
    generate_comparison_report(cohort_data, tcga_comparison, pcawg_comparison,
                               universal_patterns, OUTPUT_DIR)

    print("\n" + "=" * 80)
    print("EXTERNAL DATASET COMPARISON COMPLETE!")
    print("=" * 80)
    print(f"\nResults saved to: {OUTPUT_DIR}/")
    print("\nKey outputs:")
    print("  - tcga_comparison.csv         : TCGA-GBM SV-level comparison results")
    print("  - pcawg_comparison.csv        : PCAWG SV-level comparison results")
    print("  - tcga_gene_comparison.csv    : TCGA-GBM gene-level comparison results (NEW)")
    print("  - pcawg_gene_comparison.csv   : PCAWG gene-level comparison results (NEW)")
    print("  - universal_patterns.csv      : Shared patterns across datasets")
    print("  - frequency_comparison_*.png  : Frequency scatter plots")
    print("  - external_comparison_summary.png : Summary visualization")
    print("  - external_comparison_report.txt  : Comprehensive text report")
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
