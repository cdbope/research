#!/usr/bin/env python3
"""
Pathway Enrichment Analysis for SV-Affected Genes
==================================================

This script performs pathway enrichment analysis on genes affected by
structural variants, using multiple enrichment methods:

1. Gene Ontology (GO) enrichment
2. KEGG pathway enrichment
3. Reactome pathway enrichment
4. Disease association enrichment
5. Transcription factor target enrichment

Methods:
- gprofiler-official (g:Profiler API) - Primary method
- Local hypergeometric test (as backup)
- Export for external tools (DAVID, Enrichr)

Usage:
    python 05_pathway_enrichment.py

Requirements:
    - pip install gprofiler-official
    - Results from 03_build_matrix_and_analyze.py

Output:
    - Pathway enrichment tables
    - Visualization plots
    - Gene lists for external tools
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from scipy import stats
import requests
import json
import warnings
warnings.filterwarnings('ignore')

# Try to import gprofiler
try:
    from gprofiler import GProfiler
    GPROFILER_AVAILABLE = True
except ImportError:
    print("Warning: gprofiler-official not installed. Using local analysis only.")
    print("Install with: pip install gprofiler-official")
    GPROFILER_AVAILABLE = False

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
GENE_RECURRENCE_FILE = "/home/chbope/extension/script/svmeta/results/genes/gene_sv_counts.csv"
SAMPLE_SUMMARY_FILE = "/home/chbope/extension/script/svmeta/results/matrices/sample_summary.csv"

# Output directory
OUTPUT_DIR = "/home/chbope/extension/script/svmeta/results/pathway_enrichment"

# Analysis parameters
TOP_N_GENES = 100  # Top N most affected genes for enrichment
MIN_GENE_RECURRENCE = 3  # Minimum samples with SV for gene to be included
PVALUE_THRESHOLD = 0.05  # Significance threshold
FDR_THRESHOLD = 0.1  # FDR threshold for multiple testing
USE_ENRICHR = True  # Use Enrichr API for additional enrichment

# GBM-specific pathways (for validation)
EXPECTED_GBM_PATHWAYS = {
    'RTK_RAS_PI3K': ['EGFR', 'PDGFRA', 'MET', 'PTEN', 'PIK3CA', 'PIK3R1', 'NF1', 'KRAS', 'NRAS'],
    'TP53_pathway': ['TP53', 'MDM2', 'MDM4', 'CDKN2A', 'CDKN2B'],
    'Cell_cycle': ['CDKN2A', 'CDKN2B', 'CDK4', 'CCND1', 'CCND2', 'CCND3', 'RB1'],
    'DNA_repair': ['ATM', 'ATR', 'BRCA1', 'BRCA2', 'PARP1', 'MLH1', 'MSH2'],
    'Chromatin': ['ARID1A', 'ARID1B', 'SMARCA4', 'SETD2', 'KDM6A', 'CREBBP']
}

# Core GBM signaling pathways (for detailed analysis)
CORE_GBM_PATHWAYS = {
    'RTK/EGFR Signaling': ['EGFR', 'PDGFRA', 'PDGFRB', 'MET', 'FGFR1', 'FGFR2', 'FGFR3', 'ERBB2', 'ERBB3'],
    'PI3K-AKT-mTOR': ['PIK3CA', 'PIK3R1', 'PIK3R2', 'AKT1', 'AKT2', 'AKT3', 'MTOR', 'PTEN', 'TSC1', 'TSC2'],
    'RAS-MAPK': ['KRAS', 'NRAS', 'HRAS', 'BRAF', 'RAF1', 'MAP2K1', 'MAP2K2', 'MAPK1', 'MAPK3', 'NF1'],
    'TP53 Pathway': ['TP53', 'MDM2', 'MDM4', 'CDKN2A', 'ATM', 'ATR', 'CHEK1', 'CHEK2'],
    'RB/Cell Cycle': ['RB1', 'CDKN2A', 'CDKN2B', 'CDKN1A', 'CDKN1B', 'CDK4', 'CDK6', 'CCND1', 'CCND2', 'CCND3'],
    'Angiogenesis': ['VEGFA', 'VEGFR1', 'VEGFR2', 'HIF1A', 'EPAS1', 'VHL', 'ANGPT1', 'ANGPT2']
}

# GBM-relevant pathway keywords (19 canonical GBM pathways)
GBM_PATHWAY_KEYWORDS = [
    # 1. RTK–RAS–MAPK signaling
    'RTK', 'receptor tyrosine kinase', 'EGFR', 'PDGFRA', 'ERBB',
    'RAS', 'MAPK', 'ERK', 'MEK', 'RAF', 'BRAF',
    # 2. PI3K–AKT–mTOR pathway
    'PI3K', 'AKT', 'mTOR', 'PTEN', 'PIK3',
    # 3. p53 signaling
    'p53', 'TP53', 'MDM2', 'MDM4',
    # 4. RB / G1-S cell-cycle axis
    'RB', 'retinoblastoma', 'CDK', 'CDKN', 'cell cycle', 'G1/S', 'cyclin',
    # 5. Telomere maintenance
    'TERT', 'telomere', 'ATRX', 'DAXX', 'telomerase',
    # 6. Chromatin remodeling
    'chromatin', 'SWI/SNF', 'ARID', 'SMARCA', 'histone', 'epigenetic',
    'H3K27', 'methylation', 'acetylation', 'SETD', 'KDM', 'KMT',
    # 7-19. Other GBM pathways
    'IDH', 'MGMT', 'DNA repair', 'NF-kB', 'JAK', 'STAT',
    'Notch', 'Wnt', 'hedgehog', 'apoptosis', 'autophagy',
    'hypoxia', 'HIF', 'angiogenesis', 'VEGF', 'immune', 'PD-L1',
    'invasion', 'EMT', 'mesenchymal', 'glycolysis', 'metabolism',
    'glioma', 'glioblastoma', 'GBM'
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_gene_data():
    """Load gene data - either from recurrence CSV or simple gene list"""
    print("Loading gene data...")

    # Load gene SV counts from Step 3 results
    if not os.path.exists(GENE_RECURRENCE_FILE):
        print(f"Error: Gene file not found: {GENE_RECURRENCE_FILE}")
        print("Please run 03_build_matrix_and_analyze.py first")
        sys.exit(1)

    print(f"  Loading from: {GENE_RECURRENCE_FILE}")
    df = pd.read_csv(GENE_RECURRENCE_FILE)

    # Rename columns to match expected format
    if 'n_samples' in df.columns:
        df = df.rename(columns={'n_samples': 'num_samples', 'n_svs': 'total_svs'})

    print(f"  ✓ Loaded {len(df)} genes with SV statistics")
    return df, 'csv'


def filter_genes(gene_df, min_recurrence=3, top_n=100, data_type='csv'):
    """
    Filter genes for enrichment analysis

    Criteria:
    1. Present in at least min_recurrence samples
    2. Take top N by number of samples affected
    """
    print(f"\nFiltering genes...")
    print(f"  Min recurrence: {min_recurrence} samples")
    print(f"  Top N: {top_n} genes")

    # Filter by recurrence
    filtered = gene_df[gene_df['num_samples'] >= min_recurrence].copy()
    print(f"  After recurrence filter: {len(filtered)} genes")

    # Sort by number of samples and take top N
    filtered = filtered.sort_values('num_samples', ascending=False).head(top_n)
    print(f"  After top N filter: {len(filtered)} genes")

    return filtered


def run_gprofiler_enrichment(gene_list):
    """
    Run enrichment analysis using g:Profiler API

    Returns enrichment results for:
    - GO:BP (Biological Process)
    - GO:MF (Molecular Function)
    - GO:CC (Cellular Component)
    - KEGG (Pathways)
    - REAC (Reactome)
    """
    if not GPROFILER_AVAILABLE:
        print("  Skipping g:Profiler (not installed)")
        return None

    print("\nRunning g:Profiler enrichment analysis...")
    print(f"  Query: {len(gene_list)} genes")

    try:
        gp = GProfiler(return_dataframe=True)

        result = gp.profile(
            organism='hsapiens',
            query=gene_list,
            sources=['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'TF', 'HP'],
            user_threshold=PVALUE_THRESHOLD,
            significance_threshold_method='fdr',
            background=None,  # Use all known genes as background
            no_evidences=False
        )

        if result is not None and len(result) > 0:
            print(f"  ✓ Found {len(result)} significant terms")
            return result
        else:
            print("  No significant enrichments found")
            return None

    except Exception as e:
        print(f"  Error running g:Profiler: {e}")
        return None


def run_enrichr_analysis(gene_list):
    """
    Run Enrichr analysis with GBM-relevant databases

    Returns enrichment results from multiple pathway databases
    """
    if not USE_ENRICHR:
        print("  Skipping Enrichr (disabled in config)")
        return None

    print("\nRunning Enrichr analysis...")
    print(f"  Query: {len(gene_list)} genes")

    # GBM and cancer-relevant pathway databases
    databases = [
        'KEGG_2021_Human',
        'GO_Biological_Process_2021',
        'Reactome_2022',
        'WikiPathway_2021_Human',
        'MSigDB_Hallmark_2020',
    ]

    # Enrichr API endpoints
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    ENRICHR_QUERY_URL = 'https://maayanlab.cloud/Enrichr/enrich'

    # Submit gene list
    genes_str = '\n'.join(gene_list)
    payload = {
        'list': (None, genes_str),
        'description': (None, 'GBM_SV_Affected_Genes')
    }

    try:
        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            print(f"  Error submitting genes: {response.status_code}")
            return None

        data = json.loads(response.text)
        user_list_id = data['userListId']
        print(f"  Gene list submitted (ID: {user_list_id})")

        # Query each database
        all_results = {}
        for db in databases:
            query_url = f"{ENRICHR_QUERY_URL}?userListId={user_list_id}&backgroundType={db}"
            response = requests.get(query_url)

            if response.ok:
                results = json.loads(response.text)
                if db in results and len(results[db]) > 0:
                    # Filter by p-value
                    filtered_results = [r for r in results[db] if r[6] < PVALUE_THRESHOLD]  # r[6] is Adjusted P-value
                    if filtered_results:
                        all_results[db] = filtered_results
                        print(f"  ✓ {db}: {len(filtered_results)} significant pathways")

        if all_results:
            print(f"  ✓ Total: {sum(len(v) for v in all_results.values())} pathways from {len(all_results)} databases")
            return all_results
        else:
            print("  No significant enrichments found")
            return None

    except Exception as e:
        print(f"  Error running Enrichr: {e}")
        return None


def identify_gbm_pathways(enrichr_results):
    """
    Identify and highlight GBM-relevant pathways from Enrichr results
    """
    if not enrichr_results:
        return pd.DataFrame()

    gbm_pathways = []

    for db_name, pathways in enrichr_results.items():
        for pathway_data in pathways:
            term = pathway_data[1].lower()  # Pathway term

            # Check if pathway is GBM-relevant
            is_gbm_relevant = any(keyword.lower() in term for keyword in GBM_PATHWAY_KEYWORDS)

            if is_gbm_relevant:
                gbm_pathways.append({
                    'Database': db_name,
                    'Term': pathway_data[1],
                    'P-value': pathway_data[2],
                    'Adjusted P-value': pathway_data[6],
                    'Combined Score': pathway_data[4],
                    'Genes': pathway_data[5]
                })

    gbm_df = pd.DataFrame(gbm_pathways)
    if len(gbm_df) > 0:
        gbm_df = gbm_df.sort_values('Combined Score', ascending=False)
        print(f"\n  ✓ Identified {len(gbm_df)} GBM-relevant pathways")

    return gbm_df


def analyze_core_gbm_pathways(gene_list):
    """
    Analyze overlap with core GBM signaling pathways
    """
    print("\nAnalyzing core GBM signaling pathways...")

    gene_set = set(gene_list)
    results = []

    for pathway, pathway_genes in CORE_GBM_PATHWAYS.items():
        pathway_set = set(pathway_genes)
        overlap = gene_set.intersection(pathway_set)

        if overlap:
            results.append({
                'pathway': pathway,
                'pathway_size': len(pathway_set),
                'overlap_count': len(overlap),
                'overlapping_genes': ';'.join(sorted(overlap)),
                'enrichment_ratio': len(overlap) / len(pathway_set)
            })

    if results:
        core_df = pd.DataFrame(results)
        core_df = core_df.sort_values('enrichment_ratio', ascending=False)
        print(f"  ✓ Found overlap with {len(core_df)} core GBM pathways")
        return core_df
    else:
        print("  No overlap with core GBM pathways")
        return pd.DataFrame()


def analyze_gbm_pathways(gene_list):
    """
    Analyze overlap with known GBM pathways
    """
    print("\nAnalyzing GBM-specific pathway overlap...")

    gene_set = set(gene_list)
    results = []

    for pathway, pathway_genes in EXPECTED_GBM_PATHWAYS.items():
        pathway_set = set(pathway_genes)
        overlap = gene_set.intersection(pathway_set)

        if overlap:
            result = {
                'pathway': pathway,
                'pathway_size': len(pathway_set),
                'overlap_size': len(overlap),
                'overlap_genes': ', '.join(sorted(overlap)),
                'enrichment_ratio': len(overlap) / len(pathway_set)
            }
            results.append(result)

            print(f"  {pathway}: {len(overlap)}/{len(pathway_set)} genes")
            print(f"    Genes: {', '.join(sorted(overlap))}")

    return pd.DataFrame(results)


def create_enrichment_plots(enrichment_df, output_dir):
    """
    Create visualization of enrichment results
    """
    if enrichment_df is None or len(enrichment_df) == 0:
        print("  No enrichment data to plot")
        return

    print("\nGenerating enrichment visualizations...")

    # Plot 1: Top enriched terms by source
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    sources = ['GO:BP', 'KEGG', 'REAC', 'GO:MF']

    for idx, source in enumerate(sources):
        ax = axes[idx // 2, idx % 2]

        source_data = enrichment_df[enrichment_df['source'] == source].copy()

        if len(source_data) == 0:
            ax.text(0.5, 0.5, f'No significant\n{source} terms',
                   ha='center', va='center', fontsize=12)
            ax.set_title(f'{source} Enrichment', fontsize=14, fontweight='bold')
            ax.axis('off')
            continue

        # Take top 10
        source_data = source_data.nsmallest(10, 'p_value')

        # Calculate -log10(p-value)
        source_data['-log10(p)'] = -np.log10(source_data['p_value'])

        # Plot horizontal bar chart
        y_pos = np.arange(len(source_data))
        ax.barh(y_pos, source_data['-log10(p)'], color='steelblue', edgecolor='black')

        # Truncate long names
        labels = [term[:50] + '...' if len(term) > 50 else term
                 for term in source_data['name']]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=9)

        ax.set_xlabel('-log10(p-value)', fontsize=11)
        ax.set_title(f'{source} Enrichment (Top 10)', fontsize=14, fontweight='bold')
        ax.axvline(-np.log10(PVALUE_THRESHOLD), color='red', linestyle='--',
                  alpha=0.5, label=f'p={PVALUE_THRESHOLD}')
        ax.legend()
        ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/enrichment_overview.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: enrichment_overview.png")

    # Plot 2: Enrichment bubble plot with labels
    if len(enrichment_df) > 0:
        fig, ax = plt.subplots(figsize=(16, 12))

        # Take top 20 terms for better readability with labels
        top_terms = enrichment_df.nsmallest(20, 'p_value').copy()
        top_terms['-log10(p)'] = -np.log10(top_terms['p_value'])

        # Color by source
        source_colors = {
            'GO:BP': '#e74c3c',
            'GO:MF': '#3498db',
            'GO:CC': '#2ecc71',
            'KEGG': '#f39c12',
            'REAC': '#9b59b6',
            'TF': '#1abc9c',
            'HP': '#34495e'
        }

        colors = [source_colors.get(s, 'gray') for s in top_terms['source']]

        scatter = ax.scatter(
            top_terms['intersection_size'],
            top_terms['-log10(p)'],
            s=top_terms['intersection_size'] * 50,
            c=colors,
            alpha=0.6,
            edgecolors='black',
            linewidth=0.5
        )

        # Add term labels to bubbles
        for idx, row in top_terms.iterrows():
            # Truncate long term names
            term_name = row['name'] if len(row['name']) <= 35 else row['name'][:32] + '...'

            # Get genes for this term (first 5 genes for annotation)
            gene_label = ''
            if 'intersections' in top_terms.columns:
                try:
                    if pd.notna(row['intersections']):
                        genes = str(row['intersections']).split(',')[:5]
                        gene_label = ', '.join(genes)
                        if len(str(row['intersections']).split(',')) > 5:
                            gene_label += '...'
                except:
                    gene_label = ''

            # Position label slightly offset from bubble
            x_offset = 0.3
            y_offset = 0.5

            # Add term name
            ax.annotate(
                term_name,
                xy=(row['intersection_size'], row['-log10(p)']),
                xytext=(row['intersection_size'] + x_offset, row['-log10(p)'] + y_offset),
                fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='gray'),
                ha='left'
            )

            # Add gene label below term name (smaller font)
            if gene_label:
                ax.annotate(
                    f"Genes: {gene_label}",
                    xy=(row['intersection_size'], row['-log10(p)']),
                    xytext=(row['intersection_size'] + x_offset, row['-log10(p)'] + y_offset - 1.5),
                    fontsize=6,
                    style='italic',
                    color='#555555',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow', alpha=0.6, edgecolor='none'),
                    ha='left'
                )

        ax.set_xlabel('Gene Count', fontsize=12, fontweight='bold')
        ax.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')
        ax.set_title('Top 20 Enriched Terms with Gene Labels', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Add some padding to axes for labels
        ax.margins(x=0.15, y=0.1)

        # Legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=source_colors[s], label=s)
                          for s in source_colors.keys()
                          if s in top_terms['source'].values]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

        plt.tight_layout()
        plt.savefig(f"{output_dir}/enrichment_bubble.png", dpi=300, bbox_inches='tight')
        plt.close()

        print(f"  ✓ Saved: enrichment_bubble.png (with term and gene labels)")


def create_gbm_pathway_plot(gbm_pathway_df, output_dir):
    """
    Visualize GBM pathway overlap
    """
    if gbm_pathway_df is None or len(gbm_pathway_df) == 0:
        return

    print("\nGenerating GBM pathway visualization...")

    fig, ax = plt.subplots(figsize=(10, 6))

    # Bar plot of enrichment ratios
    y_pos = np.arange(len(gbm_pathway_df))
    colors = ['#e74c3c' if r > 0.5 else '#3498db' for r in gbm_pathway_df['enrichment_ratio']]

    bars = ax.barh(y_pos, gbm_pathway_df['enrichment_ratio'], color=colors,
                   edgecolor='black', linewidth=1.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(gbm_pathway_df['pathway'])
    ax.set_xlabel('Enrichment Ratio (Overlap / Pathway Size)', fontsize=12)
    ax.set_title('GBM Pathway Overlap', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)

    # Add gene counts as text
    for i, (idx, row) in enumerate(gbm_pathway_df.iterrows()):
        ax.text(row['enrichment_ratio'] + 0.02, i,
               f"{row['overlap_size']}/{row['pathway_size']}",
               va='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/gbm_pathway_overlap.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: gbm_pathway_overlap.png")


def create_core_gbm_pathway_plot(core_gbm_df, output_dir):
    """
    Visualize core GBM signaling pathway overlap
    """
    if core_gbm_df is None or len(core_gbm_df) == 0:
        return

    print("\nGenerating core GBM pathway visualization...")

    fig, ax = plt.subplots(figsize=(12, 8))

    # Bar plot of enrichment ratios
    y_pos = np.arange(len(core_gbm_df))
    colors = ['#e74c3c' if r > 0.3 else '#3498db' for r in core_gbm_df['enrichment_ratio']]

    bars = ax.barh(y_pos, core_gbm_df['enrichment_ratio'], color=colors,
                   edgecolor='black', linewidth=1.5, alpha=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(core_gbm_df['pathway'], fontsize=11, fontweight='bold')
    ax.set_xlabel('Enrichment Ratio (Genes with SVs / Total Pathway Genes)',
                  fontsize=12, fontweight='bold')
    ax.set_title('Core GBM Signaling Pathway Enrichment',
                 fontsize=14, fontweight='bold', pad=20)
    ax.grid(axis='x', alpha=0.3)

    # Add gene counts and enrichment ratios as text
    for i, (idx, row) in enumerate(core_gbm_df.iterrows()):
        ax.text(row['enrichment_ratio'] + 0.02, i,
               f"{row['overlap_count']}/{row['pathway_size']} ({row['enrichment_ratio']:.1%})",
               va='center', fontsize=10, fontweight='bold')

    # Add gene names as annotations
    for i, (idx, row) in enumerate(core_gbm_df.iterrows()):
        genes = row['overlapping_genes'].split(';')
        gene_text = ', '.join(genes[:5])  # First 5 genes
        if len(genes) > 5:
            gene_text += '...'
        ax.text(0.01, i, gene_text, va='center', ha='left',
               fontsize=8, style='italic', color='white',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.6))

    plt.tight_layout()
    plt.savefig(f"{output_dir}/core_gbm_pathways.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: core_gbm_pathways.png")


def export_for_external_tools(gene_list, filtered_gene_df, output_dir):
    """
    Export gene lists in formats suitable for external enrichment tools
    """
    print("\nExporting gene lists for external tools...")

    # Simple text list (for DAVID, Enrichr, etc.)
    with open(f"{output_dir}/gene_list.txt", 'w') as f:
        for gene in gene_list:
            f.write(f"{gene}\n")
    print(f"  ✓ Saved: gene_list.txt (for DAVID/Enrichr)")

    # GMT format (for GSEA)
    with open(f"{output_dir}/gene_list.gmt", 'w') as f:
        f.write(f"SV_AFFECTED_GENES\tSV_affected\t" + "\t".join(gene_list) + "\n")
    print(f"  ✓ Saved: gene_list.gmt (for GSEA)")

    # Ranked list (by recurrence) - use the already filtered dataframe
    gene_df = filtered_gene_df[filtered_gene_df['gene'].isin(gene_list)].copy()

    # Make sure we have the right column name
    rank_col = 'num_samples' if 'num_samples' in gene_df.columns else 'n_samples'

    gene_df = gene_df.sort_values(rank_col, ascending=False)

    with open(f"{output_dir}/gene_list_ranked.rnk", 'w') as f:
        for _, row in gene_df.iterrows():
            f.write(f"{row['gene']}\t{row[rank_col]}\n")
    print(f"  ✓ Saved: gene_list_ranked.rnk (ranked by recurrence)")


def generate_enrichment_report(enrichment_df, gbm_pathway_df, gene_list, output_dir):
    """
    Generate comprehensive text report
    """
    print("\nGenerating enrichment report...")

    report_file = f"{output_dir}/enrichment_report.txt"

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PATHWAY ENRICHMENT ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")

        # Summary
        f.write("SUMMARY\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total genes analyzed: {len(gene_list)}\n")

        if enrichment_df is not None:
            f.write(f"Significant terms (p < {PVALUE_THRESHOLD}): {len(enrichment_df)}\n")

            # By source
            f.write("\nTerms by source:\n")
            for source, count in enrichment_df['source'].value_counts().items():
                f.write(f"  {source}: {count}\n")
        else:
            f.write("No significant enrichments found via g:Profiler\n")

        # GBM pathways
        f.write("\n\nGBM-SPECIFIC PATHWAY OVERLAP\n")
        f.write("-" * 80 + "\n")

        if gbm_pathway_df is not None and len(gbm_pathway_df) > 0:
            for _, row in gbm_pathway_df.iterrows():
                f.write(f"\n{row['pathway']}:\n")
                f.write(f"  Overlap: {row['overlap_size']}/{row['pathway_size']} ")
                f.write(f"({row['enrichment_ratio']:.1%})\n")
                f.write(f"  Genes: {row['overlap_genes']}\n")
        else:
            f.write("No overlap with predefined GBM pathways\n")

        # Top enriched terms
        if enrichment_df is not None and len(enrichment_df) > 0:
            f.write("\n\nTOP 20 ENRICHED TERMS (ALL SOURCES)\n")
            f.write("-" * 80 + "\n\n")

            top20 = enrichment_df.nsmallest(20, 'p_value')

            for idx, row in top20.iterrows():
                f.write(f"{row['source']}: {row['name']}\n")
                f.write(f"  p-value: {row['p_value']:.2e}\n")
                if 'intersection_size' in row:
                    f.write(f"  Genes: {row['intersection_size']}\n")
                if 'intersections' in row:
                    genes = str(row['intersections'])[:100]
                    f.write(f"  Gene list: {genes}...\n")
                f.write("\n")

        # Input genes
        f.write("\n\nINPUT GENE LIST\n")
        f.write("-" * 80 + "\n")
        f.write(", ".join(gene_list))
        f.write("\n\n")

        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

    print(f"  ✓ Saved: enrichment_report.txt")


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    print("=" * 80)
    print("PATHWAY ENRICHMENT ANALYSIS")
    print("=" * 80)
    print()

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load gene data
    gene_df, data_type = load_gene_data()

    # Filter genes for enrichment
    filtered_genes = filter_genes(
        gene_df,
        min_recurrence=MIN_GENE_RECURRENCE,
        top_n=TOP_N_GENES,
        data_type=data_type
    )

    # Extract gene list
    gene_list = filtered_genes['gene'].tolist()

    print(f"\nGene list for enrichment ({len(gene_list)} genes):")
    print(f"  Top 10: {', '.join(gene_list[:10])}")

    # Save filtered gene table
    filtered_genes.to_csv(f"{OUTPUT_DIR}/input_genes.csv", index=False)
    print(f"\n✓ Saved input genes: input_genes.csv")

    # Run g:Profiler enrichment
    enrichment_df = None
    if GPROFILER_AVAILABLE:
        enrichment_df = run_gprofiler_enrichment(gene_list)

        if enrichment_df is not None:
            enrichment_df.to_csv(f"{OUTPUT_DIR}/gprofiler_results.csv", index=False)
            print(f"  ✓ Saved: gprofiler_results.csv")

    # Run Enrichr enrichment (NEW)
    enrichr_results = None
    enrichr_gbm_df = pd.DataFrame()
    if USE_ENRICHR:
        enrichr_results = run_enrichr_analysis(gene_list)

        if enrichr_results:
            # Save all Enrichr results
            all_enrichr = []
            for db_name, pathways in enrichr_results.items():
                for pathway_data in pathways:
                    all_enrichr.append({
                        'Database': db_name,
                        'Term': pathway_data[1],
                        'P-value': pathway_data[2],
                        'Z-score': pathway_data[3],
                        'Combined Score': pathway_data[4],
                        'Genes': pathway_data[5],
                        'Adjusted P-value': pathway_data[6]
                    })

            enrichr_df = pd.DataFrame(all_enrichr)
            enrichr_df.to_csv(f"{OUTPUT_DIR}/enrichr_results.csv", index=False)
            print(f"  ✓ Saved: enrichr_results.csv")

            # Identify GBM-relevant pathways from Enrichr
            enrichr_gbm_df = identify_gbm_pathways(enrichr_results)
            if len(enrichr_gbm_df) > 0:
                enrichr_gbm_df.to_csv(f"{OUTPUT_DIR}/gbm_relevant_pathways_enrichr.csv", index=False)
                print(f"  ✓ Saved: gbm_relevant_pathways_enrichr.csv")

    # Analyze core GBM signaling pathways (NEW)
    core_gbm_df = analyze_core_gbm_pathways(gene_list)
    if len(core_gbm_df) > 0:
        core_gbm_df.to_csv(f"{OUTPUT_DIR}/core_gbm_pathways.csv", index=False)
        print(f"  ✓ Saved: core_gbm_pathways.csv")

    # Analyze GBM-specific pathways (original method)
    gbm_pathway_df = analyze_gbm_pathways(gene_list)
    if len(gbm_pathway_df) > 0:
        gbm_pathway_df.to_csv(f"{OUTPUT_DIR}/gbm_pathway_overlap.csv", index=False)
        print(f"  ✓ Saved: gbm_pathway_overlap.csv")

    # Create visualizations
    if enrichment_df is not None:
        create_enrichment_plots(enrichment_df, OUTPUT_DIR)

    if len(gbm_pathway_df) > 0:
        create_gbm_pathway_plot(gbm_pathway_df, OUTPUT_DIR)

    # Create core GBM pathway visualization (NEW)
    if len(core_gbm_df) > 0:
        create_core_gbm_pathway_plot(core_gbm_df, OUTPUT_DIR)

    # Export for external tools
    export_for_external_tools(gene_list, filtered_genes, OUTPUT_DIR)

    # Generate report
    generate_enrichment_report(enrichment_df, gbm_pathway_df, gene_list, OUTPUT_DIR)

    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\nResults saved to: {OUTPUT_DIR}/")

    print("\n📊 KEY OUTPUTS:")
    print("\n  Enrichment Results:")
    if enrichment_df is not None:
        print(f"    ✓ gprofiler_results.csv         : g:Profiler enrichment ({len(enrichment_df)} terms)")
    if enrichr_results:
        total_enrichr = sum(len(v) for v in enrichr_results.values())
        print(f"    ✓ enrichr_results.csv           : Enrichr enrichment ({total_enrichr} pathways, 5 databases)")
        if len(enrichr_gbm_df) > 0:
            print(f"    ✓ gbm_relevant_pathways_enrichr.csv : GBM-specific from Enrichr ({len(enrichr_gbm_df)} pathways)")

    print("\n  GBM Pathway Analysis:")
    if len(core_gbm_df) > 0:
        print(f"    ✓ core_gbm_pathways.csv         : Core GBM signaling pathways ({len(core_gbm_df)} pathways)")
        print(f"    ✓ core_gbm_pathways.png         : Core GBM pathway visualization")
    if len(gbm_pathway_df) > 0:
        print(f"    ✓ gbm_pathway_overlap.csv       : GBM pathway overlap ({len(gbm_pathway_df)} pathways)")
        print(f"    ✓ gbm_pathway_overlap.png       : GBM pathway bar plot")

    print("\n  Visualizations:")
    print("    ✓ enrichment_overview.png       : Top enriched terms by source")
    print("    ✓ enrichment_bubble.png         : Bubble plot with gene labels (top 20)")

    print("\n  Export Files:")
    print("    ✓ gene_list.txt                 : For DAVID/Enrichr")
    print("    ✓ gene_list.gmt                 : For GSEA")
    print("    ✓ gene_list_ranked.rnk          : Ranked by SV recurrence")
    print("    ✓ enrichment_report.txt         : Comprehensive text report")

    print("\n" + "=" * 80)

    # Usage instructions
    print("\n📖 NEXT STEPS:")
    print("1. Review enrichment_report.txt for summary")
    print("2. Check visualizations (*.png files)")
    print("3. Examine core GBM pathway enrichment:")
    print(f"   - Expected high enrichment: RTK/EGFR, PI3K-AKT-mTOR, TP53, RB/Cell Cycle")
    print("4. Use gene_list.txt with external tools:")
    print("   - DAVID: https://david.ncifcrf.gov/")
    print("   - Enrichr: https://maayanlab.cloud/Enrichr/")
    print("   - g:Profiler: https://biit.cs.ut.ee/gprofiler/")
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
