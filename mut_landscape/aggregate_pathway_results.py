#!/usr/bin/env python3
"""
Aggregate Per-Sample Pathway Results
Shows which genes are associated with each pathway across all samples
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict
import ast

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)


def aggregate_pathway_genes(pathway_file):
    """
    Aggregate pathways across samples and collect all genes
    """
    print(f"\nLoading pathway data from: {pathway_file}")
    df = pd.read_csv(pathway_file, sep='\t')

    print(f"  Total pathway enrichments: {len(df)}")
    print(f"  Unique pathways: {df['Term'].nunique()}")
    print(f"  Samples analyzed: {df['Sample_ID'].nunique()}")

    # Aggregate by pathway term
    pathway_aggregation = defaultdict(lambda: {
        'samples': set(),
        'genes': set(),
        'databases': set(),
        'p_values': [],
        'combined_scores': []
    })

    for idx, row in df.iterrows():
        term = row['Term']
        sample = row['Sample_ID']
        database = row['Database']
        p_val = row['Adjusted P-value']
        score = row['Combined Score']

        # Parse genes (they're stored as string representation of list)
        genes_str = row['Genes']
        if isinstance(genes_str, str):
            try:
                # Try to parse as Python list
                genes = ast.literal_eval(genes_str)
            except:
                # If that fails, try splitting by common separators
                genes = [g.strip() for g in genes_str.replace('[', '').replace(']', '').replace("'", "").split(',')]
        else:
            genes = []

        pathway_aggregation[term]['samples'].add(sample)
        pathway_aggregation[term]['genes'].update(genes)
        pathway_aggregation[term]['databases'].add(database)
        pathway_aggregation[term]['p_values'].append(p_val)
        pathway_aggregation[term]['combined_scores'].append(score)

    # Convert to dataframe
    aggregated_data = []
    for term, data in pathway_aggregation.items():
        aggregated_data.append({
            'Pathway': term,
            'N_Samples': len(data['samples']),
            'Samples': ';'.join(sorted(data['samples'])),
            'N_Genes': len(data['genes']),
            'Genes': ';'.join(sorted(data['genes'])),
            'Databases': ';'.join(sorted(data['databases'])),
            'Mean_Adj_P_value': np.mean(data['p_values']),
            'Min_Adj_P_value': np.min(data['p_values']),
            'Mean_Combined_Score': np.mean(data['combined_scores']),
            'Max_Combined_Score': np.max(data['combined_scores'])
        })

    agg_df = pd.DataFrame(aggregated_data)
    agg_df = agg_df.sort_values('N_Samples', ascending=False)

    return agg_df


def create_pathway_gene_visualizations(agg_df, output_dir):
    """Create visualizations showing pathway-gene associations"""

    print(f"\n{'='*80}")
    print("CREATING PATHWAY-GENE ASSOCIATION VISUALIZATIONS")
    print(f"{'='*80}\n")

    # 1. Top pathways with gene counts
    print("[1/4] Top recurrent pathways with gene associations...")
    top_pathways = agg_df.head(30)

    fig, ax = plt.subplots(figsize=(16, 12))

    y_pos = np.arange(len(top_pathways))

    # Create bars colored by number of samples
    colors = plt.cm.RdYlGn_r(top_pathways['N_Samples'].values / top_pathways['N_Samples'].max())
    bars = ax.barh(y_pos, top_pathways['N_Samples'].values,
                   color=colors, edgecolor='black', alpha=0.8)

    # Pathway labels
    pathway_labels = [term[:70] + '...' if len(term) > 70 else term
                     for term in top_pathways['Pathway'].values]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pathway_labels, fontsize=9)
    ax.set_xlabel('Number of Samples with Enrichment', fontsize=12, fontweight='bold')
    ax.set_title('Top 30 Recurrent GBM Pathways\n(with gene associations)',
                fontsize=14, fontweight='bold', pad=20, color='darkred')
    ax.grid(axis='x', alpha=0.3)

    # Add sample count and gene count labels
    for i, (n_samples, n_genes) in enumerate(zip(top_pathways['N_Samples'].values,
                                                   top_pathways['N_Genes'].values)):
        ax.text(n_samples, i, f' {n_samples} samples\n({n_genes} genes)',
               va='center', fontsize=8, fontweight='bold')

    plt.tight_layout()
    output_path = output_dir / 'top_pathways_with_genes.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"   Saved: {output_path}")
    plt.close()

    # 2. Gene count vs sample count scatter
    print("[2/4] Gene-sample relationship scatter plot...")
    fig, ax = plt.subplots(figsize=(12, 8))

    scatter = ax.scatter(agg_df['N_Samples'], agg_df['N_Genes'],
                        c=agg_df['Mean_Combined_Score'], cmap='viridis',
                        s=100, alpha=0.6, edgecolors='black')

    ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes in Pathway', fontsize=12, fontweight='bold')
    ax.set_title('Pathway Recurrence vs Gene Count',
                fontsize=14, fontweight='bold', pad=20)
    ax.grid(alpha=0.3)

    # Color bar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Mean Combined Score', fontsize=11, fontweight='bold')

    # Annotate top pathways
    top_5 = agg_df.head(5)
    for idx, row in top_5.iterrows():
        pathway_label = row['Pathway'][:30] + '...' if len(row['Pathway']) > 30 else row['Pathway']
        ax.annotate(pathway_label,
                   xy=(row['N_Samples'], row['N_Genes']),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5))

    plt.tight_layout()
    output_path = output_dir / 'pathway_gene_sample_relationship.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"   Saved: {output_path}")
    plt.close()

    # 3. Top pathway-gene table
    print("[3/4] Creating pathway-gene association table...")
    fig, ax = plt.subplots(figsize=(18, 14))
    ax.axis('tight')
    ax.axis('off')

    top_15 = agg_df.head(15)
    table_data = []
    for idx, row in top_15.iterrows():
        pathway = row['Pathway'][:45] + '...' if len(row['Pathway']) > 45 else row['Pathway']
        genes = row['Genes'][:60] + '...' if len(row['Genes']) > 60 else row['Genes']
        table_data.append([
            pathway,
            str(row['N_Samples']),
            str(row['N_Genes']),
            genes
        ])

    table = ax.table(cellText=table_data,
                    colLabels=['Pathway', '#Samples', '#Genes', 'Associated Genes'],
                    cellLoc='left',
                    loc='center',
                    colWidths=[0.35, 0.08, 0.08, 0.49])

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2.5)

    # Style header
    for i in range(4):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(weight='bold', color='white', fontsize=10)

    # Alternate row colors
    for i in range(1, len(table_data) + 1):
        for j in range(4):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#ecf0f1')

    ax.set_title('Top 15 Recurrent Pathways and Associated Genes',
                fontsize=16, fontweight='bold', pad=30)

    plt.tight_layout()
    output_path = output_dir / 'pathway_gene_association_table.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"   Saved: {output_path}")
    plt.close()

    # 4. Most frequent genes across pathways
    print("[4/4] Most frequent genes across pathways...")
    all_genes = []
    for genes_str in agg_df['Genes']:
        genes = genes_str.split(';')
        all_genes.extend([g for g in genes if g])

    gene_counts = Counter(all_genes).most_common(40)

    if gene_counts:
        fig, ax = plt.subplots(figsize=(14, 12))
        genes, counts = zip(*gene_counts)

        bars = ax.barh(range(len(genes)), counts,
                      color='#e74c3c', edgecolor='black', alpha=0.8)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=10, fontweight='bold')
        ax.set_xlabel('Number of Enriched Pathways', fontsize=12, fontweight='bold')
        ax.set_title('Top 40 Genes Most Frequently in Enriched Pathways\n(across all samples)',
                    fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='x', alpha=0.3)

        for i, count in enumerate(counts):
            ax.text(count, i, f' {count}', va='center', fontsize=9, fontweight='bold')

        plt.tight_layout()
        output_path = output_dir / 'most_frequent_genes_in_pathways.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"   Saved: {output_path}")
        plt.close()

    print(f"\n{'='*80}")
    print("All visualizations generated successfully!")
    print(f"{'='*80}")


def main():
    """Main aggregation workflow"""
    print("="*80)
    print("AGGREGATE PATHWAY-GENE ASSOCIATIONS")
    print("="*80)

    # Paths
    pathway_file = "/home/chbope/extension/script/mut_landscape/results/gbm_pathway_analysis_v2_persample/gbm_pathways_per_sample.tsv"
    output_dir = Path("/home/chbope/extension/script/mut_landscape/results/gbm_pathway_analysis_v2_persample")

    # Aggregate pathways
    agg_df = aggregate_pathway_genes(pathway_file)

    print(f"\n{'='*80}")
    print("AGGREGATION SUMMARY")
    print(f"{'='*80}")
    print(f"Total unique pathways: {len(agg_df)}")
    print(f"Most recurrent pathway: {agg_df.iloc[0]['Pathway']}")
    print(f"  - Enriched in {agg_df.iloc[0]['N_Samples']} samples")
    print(f"  - Involves {agg_df.iloc[0]['N_Genes']} genes: {agg_df.iloc[0]['Genes'][:100]}...")
    print(f"\nPathway with most genes: {agg_df.iloc[agg_df['N_Genes'].argmax()]['Pathway']}")
    print(f"  - Involves {agg_df['N_Genes'].max()} genes")
    print(f"  - Enriched in {agg_df.iloc[agg_df['N_Genes'].argmax()]['N_Samples']} samples")

    # Save aggregated results
    output_file = output_dir / 'aggregated_pathways_with_genes.tsv'
    agg_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved aggregated pathways: {output_file}")

    # Create visualizations
    create_pathway_gene_visualizations(agg_df, output_dir)

    print(f"\n{'='*80}")
    print("AGGREGATION COMPLETE!")
    print(f"{'='*80}")
    print(f"\nResults: {output_dir}")
    print(f"  - aggregated_pathways_with_genes.tsv: Pathway-gene associations")
    print(f"  - 4 new visualization PNG files")
    print("="*80)


if __name__ == "__main__":
    main()
