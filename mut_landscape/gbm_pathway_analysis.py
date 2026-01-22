#!/usr/bin/env python3
"""
GBM-Focused Pathway Enrichment Analysis
Analyzes pathway enrichment for glioblastoma-related genes from somatic variants
Includes GBM-specific pathway databases and visualization
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from collections import Counter
import requests
import json

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

class GBMPathwayAnalyzer:
    """
    GBM-specific pathway enrichment analysis for mutated genes
    """

    def __init__(self, genes, background_genes=None):
        """
        Initialize GBM pathway analyzer

        Args:
            genes: List of genes to analyze (from VCF/mutations)
            background_genes: Optional background genes (e.g., cancer gene panel)
        """
        self.genes = list(set(genes))
        self.background_genes = background_genes
        self.enrichment_results = {}

        # GBM-relevant pathways to highlight (19 canonical GBM pathways)
        self.gbm_key_pathways = [
            # 1. RTK–RAS–MAPK signaling
            'RTK', 'receptor tyrosine kinase', 'EGFR', 'PDGFRA', 'ERBB',
            'RAS', 'MAPK', 'ERK', 'MEK', 'RAF', 'BRAF',
            # 2. PI3K–AKT–mTOR pathway
            'PI3K', 'AKT', 'mTOR', 'PTEN', 'PIK3',
            # 3. p53 signaling
            'p53', 'TP53', 'MDM2', 'MDM4',
            # 4. RB / G1-S cell-cycle axis
            'RB', 'retinoblastoma', 'CDK', 'CDKN', 'cell cycle', 'G1/S', 'cyclin',
            # 5. Telomere maintenance (TERT / ATRX)
            'TERT', 'telomere', 'ATRX', 'DAXX', 'telomerase',
            # 6. Chromatin remodeling (SWI/SNF, epigenetic regulators)
            'chromatin', 'SWI/SNF', 'ARID', 'SMARCA', 'histone', 'epigenetic',
            'H3K27', 'methylation', 'acetylation', 'SETD', 'KDM', 'KMT',
            # 7. IDH pathway (oncometabolite)
            'IDH', '2-hydroxyglutarate', 'oncometabolite',
            # 8. MGMT & TMZ resistance
            'MGMT', 'temozolomide', 'TMZ', 'O6-methylguanine', 'alkylating',
            # 9. DNA repair (HR, NHEJ, MMR, BER)
            'DNA repair', 'homologous recombination', 'NHEJ', 'mismatch repair',
            'MMR', 'base excision', 'BER', 'BRCA', 'ATM', 'ATR', 'CHEK',
            # 10. NF-κB signaling
            'NF-kB', 'NF-kappaB', 'NFkB', 'NFKB', 'RELA', 'RELB', 'IKK',
            # 11. JAK–STAT3 signaling
            'JAK', 'STAT', 'STAT3', 'cytokine', 'interleukin',
            # 12. Notch / Wnt / SHH stemness pathways
            'Notch', 'Wnt', 'beta-catenin', 'hedgehog', 'SHH', 'GLI',
            'stemness', 'stem cell', 'progenitor',
            # 13. Apoptosis & necroptosis
            'apoptosis', 'necroptosis', 'cell death', 'caspase', 'BCL',
            'BAX', 'RIPK', 'programmed cell death',
            # 14. Autophagy regulation
            'autophagy', 'ATG', 'mTORC1', 'ULK', 'LC3', 'BECN',
            # 15. Hypoxia / HIF signaling
            'hypoxia', 'HIF', 'oxygen', 'normoxia',
            # 16. Angiogenesis / VEGF
            'angiogenesis', 'VEGF', 'vascular', 'endothelial',
            # 17. Immune evasion (PD-L1, CD47, TAM/M2)
            'PD-L1', 'PD-1', 'PDCD1', 'CD274', 'CD47', 'checkpoint',
            'immune evasion', 'immunosuppression', 'macrophage', 'TAM',
            'M2', 'T cell', 'immune',
            # 18. Invasion / EMT / mesenchymal transition
            'invasion', 'migration', 'EMT', 'mesenchymal', 'epithelial',
            'metastasis', 'MMP', 'matrix metalloproteinase',
            # 19. Metabolic reprogramming (glycolysis, glutamine, mitochondria)
            'glycolysis', 'Warburg', 'glucose', 'glutamine', 'glutaminolysis',
            'mitochondria', 'OXPHOS', 'lactate', 'metabolism',
            # General glioma terms
            'glioma', 'glioblastoma', 'astrocytoma', 'GBM'
        ]

        print(f"Initialized GBM pathway analyzer with {len(self.genes)} unique genes")

    def run_enrichr_analysis(self):
        """
        Run Enrichr analysis with GBM-relevant databases
        """
        # GBM and cancer-relevant pathway databases
        databases = [
            'KEGG_2021_Human',                    # General pathways
            'GO_Biological_Process_2021',         # GO biological processes
            'Reactome_2022',                      # Reactome pathways
            'WikiPathway_2021_Human',             # WikiPathways
            'MSigDB_Hallmark_2020',               # Hallmark gene sets
            'BioPlanet_2019',                     # BioPlanet pathways
            'OMIM_Disease',                       # Disease associations
            'DisGeNET',                           # Disease-gene associations
            'Jensen_DISEASES',                    # Text-mining disease associations
        ]

        print(f"\nRunning Enrichr analysis with {len(databases)} databases...")
        print("Focus: GBM/Glioma-relevant pathways")

        # Enrichr API endpoints
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
        ENRICHR_QUERY_URL = 'https://maayanlab.cloud/Enrichr/enrich'

        # Submit gene list
        genes_str = '\n'.join(self.genes)
        payload = {
            'list': (None, genes_str),
            'description': (None, 'GBM_Somatic_Mutations')
        }

        try:
            response = requests.post(ENRICHR_URL, files=payload)
            if not response.ok:
                print(f"Error submitting genes: {response.status_code}")
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
                        all_results[db] = results[db]
                        print(f"  ✓ {db}: {len(results[db])} pathways")
                    else:
                        print(f"  - {db}: No significant pathways")
                else:
                    print(f"  ✗ {db}: Query failed")

            self.enrichment_results = all_results
            return all_results

        except Exception as e:
            print(f"Error: {e}")
            return None

    def identify_gbm_pathways(self, results):
        """
        Identify and highlight GBM-relevant pathways
        """
        gbm_pathways = []

        for db_name, pathways in results.items():
            for pathway_data in pathways:
                term = pathway_data[1].lower()  # Pathway term

                # Check if pathway is GBM-relevant
                is_gbm_relevant = any(keyword.lower() in term for keyword in self.gbm_key_pathways)

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
            print(f"\n{'='*60}")
            print(f"IDENTIFIED {len(gbm_df)} GBM-RELEVANT PATHWAYS")
            print(f"{'='*60}")

        return gbm_df

    def parse_results(self, max_pathways=30, p_threshold=0.05):
        """Parse and filter enrichment results"""
        parsed_results = {}

        for db_name, pathways in self.enrichment_results.items():
            df = pd.DataFrame(pathways, columns=[
                'Rank', 'Term', 'P-value', 'Z-score', 'Combined Score',
                'Overlapping Genes', 'Adjusted P-value', 'Old P-value',
                'Old Adjusted P-value'
            ])

            # Filter by p-value
            df_filtered = df[df['Adjusted P-value'] < p_threshold].copy()
            df_filtered = df_filtered.head(max_pathways)

            if len(df_filtered) > 0:
                parsed_results[db_name] = df_filtered
                print(f"\n{db_name}: {len(df_filtered)} significant pathways")

        return parsed_results

    def create_gbm_visualizations(self, results, gbm_pathways, output_dir):
        """Create GBM-focused visualizations"""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

        print(f"\n{'='*60}")
        print("CREATING GBM-FOCUSED VISUALIZATIONS")
        print(f"{'='*60}\n")

        # 1. GBM-specific pathways
        if len(gbm_pathways) > 0:
            print("[1/6] GBM-relevant pathways...")
            self._plot_gbm_pathways(gbm_pathways, output_dir)

        # 2. Core GBM signaling pathways
        print("[2/6] Core GBM signaling pathways...")
        self._plot_core_gbm_pathways(results, output_dir)

        # 3. Top pathways across all databases
        print("[3/6] Top enriched pathways...")
        self._plot_top_pathways(results, output_dir)

        # 4. Pathway categories
        print("[4/6] Pathway category distribution...")
        self._plot_pathway_categories(results, output_dir)

        # 5. Database comparison
        print("[5/6] Database comparison...")
        self._plot_database_comparison(results, output_dir)

        # 6. Gene-pathway network
        print("[6/6] Gene-pathway network...")
        self._plot_gene_pathway_network(results, output_dir)

    def _plot_gbm_pathways(self, gbm_df, output_dir, top_n=30):
        """Plot GBM-specific pathways"""
        if len(gbm_df) == 0:
            return

        df_plot = gbm_df.head(top_n)

        fig, ax = plt.subplots(figsize=(14, 10))

        # Color by -log10(p-value)
        neg_log_p = -np.log10(df_plot['Adjusted P-value'].values)
        colormap = plt.cm.RdYlGn_r

        bars = ax.barh(range(len(df_plot)), df_plot['Combined Score'].values)

        for bar, log_p in zip(bars, neg_log_p):
            color_val = min(log_p / 10, 1)
            bar.set_color(colormap(color_val))
            bar.set_edgecolor('black')
            bar.set_linewidth(0.5)

        # Pathway labels
        pathway_labels = [term[:60] + '...' if len(term) > 60 else term
                         for term in df_plot['Term'].values]

        ax.set_yticks(range(len(df_plot)))
        ax.set_yticklabels(pathway_labels, fontsize=9)
        ax.set_xlabel('Enrichment Score', fontsize=12, fontweight='bold')
        ax.set_title('GBM-Relevant Enriched Pathways',
                    fontsize=14, fontweight='bold', pad=30,
                    color='darkred')
        ax.grid(axis='x', alpha=0.3)

        # Add score and p-value labels
        for i, (score, pval) in enumerate(zip(df_plot['Combined Score'].values,
                                              df_plot['Adjusted P-value'].values)):
            ax.text(score, i, f' {score:.1f}\n(p={pval:.1e})',
                   va='center', fontsize=7, fontweight='bold')

        plt.tight_layout()
        output_path = output_dir / 'gbm_relevant_pathways.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close()

    def _plot_core_gbm_pathways(self, results, output_dir):
        """Plot core GBM signaling pathways (RTK, PI3K-AKT-mTOR, RAS-MAPK, p53, RB)"""
        core_pathways = {
            'RTK/EGFR Signaling': ['RTK', 'EGFR', 'PDGFR', 'receptor tyrosine kinase', 'ErbB'],
            'PI3K-AKT-mTOR': ['PI3K', 'AKT', 'mTOR', 'PTEN'],
            'RAS-MAPK': ['RAS', 'MAPK', 'ERK', 'MEK', 'RAF'],
            'TP53 Pathway': ['p53', 'TP53', 'MDM', 'apoptosis'],
            'RB/Cell Cycle': ['RB', 'retinoblastoma', 'CDK', 'CDKN', 'cell cycle', 'G1/S'],
            'Angiogenesis': ['VEGF', 'angiogenesis', 'hypoxia', 'HIF']
        }

        pathway_scores = {cat: [] for cat in core_pathways.keys()}

        for db_name, df in results.items():
            for _, row in df.iterrows():
                term_lower = row['Term'].lower()
                for category, keywords in core_pathways.items():
                    if any(kw.lower() in term_lower for kw in keywords):
                        pathway_scores[category].append(row['Combined Score'])

        # Calculate average scores
        avg_scores = {cat: np.mean(scores) if scores else 0
                     for cat, scores in pathway_scores.items()}
        counts = {cat: len(scores) for cat, scores in pathway_scores.items()}

        # Plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Average enrichment scores
        categories = list(avg_scores.keys())
        scores = list(avg_scores.values())

        bars1 = ax1.barh(range(len(categories)), scores,
                        color='#e74c3c', edgecolor='black', alpha=0.8)
        ax1.set_yticks(range(len(categories)))
        ax1.set_yticklabels(categories, fontsize=10, fontweight='bold')
        ax1.set_xlabel('Average Enrichment Score', fontsize=11, fontweight='bold')
        ax1.set_title('Core GBM Pathway Enrichment',
                     fontsize=13, fontweight='bold', pad=15)
        ax1.grid(axis='x', alpha=0.3)

        for i, score in enumerate(scores):
            if score > 0:
                ax1.text(score, i, f' {score:.1f}', va='center',
                        fontsize=9, fontweight='bold')

        # Pathway counts
        count_vals = list(counts.values())
        bars2 = ax2.bar(range(len(categories)), count_vals,
                       color='#3498db', edgecolor='black', alpha=0.8)
        ax2.set_xticks(range(len(categories)))
        ax2.set_xticklabels(categories, rotation=45, ha='right', fontsize=9)
        ax2.set_ylabel('Number of Enriched Pathways', fontsize=11, fontweight='bold')
        ax2.set_title('GBM Pathway Category Counts',
                     fontsize=13, fontweight='bold', pad=15)
        ax2.grid(axis='y', alpha=0.3)

        for i, count in enumerate(count_vals):
            if count > 0:
                ax2.text(i, count, f'{count}', ha='center', va='bottom',
                        fontsize=10, fontweight='bold')

        plt.tight_layout()
        output_path = output_dir / 'gbm_core_pathways.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close()

    def _plot_top_pathways(self, results, output_dir, top_n=30):
        """Plot top pathways across all databases"""
        all_pathways = []
        for db_name, df in results.items():
            for _, row in df.iterrows():
                all_pathways.append({
                    'Database': db_name.replace('_', ' '),
                    'Pathway': row['Term'],
                    'Adjusted_P_value': row['Adjusted P-value'],
                    'Combined_Score': row['Combined Score'],
                })

        if not all_pathways:
            return

        pathway_df = pd.DataFrame(all_pathways)
        pathway_df = pathway_df.sort_values('Combined_Score', ascending=False).head(top_n)

        fig, ax = plt.subplots(figsize=(14, 10))

        y_pos = range(len(pathway_df))
        neg_log_p = -np.log10(pathway_df['Adjusted_P_value'].values)

        bars = ax.barh(y_pos, pathway_df['Combined_Score'].values)

        colormap = plt.cm.RdYlGn_r
        for bar, log_p in zip(bars, neg_log_p):
            color_val = min(log_p / 10, 1)
            bar.set_color(colormap(color_val))
            bar.set_edgecolor('black')
            bar.set_linewidth(0.5)

        pathway_labels = [term[:65] + '...' if len(term) > 65 else term
                         for term in pathway_df['Pathway'].values]

        ax.set_yticks(y_pos)
        ax.set_yticklabels(pathway_labels, fontsize=9)
        ax.set_xlabel('Enrichment Score', fontsize=12, fontweight='bold')
        ax.set_title(f'Top {top_n} Enriched Pathways (All Databases)',
                    fontsize=14, fontweight='bold', pad=30)
        ax.grid(axis='x', alpha=0.3)

        for i, score in enumerate(pathway_df['Combined_Score'].values):
            ax.text(score, i, f' {score:.1f}', va='center',
                   fontsize=8, fontweight='bold')

        plt.tight_layout()
        output_path = output_dir / 'top_pathways_all.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close()

    def _plot_pathway_categories(self, results, output_dir):
        """Categorize pathways into functional groups"""
        categories = {
            'Signal Transduction': ['signal', 'signaling', 'transduction', 'pathway'],
            'Cell Cycle': ['cell cycle', 'mitosis', 'G1', 'G2', 'M phase', 'S phase'],
            'DNA Repair': ['DNA repair', 'damage', 'replication', 'recombination'],
            'Apoptosis': ['apoptosis', 'cell death', 'programmed cell death'],
            'Metabolism': ['metabolic', 'metabolism', 'glycolysis', 'TCA'],
            'Immune': ['immune', 'inflammation', 'cytokine', 'interferon'],
            'Development': ['development', 'differentiation', 'morphogenesis'],
            'Other': []
        }

        category_counts = Counter()

        for db_name, df in results.items():
            for _, row in df.iterrows():
                term_lower = row['Term'].lower()
                categorized = False

                for cat, keywords in categories.items():
                    if cat != 'Other' and any(kw in term_lower for kw in keywords):
                        category_counts[cat] += 1
                        categorized = True
                        break

                if not categorized:
                    category_counts['Other'] += 1

        # Plot
        fig, ax = plt.subplots(figsize=(12, 8))

        cats = list(category_counts.keys())
        vals = [category_counts[c] for c in cats]

        # Sort by size for better visualization
        sorted_data = sorted(zip(cats, vals), key=lambda x: x[1], reverse=True)
        cats, vals = zip(*sorted_data)

        colors = plt.cm.Set3(range(len(cats)))

        # Only show percentage for slices > 2%, use legend for labels
        def autopct_format(pct):
            return f'{pct:.1f}%' if pct > 2 else ''

        wedges, texts, autotexts = ax.pie(vals, autopct=autopct_format,
                                            colors=colors, startangle=90,
                                            textprops={'fontsize': 11, 'fontweight': 'bold'},
                                            pctdistance=0.85)

        # Create legend with counts
        legend_labels = [f'{cat} (n={val})' for cat, val in zip(cats, vals)]
        ax.legend(wedges, legend_labels,
                 title="Pathway Categories",
                 loc="center left",
                 bbox_to_anchor=(1, 0, 0.5, 1),
                 fontsize=10,
                 title_fontsize=11)

        ax.set_title('Pathway Functional Categories',
                    fontsize=14, fontweight='bold', pad=30)

        plt.tight_layout()
        output_path = output_dir / 'pathway_categories.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close()

    def _plot_database_comparison(self, results, output_dir):
        """Compare pathway counts across databases"""
        db_counts = {db.replace('_', ' '): len(df) for db, df in results.items()}

        fig, ax = plt.subplots(figsize=(10, 6))

        dbs = list(db_counts.keys())
        counts = list(db_counts.values())

        bars = ax.bar(range(len(dbs)), counts, color='#16a085',
                     edgecolor='black', alpha=0.8)
        ax.set_xticks(range(len(dbs)))
        ax.set_xticklabels(dbs, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel('Number of Significant Pathways', fontsize=11, fontweight='bold')
        ax.set_title('Pathway Enrichment by Database', fontsize=13, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)

        for i, (bar, count) in enumerate(zip(bars, counts)):
            ax.text(i, count, f'{count}', ha='center', va='bottom',
                   fontsize=10, fontweight='bold')

        plt.tight_layout()
        output_path = output_dir / 'database_comparison.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close()

    def _plot_gene_pathway_network(self, results, output_dir, top_n=30):
        """Plot gene-pathway associations"""
        # Collect top pathways and their genes
        pathway_data = []

        for db_name, df in results.items():
            for _, row in df.head(5).iterrows():  # Top 5 from each DB
                pathway_data.append({
                    'Pathway': row['Term'][:40],
                    'Score': row['Combined Score'],
                    'Genes': row['Overlapping Genes'] if isinstance(row['Overlapping Genes'], list)
                            else str(row['Overlapping Genes']).split(';')[:5]
                })

        # Sort by score and take top N
        pathway_data.sort(key=lambda x: x['Score'], reverse=True)
        pathway_data = pathway_data[:top_n]

        # Create table visualization
        fig, ax = plt.subplots(figsize=(14, 18))
        ax.axis('tight')
        ax.axis('off')

        table_data = []
        for p in pathway_data:
            genes_str = ', '.join(p['Genes'][:5]) if isinstance(p['Genes'], list) else str(p['Genes'])[:50]
            table_data.append([p['Pathway'], f"{p['Score']:.1f}", genes_str])

        table = ax.table(cellText=table_data,
                        colLabels=['Pathway', 'Score', 'Top Genes'],
                        cellLoc='left',
                        loc='center',
                        colWidths=[0.4, 0.1, 0.5])

        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)

        # Style header
        for i in range(3):
            table[(0, i)].set_facecolor('#3498db')
            table[(0, i)].set_text_props(weight='bold', color='white')

        # Alternate row colors
        for i in range(1, len(table_data) + 1):
            for j in range(3):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#ecf0f1')

        ax.set_title('Top Pathways and Associated Genes',
                    fontsize=14, fontweight='bold', pad=30)

        plt.tight_layout()
        output_path = output_dir / 'gene_pathway_table.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close()

    def save_results(self, results, gbm_pathways, output_dir):
        """Save all results to files"""
        output_dir = Path(output_dir)

        # Save all pathway results
        for db_name, df in results.items():
            safe_name = db_name.replace(' ', '_').replace('/', '_')
            output_file = output_dir / f'pathway_{safe_name}.tsv'
            df.to_csv(output_file, sep='\t', index=False)

        # Save GBM-specific pathways
        if len(gbm_pathways) > 0:
            gbm_file = output_dir / 'gbm_relevant_pathways.tsv'
            gbm_pathways.to_csv(gbm_file, sep='\t', index=False)
            print(f"\n  Saved GBM pathways: {gbm_file}")


def load_occ_panel_genes(bed_file):
    """Load gene list from OCC panel BED file"""
    print(f"\nLoading OCC panel genes from: {bed_file}")

    target_genes = []
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 10:
                attributes = fields[9]
                for attr in attributes.split(';'):
                    if 'gene_name=' in attr:
                        gene_name = attr.split('gene_name=')[1]
                        target_genes.append(gene_name)
                        break

    target_genes = sorted(set(target_genes))
    print(f"  Loaded {len(target_genes)} genes from OCC panel")
    return target_genes


def load_genes_from_variants(variant_file, occ_genes):
    """Load genes from variant file, filtered to OCC panel"""
    print(f"\nLoading mutated genes from: {variant_file}")

    try:
        df = pd.read_csv(variant_file)
    except:
        df = pd.read_csv(variant_file, sep='\t')

    # Find gene column
    gene_col = None
    for col in df.columns:
        if 'gene' in col.lower() and 'refgene' in col.lower():
            gene_col = col
            break

    if gene_col is None:
        raise ValueError("Could not find gene column")

    # Get unique genes from variants
    all_mutated_genes = df[gene_col].dropna().unique().tolist()
    print(f"  Total mutated genes in data: {len(all_mutated_genes)}")

    # Filter to only OCC panel genes
    genes = [g for g in all_mutated_genes if g in occ_genes]
    print(f"  OCC panel genes with mutations: {len(genes)}")

    # Show gene counts
    df_occ = df[df[gene_col].isin(genes)]
    gene_counts = df_occ[gene_col].value_counts().head(15)
    print(f"\n  Top 15 mutated OCC panel genes:")
    for gene, count in gene_counts.items():
        print(f"    {gene}: {count} mutations")

    return genes


def main():
    """Main GBM pathway analysis workflow"""
    print("="*80)
    print("GBM PATHWAY ENRICHMENT ANALYSIS - OCC PANEL GENES")
    print("="*80)

    # Paths
    variant_file = "/home/chbope/extension/script/mut_landscape/results/combined_variants_all_samples.csv"
    occ_bed_file = "/home/chbope/extension/nWGS_manuscript_data/data/reference/OCC.protein_coding.bed"
    output_dir = Path("/home/chbope/extension/script/mut_landscape/results/gbm_pathway_analysisclairs")
    output_dir.mkdir(exist_ok=True)

    print(f"\nVariant file: {variant_file}")
    print(f"OCC panel: {occ_bed_file}")
    print(f"Output: {output_dir}")

    # Load OCC panel genes
    occ_genes = load_occ_panel_genes(occ_bed_file)

    # Load mutated genes (filtered to OCC panel)
    genes = load_genes_from_variants(variant_file, occ_genes)

    # Initialize analyzer
    analyzer = GBMPathwayAnalyzer(genes)

    # Run enrichment
    print("\n" + "="*80)
    print("RUNNING ENRICHMENT ANALYSIS")
    print("="*80)
    enrichr_results = analyzer.run_enrichr_analysis()

    if enrichr_results:
        # Parse results
        print("\n" + "="*80)
        print("PARSING RESULTS")
        print("="*80)
        parsed_results = analyzer.parse_results(max_pathways=30, p_threshold=0.05)

        if parsed_results:
            # Identify GBM pathways
            gbm_pathways = analyzer.identify_gbm_pathways(enrichr_results)

            # Create visualizations
            analyzer.create_gbm_visualizations(parsed_results, gbm_pathways, output_dir)

            # Save results
            print("\n" + "="*80)
            print("SAVING RESULTS")
            print("="*80)
            analyzer.save_results(parsed_results, gbm_pathways, output_dir)

            print("\n" + "="*80)
            print("GBM PATHWAY ANALYSIS COMPLETE!")
            print("="*80)
            print(f"\nResults: {output_dir}")
            print(f"Databases: {len(parsed_results)}")
            print(f"Total pathways: {sum(len(df) for df in parsed_results.values())}")
            print(f"GBM-relevant: {len(gbm_pathways)}")
            print("="*80)
        else:
            print("\nNo significant pathways found")
    else:
        print("\nEnrichment analysis failed")


if __name__ == "__main__":
    main()
