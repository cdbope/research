#!/usr/bin/env python3
"""
GBM-Focused Pathway Enrichment Analysis - Version 2 (Per-Sample)
Analyzes pathway enrichment for each individual GBM sample separately
Uses PCGR-annotated variant files for each sample
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict
import requests
import json
from glob import glob
import time

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

class GBMPathwayAnalyzerPerSample:
    """
    Per-sample GBM pathway enrichment analysis
    """

    def __init__(self):
        """Initialize per-sample pathway analyzer"""

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

        self.all_sample_results = []

        # Pathway databases
        self.databases = [
            'KEGG_2021_Human',
            'GO_Biological_Process_2021',
            'Reactome_2022',
            'WikiPathway_2021_Human',
            'MSigDB_Hallmark_2020',
            'BioPlanet_2019',
            'OMIM_Disease',
            'DisGeNET',
            'Jensen_DISEASES',
        ]

    def load_sample_genes(self, sample_file, occ_genes):
        """Load genes from a single sample PCGR file"""
        try:
            df = pd.read_csv(sample_file, sep='\t')
        except:
            try:
                df = pd.read_csv(sample_file)
            except Exception as e:
                print(f"    Error reading file: {e}")
                return []

        # Get unique genes
        if 'Gene.refGene' in df.columns:
            genes = df['Gene.refGene'].dropna().unique().tolist()
            # Filter to OCC panel genes
            genes = [g for g in genes if g in occ_genes]
            return genes
        else:
            print(f"    Warning: No Gene.refGene column found")
            return []

    def run_enrichr_for_sample(self, sample_id, genes):
        """Run Enrichr analysis for a single sample"""
        if len(genes) < 3:
            print(f"    {sample_id}: Too few genes ({len(genes)}), skipping")
            return None

        # Enrichr API endpoints
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
        ENRICHR_QUERY_URL = 'https://maayanlab.cloud/Enrichr/enrich'

        # Submit gene list
        genes_str = '\n'.join(genes)
        payload = {
            'list': (None, genes_str),
            'description': (None, f'GBM_Sample_{sample_id}')
        }

        try:
            response = requests.post(ENRICHR_URL, files=payload, timeout=30)
            if not response.ok:
                print(f"    {sample_id}: Error submitting genes ({response.status_code})")
                return None

            data = json.loads(response.text)
            user_list_id = data['userListId']

            # Query databases
            sample_results = {
                'sample_id': sample_id,
                'n_genes': len(genes),
                'genes': genes,
                'pathways': {}
            }

            for db in self.databases:
                query_url = f"{ENRICHR_QUERY_URL}?userListId={user_list_id}&backgroundType={db}"
                response = requests.get(query_url, timeout=30)

                if response.ok:
                    results = json.loads(response.text)
                    if db in results and len(results[db]) > 0:
                        # Parse and filter pathways
                        pathway_list = []
                        for pathway_data in results[db]:
                            if pathway_data[6] < 0.05:  # Adjusted p-value < 0.05
                                pathway_list.append({
                                    'Term': pathway_data[1],
                                    'P-value': pathway_data[2],
                                    'Adjusted P-value': pathway_data[6],
                                    'Combined Score': pathway_data[4],
                                    'Genes': pathway_data[5]
                                })

                        if pathway_list:
                            sample_results['pathways'][db] = pathway_list

            # Add delay to avoid rate limiting
            time.sleep(0.5)

            return sample_results

        except Exception as e:
            print(f"    {sample_id}: Error - {e}")
            return None

    def identify_gbm_pathways_per_sample(self, sample_results):
        """Identify GBM-relevant pathways for a sample"""
        gbm_pathways = []

        for db_name, pathways in sample_results['pathways'].items():
            for pathway in pathways:
                term = pathway['Term'].lower()

                # Check if pathway is GBM-relevant
                is_gbm_relevant = any(keyword.lower() in term for keyword in self.gbm_key_pathways)

                if is_gbm_relevant:
                    gbm_pathways.append({
                        'Sample_ID': sample_results['sample_id'],
                        'Database': db_name,
                        'Term': pathway['Term'],
                        'P-value': pathway['P-value'],
                        'Adjusted P-value': pathway['Adjusted P-value'],
                        'Combined Score': pathway['Combined Score'],
                        'Genes': pathway['Genes']
                    })

        return gbm_pathways

    def analyze_all_samples(self, pcgr_dir, occ_genes, output_dir, max_samples=None):
        """Analyze pathway enrichment for all samples"""
        # Find all PCGR files
        pcgr_files = sorted(glob(f"{pcgr_dir}/*_annotateandfilter_clairsto.csv"))

        if max_samples:
            pcgr_files = pcgr_files[:max_samples]

        print(f"\nFound {len(pcgr_files)} PCGR sample files")
        print(f"Processing samples for pathway enrichment...\n")

        all_gbm_pathways = []
        sample_summary = []

        for i, sample_file in enumerate(pcgr_files, 1):
            sample_id = Path(sample_file).name.replace('_annotateandfilter_clairsto.csv', '')
            print(f"[{i}/{len(pcgr_files)}] {sample_id}...")

            # Load genes for this sample
            genes = self.load_sample_genes(sample_file, occ_genes)

            if len(genes) == 0:
                print(f"    No OCC panel genes found, skipping")
                continue

            print(f"    Found {len(genes)} OCC panel genes")

            # Run enrichment for this sample
            sample_results = self.run_enrichr_for_sample(sample_id, genes)

            if sample_results is None:
                continue

            # Store results
            self.all_sample_results.append(sample_results)

            # Count pathways
            n_pathways = sum(len(p) for p in sample_results['pathways'].values())

            # Identify GBM pathways
            gbm_pathways = self.identify_gbm_pathways_per_sample(sample_results)
            all_gbm_pathways.extend(gbm_pathways)

            print(f"    Enrichment: {n_pathways} total pathways, {len(gbm_pathways)} GBM-relevant")

            # Summary
            sample_summary.append({
                'Sample_ID': sample_id,
                'N_Genes': len(genes),
                'N_Total_Pathways': n_pathways,
                'N_GBM_Pathways': len(gbm_pathways)
            })

        print(f"\n{'='*80}")
        print(f"COMPLETED: {len(self.all_sample_results)}/{len(pcgr_files)} samples analyzed")
        print(f"Total GBM-relevant pathway enrichments: {len(all_gbm_pathways)}")
        print(f"{'='*80}\n")

        return pd.DataFrame(all_gbm_pathways), pd.DataFrame(sample_summary)

    def create_summary_visualizations(self, gbm_df, summary_df, output_dir):
        """Create summary visualizations across all samples"""
        print(f"\n{'='*80}")
        print("CREATING SUMMARY VISUALIZATIONS")
        print(f"{'='*80}\n")

        # 1. Pathway recurrence across samples
        print("[1/5] Pathway recurrence across samples...")
        pathway_counts = gbm_df['Term'].value_counts().head(30)

        fig, ax = plt.subplots(figsize=(14, 10))
        bars = ax.barh(range(len(pathway_counts)), pathway_counts.values,
                      color='#3498db', edgecolor='black', alpha=0.8)

        pathway_labels = [term[:60] + '...' if len(term) > 60 else term
                         for term in pathway_counts.index]
        ax.set_yticks(range(len(pathway_counts)))
        ax.set_yticklabels(pathway_labels, fontsize=9)
        ax.set_xlabel('Number of Samples with Pathway Enrichment', fontsize=12, fontweight='bold')
        ax.set_title('Most Recurrent GBM Pathways Across Samples',
                    fontsize=14, fontweight='bold', pad=20, color='darkred')
        ax.grid(axis='x', alpha=0.3)

        for i, count in enumerate(pathway_counts.values):
            ax.text(count, i, f' {count}', va='center', fontsize=9, fontweight='bold')

        plt.tight_layout()
        output_path = output_dir / 'pathway_recurrence_across_samples.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"   Saved: {output_path}")
        plt.close()

        # 2. Per-sample pathway enrichment distribution
        print("[2/5] Per-sample pathway enrichment distribution...")
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

        # Total pathways per sample
        sorted_summary = summary_df.sort_values('N_Total_Pathways', ascending=False)
        ax1.bar(range(len(sorted_summary)), sorted_summary['N_Total_Pathways'].values,
               color='#2ecc71', edgecolor='black', alpha=0.7)
        ax1.set_xlabel('Sample Rank', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Number of Enriched Pathways', fontsize=12, fontweight='bold')
        ax1.set_title('Total Pathway Enrichment per Sample',
                     fontsize=13, fontweight='bold', pad=15)
        ax1.axhline(sorted_summary['N_Total_Pathways'].mean(), color='red',
                   linestyle='--', label=f'Mean: {sorted_summary["N_Total_Pathways"].mean():.1f}')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)

        # GBM pathways per sample
        sorted_summary_gbm = summary_df.sort_values('N_GBM_Pathways', ascending=False)
        ax2.bar(range(len(sorted_summary_gbm)), sorted_summary_gbm['N_GBM_Pathways'].values,
               color='#e74c3c', edgecolor='black', alpha=0.7)
        ax2.set_xlabel('Sample Rank', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Number of GBM-Relevant Pathways', fontsize=12, fontweight='bold')
        ax2.set_title('GBM-Relevant Pathway Enrichment per Sample',
                     fontsize=13, fontweight='bold', pad=15)
        ax2.axhline(sorted_summary_gbm['N_GBM_Pathways'].mean(), color='blue',
                   linestyle='--', label=f'Mean: {sorted_summary_gbm["N_GBM_Pathways"].mean():.1f}')
        ax2.legend()
        ax2.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        output_path = output_dir / 'per_sample_pathway_distribution.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"   Saved: {output_path}")
        plt.close()

        # 3. Database contribution across samples
        print("[3/5] Database contribution analysis...")
        db_counts = gbm_df['Database'].value_counts()

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

        # Pie chart
        colors = plt.cm.Set3(range(len(db_counts)))
        wedges, texts, autotexts = ax1.pie(db_counts.values,
                                            autopct='%1.1f%%',
                                            startangle=90, colors=colors,
                                            textprops={'fontsize': 10, 'fontweight': 'bold'})
        ax1.set_title('GBM Pathway Enrichments by Database',
                     fontsize=13, fontweight='bold', pad=15)

        legend_labels = [f'{db.replace("_", " ")} (n={count})'
                        for db, count in db_counts.items()]
        ax1.legend(wedges, legend_labels,
                  loc="center left",
                  bbox_to_anchor=(1, 0, 0.5, 1),
                  fontsize=9)

        # Bar chart
        bars = ax2.bar(range(len(db_counts)), db_counts.values,
                      color=colors, edgecolor='black', alpha=0.8)
        ax2.set_xticks(range(len(db_counts)))
        ax2.set_xticklabels([db.replace('_', ' ') for db in db_counts.index],
                           rotation=45, ha='right', fontsize=9)
        ax2.set_ylabel('Number of Pathway Enrichments', fontsize=11, fontweight='bold')
        ax2.set_title('Database Contribution to GBM Pathways',
                     fontsize=13, fontweight='bold', pad=15)
        ax2.grid(axis='y', alpha=0.3)

        for i, count in enumerate(db_counts.values):
            ax2.text(i, count, f'{count}', ha='center', va='bottom',
                    fontsize=9, fontweight='bold')

        plt.tight_layout()
        output_path = output_dir / 'database_contribution_across_samples.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"   Saved: {output_path}")
        plt.close()

        # 4. Gene involvement in pathways
        print("[4/5] Gene involvement analysis...")
        all_genes = []
        for genes_str in gbm_df['Genes'].dropna():
            if isinstance(genes_str, str):
                genes = genes_str.split(';')
                all_genes.extend([g.strip() for g in genes if g.strip()])

        gene_counts = Counter(all_genes).most_common(30)

        if gene_counts:
            fig, ax = plt.subplots(figsize=(12, 10))
            genes, counts = zip(*gene_counts)

            bars = ax.barh(range(len(genes)), counts,
                          color='#9b59b6', edgecolor='black', alpha=0.8)
            ax.set_yticks(range(len(genes)))
            ax.set_yticklabels(genes, fontsize=10, fontweight='bold')
            ax.set_xlabel('Number of Pathway Enrichments', fontsize=12, fontweight='bold')
            ax.set_title('Top Genes Involved in GBM Pathway Enrichments',
                        fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='x', alpha=0.3)

            for i, count in enumerate(counts):
                ax.text(count, i, f' {count}', va='center', fontsize=9, fontweight='bold')

            plt.tight_layout()
            output_path = output_dir / 'top_genes_in_pathways.png'
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"   Saved: {output_path}")
            plt.close()

        # 5. Summary statistics
        print("[5/5] Summary statistics...")
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('tight')
        ax.axis('off')

        stats_data = [
            ['Total Samples Analyzed', str(len(summary_df))],
            ['Total GBM Pathway Enrichments', str(len(gbm_df))],
            ['Unique GBM Pathways', str(gbm_df['Term'].nunique())],
            ['Mean Pathways per Sample', f"{summary_df['N_Total_Pathways'].mean():.1f}"],
            ['Mean GBM Pathways per Sample', f"{summary_df['N_GBM_Pathways'].mean():.1f}"],
            ['Max Pathways in Single Sample', str(summary_df['N_Total_Pathways'].max())],
            ['Min Pathways in Single Sample', str(summary_df['N_Total_Pathways'].min())],
            ['Most Recurrent Pathway', pathway_counts.index[0] if len(pathway_counts) > 0 else 'N/A'],
            ['Pathway Recurrence (top)', str(pathway_counts.values[0]) + ' samples' if len(pathway_counts) > 0 else 'N/A'],
        ]

        table = ax.table(cellText=stats_data,
                        colLabels=['Metric', 'Value'],
                        cellLoc='left',
                        loc='center',
                        colWidths=[0.6, 0.4])

        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 2.5)

        # Style header
        for i in range(2):
            table[(0, i)].set_facecolor('#3498db')
            table[(0, i)].set_text_props(weight='bold', color='white', fontsize=12)

        # Alternate row colors
        for i in range(1, len(stats_data) + 1):
            for j in range(2):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#ecf0f1')

        ax.set_title('Per-Sample Pathway Enrichment Summary',
                    fontsize=16, fontweight='bold', pad=30)

        plt.tight_layout()
        output_path = output_dir / 'summary_statistics_table.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"   Saved: {output_path}")
        plt.close()

        print(f"\n{'='*80}")
        print("All visualizations generated successfully!")
        print(f"{'='*80}")


def load_occ_panel_genes(bed_file):
    """Load gene list from OCC panel BED file"""
    print(f"Loading OCC panel genes from: {bed_file}")

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
    print(f"  Loaded {len(target_genes)} genes from OCC panel\n")
    return target_genes


def main():
    """Main per-sample GBM pathway analysis workflow"""
    print("="*80)
    print("GBM PATHWAY ENRICHMENT ANALYSIS - PER-SAMPLE (VERSION 2)")
    print("="*80)

    # Paths
    pcgr_dir = "/home/chbope/extension/data/200GMBs/pcgr"
    occ_bed_file = "/home/chbope/extension/nWGS_manuscript_data/data/reference/OCC.protein_coding.bed"
    output_dir = Path("/home/chbope/extension/script/mut_landscape/results/gbm_pathway_analysis_v2_persample")
    output_dir.mkdir(exist_ok=True)

    print(f"\nPCGR directory: {pcgr_dir}")
    print(f"OCC panel: {occ_bed_file}")
    print(f"Output: {output_dir}")

    # Load OCC panel genes
    occ_genes = load_occ_panel_genes(occ_bed_file)

    # Initialize analyzer
    analyzer = GBMPathwayAnalyzerPerSample()

    # Analyze all samples (set max_samples to limit for testing, None for all)
    gbm_df, summary_df = analyzer.analyze_all_samples(
        pcgr_dir,
        occ_genes,
        output_dir,
        max_samples=None  # Set to e.g., 20 for testing, None for all samples
    )

    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)

    if len(gbm_df) > 0:
        # Save GBM pathways
        gbm_file = output_dir / 'gbm_pathways_per_sample.tsv'
        gbm_df.to_csv(gbm_file, sep='\t', index=False)
        print(f"  Saved GBM pathways: {gbm_file}")

        # Save summary
        summary_file = output_dir / 'sample_summary.tsv'
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"  Saved sample summary: {summary_file}")

        # Create visualizations
        analyzer.create_summary_visualizations(gbm_df, summary_df, output_dir)

        print("\n" + "="*80)
        print("PER-SAMPLE PATHWAY ANALYSIS COMPLETE!")
        print("="*80)
        print(f"\nResults: {output_dir}")
        print(f"Samples analyzed: {len(summary_df)}")
        print(f"Total GBM pathway enrichments: {len(gbm_df)}")
        print(f"Unique GBM pathways: {gbm_df['Term'].nunique()}")
        print("="*80)

    else:
        print("\nNo pathways found!")


if __name__ == "__main__":
    main()
