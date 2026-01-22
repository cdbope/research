#!/usr/bin/env python3
"""
Create Visual Benchmark Report with Gene Illustrations
"""

import pandas as pd
import sys
from pathlib import Path
from collections import Counter

def load_sample_data(sample_id, ds_dir, clair_dir):
    """Load detailed variant data for a sample"""
    ds_file = ds_dir / f"{sample_id}_annotateandfilter_deep_somatic.csv"
    clair_file = clair_dir / f"{sample_id}_annotateandfilter_clairsto.csv"

    ds_df = pd.read_csv(ds_file, sep='\t') if ds_file.exists() else None
    clair_df = pd.read_csv(clair_file, sep='\t') if clair_file.exists() else None

    return ds_df, clair_df

def extract_genes(df):
    """Extract gene list from dataframe"""
    if df is None or 'Gene.refGene' not in df.columns:
        return []

    genes = []
    for gene_str in df['Gene.refGene'].dropna():
        genes.extend([g.strip() for g in str(gene_str).split(',')])
    return genes

def create_visual_report():
    benchmark_dir = Path("/home/chbope/extension/script/deepsomatic/benchmark")
    ds_dir = benchmark_dir / "deepsomatic"
    clair_dir = benchmark_dir / "clairsto"
    output_dir = benchmark_dir / "comparison_results"

    # Load summary data
    summary_df = pd.read_csv(output_dir / "per_sample_comparison.csv")

    print("="*100)
    print("VISUAL BENCHMARK REPORT: ClairS-TO vs DeepSomatic")
    print("="*100)
    print()

    # Collect all genes
    all_ds_genes = []
    all_clair_genes = []

    for _, row in summary_df.iterrows():
        sample_id = row['sample_id']
        ds_df, clair_df = load_sample_data(sample_id, ds_dir, clair_dir)

        all_ds_genes.extend(extract_genes(ds_df))
        all_clair_genes.extend(extract_genes(clair_df))

    # Count gene frequencies
    ds_gene_counts = Counter(all_ds_genes)
    clair_gene_counts = Counter(all_clair_genes)

    # Get all unique genes
    all_genes = set(ds_gene_counts.keys()) | set(clair_gene_counts.keys())

    # Create gene comparison
    output = []

    output.append("\n" + "="*100)
    output.append("1. VARIANT DETECTION VISUAL COMPARISON")
    output.append("="*100)
    output.append("")

    # Bar chart for total variants
    ds_total = summary_df['ds_variants'].sum()
    clair_total = summary_df['clair_variants'].sum()
    max_val = max(ds_total, clair_total)

    output.append("Total Variants Detected Across All Samples:")
    output.append("")

    ds_bar = "█" * int((ds_total / max_val) * 50)
    clair_bar = "█" * int((clair_total / max_val) * 50)

    output.append(f"  DeepSomatic:  {ds_bar} {ds_total}")
    output.append(f"  ClairS-TO:    {clair_bar} {clair_total}")
    output.append("")
    output.append(f"  Difference:   DeepSomatic detected {ds_total - clair_total} more variants (+{(ds_total/clair_total - 1)*100:.1f}%)")

    # Per-sample visualization
    output.append("\n" + "="*100)
    output.append("2. PER-SAMPLE VARIANT DETECTION")
    output.append("="*100)
    output.append("")

    output.append(f"{'Sample':<12} {'DeepSomatic':<25} {'ClairS-TO':<25} {'Concordance':<15}")
    output.append("-"*100)

    for _, row in summary_df.iterrows():
        sample = row['sample_id']
        ds_cnt = row['ds_variants']
        clair_cnt = row['clair_variants']
        conc = row['concordance']

        ds_bar = "█" * ds_cnt
        clair_bar = "█" * clair_cnt

        output.append(f"{sample:<12} {ds_bar:<25} {clair_bar:<25} {conc:>6.1f}%")

    # Gene frequency analysis
    output.append("\n" + "="*100)
    output.append("3. CANCER GENES DETECTED (Top 20 Most Frequent)")
    output.append("="*100)
    output.append("")

    # Combine and sort genes
    all_gene_data = {}
    for gene in all_genes:
        ds_count = ds_gene_counts.get(gene, 0)
        clair_count = clair_gene_counts.get(gene, 0)
        total = ds_count + clair_count
        all_gene_data[gene] = (ds_count, clair_count, total)

    # Sort by total frequency
    sorted_genes = sorted(all_gene_data.items(), key=lambda x: x[1][2], reverse=True)[:20]

    output.append(f"{'Gene':<15} {'DeepSomatic':<30} {'ClairS-TO':<30} {'Total':<10}")
    output.append("-"*100)

    max_count = max([data[2] for _, data in sorted_genes])

    for gene, (ds_cnt, clair_cnt, total) in sorted_genes:
        ds_bar = "█" * int((ds_cnt / max_count) * 20) if ds_cnt > 0 else ""
        clair_bar = "█" * int((clair_cnt / max_count) * 20) if clair_cnt > 0 else ""

        output.append(f"{gene:<15} {ds_bar:<20} ({ds_cnt:>2})  {clair_bar:<20} ({clair_cnt:>2})  {total:>3}")

    # Shared vs unique genes
    output.append("\n" + "="*100)
    output.append("4. GENE DETECTION VENN DIAGRAM")
    output.append("="*100)
    output.append("")

    ds_genes = set(ds_gene_counts.keys())
    clair_genes = set(clair_gene_counts.keys())
    shared_genes = ds_genes & clair_genes
    ds_only_genes = ds_genes - clair_genes
    clair_only_genes = clair_genes - ds_genes

    output.append("                    ┌─────────────────────────────────┐")
    output.append("                    │      DeepSomatic Only           │")
    output.append(f"                    │         {len(ds_only_genes):>3} genes               │")
    output.append("    ┌───────────────┼─────────────────────────────────┼───────────────┐")
    output.append("    │               │                                 │               │")
    output.append("    │  DeepSomatic  │         Shared Genes            │  ClairS-TO    │")
    output.append(f"    │   {len(ds_genes):>3} genes   │          {len(shared_genes):>3} genes            │  {len(clair_genes):>3} genes   │")
    output.append("    │               │                                 │               │")
    output.append("    └───────────────┼─────────────────────────────────┼───────────────┘")
    output.append(f"                    │       ClairS-TO Only            │")
    output.append(f"                    │          {len(clair_only_genes):>3} genes              │")
    output.append("                    └─────────────────────────────────┘")
    output.append("")

    # Show DeepSomatic-only genes
    if ds_only_genes:
        output.append("DeepSomatic-Only Genes (not detected by ClairS-TO):")
        gene_list = ', '.join(sorted(ds_only_genes))
        output.append(f"  {gene_list}")
        output.append("")

    # Show ClairS-TO-only genes
    if clair_only_genes:
        output.append("ClairS-TO-Only Genes (not detected by DeepSomatic):")
        gene_list = ', '.join(sorted(clair_only_genes))
        output.append(f"  {gene_list}")
        output.append("")

    # Cancer gene categories
    output.append("\n" + "="*100)
    output.append("5. CANCER GENE CATEGORIES")
    output.append("="*100)
    output.append("")

    # Define cancer gene categories
    categories = {
        'Oncogenes': ['KRAS', 'NRAS', 'HRAS', 'BRAF', 'MYC', 'EGFR', 'ERBB2', 'PIK3CA', 'AKT1', 'MET', 'ALK', 'RET', 'ROS1'],
        'Tumor Suppressors': ['TP53', 'PTEN', 'RB1', 'APC', 'VHL', 'NF1', 'TSC1', 'TSC2', 'CDKN2A', 'BRCA1', 'BRCA2', 'STK11'],
        'DNA Repair': ['MLH1', 'MSH2', 'MSH6', 'PMS2', 'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2'],
        'Epigenetic': ['DNMT3A', 'TET2', 'IDH1', 'IDH2', 'ARID1A', 'SMARCA4', 'KMT2A', 'EZH2'],
        'Cell Cycle': ['CDKN2A', 'CDKN2B', 'CDK4', 'CDK6', 'CCND1', 'RB1', 'TP53'],
        'Growth Factor': ['EGFR', 'ERBB2', 'MET', 'FGFR1', 'FGFR2', 'FGFR3', 'PDGFRA', 'KIT']
    }

    for category, genes in categories.items():
        ds_in_cat = [g for g in genes if g in ds_genes]
        clair_in_cat = [g for g in genes if g in clair_genes]

        output.append(f"{category}:")
        output.append(f"  DeepSomatic: {len(ds_in_cat)}/{len(genes)} - {', '.join(ds_in_cat) if ds_in_cat else 'None'}")
        output.append(f"  ClairS-TO:   {len(clair_in_cat)}/{len(genes)} - {', '.join(clair_in_cat) if clair_in_cat else 'None'}")
        output.append("")

    # VAF distribution
    output.append("\n" + "="*100)
    output.append("6. VAF (VARIANT ALLELE FREQUENCY) DISTRIBUTION")
    output.append("="*100)
    output.append("")

    output.append("Minimum VAF Detected (Lower is Better):")
    output.append("")

    ds_min = summary_df['ds_min_vaf'].min()
    clair_min = summary_df['clair_min_vaf'].min()

    output.append("  VAF Scale:  0%    5%    10%   15%   20%   25%   30%   35%   40%")
    output.append("              |-----|-----|-----|-----|-----|-----|-----|-----|")

    ds_pos = int((ds_min / 0.40) * 60)
    clair_pos = int((clair_min / 0.40) * 60)

    output.append(f"  DeepSomatic {'─'*ds_pos}▼ {ds_min*100:.2f}%")
    output.append(f"  ClairS-TO   {'─'*clair_pos}▼ {clair_min*100:.2f}%")
    output.append("")
    output.append(f"  DeepSomatic can detect {clair_min - ds_min:.4f} ({(clair_min - ds_min)*100:.2f}%) lower VAF ✓")

    # Low VAF variants
    output.append("\n" + "="*100)
    output.append("7. LOW VAF VARIANTS DISTRIBUTION (<30%)")
    output.append("="*100)
    output.append("")

    ds_low = summary_df['ds_low_vaf'].sum()
    clair_low = summary_df['clair_low_vaf'].sum()
    max_low = max(ds_low, clair_low)

    ds_low_bar = "█" * int((ds_low / max_low) * 50) if ds_low > 0 else ""
    clair_low_bar = "█" * int((clair_low / max_low) * 50) if clair_low > 0 else ""

    output.append(f"  DeepSomatic:  {ds_low_bar} {ds_low} variants")
    output.append(f"  ClairS-TO:    {clair_low_bar} {clair_low} variants")
    output.append("")
    output.append(f"  DeepSomatic detected {ds_low - clair_low} more low VAF variants (+{(ds_low/clair_low - 1)*100:.1f}%)")

    # Per-sample gene comparison
    output.append("\n" + "="*100)
    output.append("8. PER-SAMPLE GENE COMPARISON")
    output.append("="*100)
    output.append("")

    for _, row in summary_df.iterrows():
        sample_id = row['sample_id']
        ds_df, clair_df = load_sample_data(sample_id, ds_dir, clair_dir)

        ds_sample_genes = set(extract_genes(ds_df))
        clair_sample_genes = set(extract_genes(clair_df))
        shared_sample_genes = ds_sample_genes & clair_sample_genes
        ds_only_sample = ds_sample_genes - clair_sample_genes
        clair_only_sample = clair_sample_genes - ds_sample_genes

        output.append(f"\n{sample_id}:")
        output.append(f"  {'─'*90}")

        # Visual gene count comparison
        max_genes = max(len(ds_sample_genes), len(clair_sample_genes), 1)
        ds_gene_bar = "█" * int((len(ds_sample_genes) / max_genes) * 30)
        clair_gene_bar = "█" * int((len(clair_sample_genes) / max_genes) * 30)

        output.append(f"  DeepSomatic: {ds_gene_bar:<30} {row['ds_variants']} variants, {len(ds_sample_genes)} genes")
        output.append(f"  ClairS-TO:   {clair_gene_bar:<30} {row['clair_variants']} variants, {len(clair_sample_genes)} genes")

        # Concordance bar
        conc_bar_filled = "█" * int(row['concordance'] / 5)
        conc_bar_empty = "░" * (20 - int(row['concordance'] / 5))
        output.append(f"  Concordance: {conc_bar_filled}{conc_bar_empty} {row['concordance']:.1f}%")
        output.append("")

        # Mini Venn diagram for genes
        output.append(f"  Gene Overlap Diagram:")
        output.append(f"  ┌────────────────┬────────────────┬────────────────┐")
        output.append(f"  │  DS Only: {len(ds_only_sample):>2}   │  Shared: {len(shared_sample_genes):>2}    │  Clair Only: {len(clair_only_sample):>2} │")
        output.append(f"  └────────────────┴────────────────┴────────────────┘")
        output.append("")

        if shared_sample_genes:
            output.append(f"  ✓ Shared genes ({len(shared_sample_genes)}): {', '.join(sorted(shared_sample_genes))}")
        else:
            output.append(f"  ✗ Shared genes: None")

        if ds_only_sample:
            output.append(f"  ◆ DeepSomatic only ({len(ds_only_sample)}): {', '.join(sorted(ds_only_sample))}")

        if clair_only_sample:
            output.append(f"  ◇ ClairS-TO only ({len(clair_only_sample)}): {', '.join(sorted(clair_only_sample))}")

    # Concordance distribution
    output.append("\n" + "="*100)
    output.append("9. CONCORDANCE DISTRIBUTION")
    output.append("="*100)
    output.append("")

    conc_bins = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]
    output.append("Concordance Range Distribution:")
    output.append("")

    for low, high in conc_bins:
        count = ((summary_df['concordance'] >= low) & (summary_df['concordance'] < high)).sum()
        bar = "█" * count
        pct = (count / len(summary_df)) * 100
        output.append(f"  {low:>3}-{high:<3}%:  {bar:<15} {count} samples ({pct:.1f}%)")

    # Summary scorecard
    output.append("\n" + "="*100)
    output.append("10. FINAL SCORECARD")
    output.append("="*100)
    output.append("")

    metrics = {
        'Total Variants': (summary_df['ds_variants'].sum(), summary_df['clair_variants'].sum()),
        'Avg per Sample': (summary_df['ds_variants'].mean(), summary_df['clair_variants'].mean()),
        'Unique Genes': (len(ds_genes), len(clair_genes)),
        'Cancer Genes': (summary_df['ds_cancer_genes'].sum(), summary_df['clair_cancer_genes'].sum()),
        'Pathogenic': (summary_df['ds_pathogenic'].sum(), summary_df['clair_pathogenic'].sum()),
        'COSMIC Hits': (summary_df['ds_cosmic'].sum(), summary_df['clair_cosmic'].sum()),
        'Low VAF (<30%)': (ds_low, clair_low)
    }

    ds_wins = 0
    clair_wins = 0

    output.append(f"{'Metric':<20} {'DeepSomatic':>15} {'ClairS-TO':>15} {'Winner':>15}")
    output.append("-"*70)

    for metric, (ds_val, clair_val) in metrics.items():
        winner = '✓ DS' if ds_val > clair_val else '✓ Clair' if clair_val > ds_val else 'Tie'
        if ds_val > clair_val:
            ds_wins += 1
        elif clair_val > ds_val:
            clair_wins += 1

        if isinstance(ds_val, float):
            output.append(f"{metric:<20} {ds_val:>15.1f} {clair_val:>15.1f} {winner:>15}")
        else:
            output.append(f"{metric:<20} {ds_val:>15} {clair_val:>15} {winner:>15}")

    output.append("")
    output.append(f"FINAL SCORE: DeepSomatic {ds_wins} - {clair_wins} ClairS-TO")
    output.append("")

    winner = "DeepSomatic" if ds_wins > clair_wins else "ClairS-TO" if clair_wins > ds_wins else "TIE"
    output.append(f"🏆 WINNER: {winner}")

    output.append("\n" + "="*100)

    # Write to file
    report_text = '\n'.join(output)

    with open(output_dir / "VISUAL_BENCHMARK_REPORT.txt", 'w') as f:
        f.write(report_text)

    print(report_text)

    print(f"\n\nVisual report saved to: {output_dir / 'VISUAL_BENCHMARK_REPORT.txt'}")
    print("="*100)

    return 0

if __name__ == "__main__":
    sys.exit(create_visual_report())
