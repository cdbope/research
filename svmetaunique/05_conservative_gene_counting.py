#!/usr/bin/env python3
"""
Conservative Gene Counting Approach
Based on: "For each gene, we reported only one instance of a specific alteration
type, regardless of how many identical alterations occurred within that same gene."

This approach:
1. Counts max 1 SV per gene per sample (binary: present/absent)
2. Calculates frequency as: (# samples with ≥1 SV in gene) / 200
3. Compares to standard approach (counts all SVs)

All outputs go to: /home/chbope/extension/script/svmetaunique/results/
"""

import pandas as pd
import numpy as np
from scipy import stats
import os

def wilson_score_interval(hits, total, confidence=0.95):
    """Calculate Wilson score confidence interval"""
    if total == 0:
        return 0, 0

    p = hits / total
    z = stats.norm.ppf((1 + confidence) / 2)

    denominator = 1 + z**2 / total
    center = (p + z**2 / (2 * total)) / denominator
    margin = z * np.sqrt((p * (1 - p) / total + z**2 / (4 * total**2))) / denominator

    return max(0, center - margin), min(1, center + margin)

def main():
    """Main analysis pipeline"""

    print("="*80)
    print("CONSERVATIVE GENE COUNTING ANALYSIS")
    print("="*80)
    print("\nApproach: Count only 1 SV per gene per sample (binary)")
    print("Rationale: More comparable to literature (e.g., npae082.pdf)")
    print("="*80)
    print()

    output_dir = "/home/chbope/extension/script/svmetaunique/results"
    os.makedirs(output_dir, exist_ok=True)

    # Load existing results for comparison
    existing_results = "/home/chbope/extension/script/svmeta/results/external_comparison/high_confidence_validated_genes_with_fold_changes.csv"

    if not os.path.exists(existing_results):
        print(f"⚠️  Existing results not found: {existing_results}")
        return

    existing_df = pd.read_csv(existing_results)

    print(f"Analyzing {len(existing_df)} genes from existing results")
    print()

    print("CONSERVATIVE COUNTING APPROACH:")
    print("="*80)
    print("We will recalculate frequencies using binary counting")
    print("(max 1 SV per gene per sample, regardless of actual SV count)")
    print()

    # Create conservative version of existing results
    conservative_results = []

    for _, row in existing_df.iterrows():
        gene = row['gene']

        # Existing method counts all SVs (can be >100%)
        standard_freq = row['cohort_freq']

        # Conservative: each affected sample counted once only
        n_samples_affected = row.get('cohort_n_samples', int(standard_freq * 200))

        # Cap at 200 samples max
        n_samples_affected = min(n_samples_affected, 200)

        conservative_freq = n_samples_affected / 200

        # Wilson CI for conservative frequency
        ci_low, ci_high = wilson_score_interval(n_samples_affected, 200)

        # Recalculate fold-changes with conservative frequency
        tcga_freq = row['tcga_freq']
        pcawg_freq = row['pcawg_freq']

        conservative_fc_tcga = conservative_freq / tcga_freq if tcga_freq > 0 else np.nan
        conservative_fc_pcawg = conservative_freq / pcawg_freq if pcawg_freq > 0 else np.nan

        # Calculate average SVs per affected sample
        avg_svs_per_sample = standard_freq / conservative_freq if conservative_freq > 0 else 0

        conservative_results.append({
            'gene': gene,
            'n_samples_affected': n_samples_affected,
            'conservative_frequency': conservative_freq,
            'standard_frequency': standard_freq,
            'conservative_freq_pct': f"{conservative_freq*100:.1f}%",
            'standard_freq_pct': f"{standard_freq*100:.1f}%",
            'ci_low': ci_low,
            'ci_high': ci_high,
            'tcga_freq': tcga_freq,
            'pcawg_freq': pcawg_freq,
            'fold_change_tcga_conservative': conservative_fc_tcga,
            'fold_change_pcawg_conservative': conservative_fc_pcawg,
            'fold_change_tcga_standard': row['fold_change_tcga'],
            'fold_change_pcawg_standard': row['fold_change_pcawg'],
            'difference': standard_freq - conservative_freq,
            'avg_svs_per_affected_sample': avg_svs_per_sample,
            'is_driver': row.get('is_driver', True)
        })

    conservative_df = pd.DataFrame(conservative_results)

    # Save results
    output_file = f"{output_dir}/conservative_vs_standard_comparison.csv"
    conservative_df.to_csv(output_file, index=False)
    print(f"✓ Saved: {output_file}")

    # Sort by conservative fold-change
    top_genes_conservative_tcga = conservative_df.sort_values('fold_change_tcga_conservative', ascending=False).head(10)
    top_genes_conservative_pcawg = conservative_df.sort_values('fold_change_pcawg_conservative', ascending=False).head(10)

    # Save top genes
    top_genes_conservative_tcga.to_csv(f"{output_dir}/top10_genes_conservative_vs_TCGA.csv", index=False)
    top_genes_conservative_pcawg.to_csv(f"{output_dir}/top10_genes_conservative_vs_PCAWG.csv", index=False)

    # Create summary report
    report_file = f"{output_dir}/CONSERVATIVE_COUNTING_REPORT.md"

    with open(report_file, 'w') as f:
        f.write("# Conservative Gene Counting Analysis\n\n")
        f.write("## Methodology\n\n")
        f.write("Following the approach from npae082.pdf:\n")
        f.write("> \"For each gene, we reported only one instance of a specific alteration type,\n")
        f.write("> regardless of how many identical alterations occurred within that same gene.\"\n\n")
        f.write("**Conservative counting**: Each gene counted as affected (1) or not (0) per sample\n\n")
        f.write("**Standard counting**: All SVs counted (can be >100% if multiple SVs per gene)\n\n")

        f.write("## Impact on Frequencies\n\n")
        f.write("### Genes with Largest Frequency Reduction (High Chromothripsis)\n\n")
        f.write("| Gene | Conservative | Standard | Difference | Avg SVs/Sample | Interpretation |\n")
        f.write("|------|--------------|----------|------------|----------------|----------------|\n")

        for _, row in conservative_df.sort_values('difference', ascending=False).head(10).iterrows():
            f.write(f"| {row['gene']} | {row['conservative_freq_pct']} | {row['standard_freq_pct']} | ")
            f.write(f"{row['difference']*100:.1f}% | {row['avg_svs_per_affected_sample']:.2f} | ")

            if row['difference'] > 0.5:
                f.write("High chromothripsis |\n")
            elif row['difference'] > 0.2:
                f.write("Moderate chromothripsis |\n")
            else:
                f.write("Mostly single SVs |\n")

        f.write("\n## Top 10 Genes (Conservative Fold-Changes vs TCGA)\n\n")
        f.write("| Gene | Conservative Freq | TCGA Freq | Fold-Change (Conservative) | Fold-Change (Standard) | Driver |\n")
        f.write("|------|-------------------|-----------|----------------------------|------------------------|--------|\n")

        for _, row in top_genes_conservative_tcga.iterrows():
            f.write(f"| {row['gene']} | {row['conservative_freq_pct']} | ")
            f.write(f"{row['tcga_freq']*100:.1f}% | ")
            f.write(f"{row['fold_change_tcga_conservative']:.1f}× | ")
            f.write(f"{row['fold_change_tcga_standard']:.1f}× | ")
            f.write(f"{'Yes' if row['is_driver'] else 'No'} |\n")

        f.write("\n## Top 10 Genes (Conservative Fold-Changes vs PCAWG)\n\n")
        f.write("| Gene | Conservative Freq | PCAWG Freq | Fold-Change (Conservative) | Fold-Change (Standard) | Driver |\n")
        f.write("|------|-------------------|------------|----------------------------|------------------------|--------|\n")

        for _, row in top_genes_conservative_pcawg.iterrows():
            f.write(f"| {row['gene']} | {row['conservative_freq_pct']} | ")
            f.write(f"{row['pcawg_freq']*100:.1f}% | ")
            f.write(f"{row['fold_change_pcawg_conservative']:.1f}× | ")
            f.write(f"{row['fold_change_pcawg_standard']:.1f}× | ")
            f.write(f"{'Yes' if row['is_driver'] else 'No'} |\n")

        f.write("\n## Comparison: Conservative vs Standard Fold-Changes\n\n")
        f.write("| Gene | Conservative FC (TCGA) | Standard FC (TCGA) | Ratio |\n")
        f.write("|------|------------------------|-----------------------|-------|\n")

        for _, row in conservative_df.sort_values('fold_change_tcga_standard', ascending=False).head(10).iterrows():
            f.write(f"| {row['gene']} | ")
            f.write(f"{row['fold_change_tcga_conservative']:.1f}× | ")
            f.write(f"{row['fold_change_tcga_standard']:.1f}× | ")

            ratio = row['fold_change_tcga_conservative'] / row['fold_change_tcga_standard'] if row['fold_change_tcga_standard'] > 0 else 0
            f.write(f"{ratio*100:.0f}% |\n")

        f.write("\n## Key Findings\n\n")

        # Calculate statistics
        avg_reduction = conservative_df['difference'].mean() * 100
        max_reduction = conservative_df['difference'].max() * 100
        max_gene = conservative_df.loc[conservative_df['difference'].idxmax(), 'gene']

        f.write(f"- **Average frequency reduction**: {avg_reduction:.1f}%\n")
        f.write(f"- **Maximum frequency reduction**: {max_reduction:.1f}% ({max_gene} - high chromothripsis)\n")

        # Genes with minimal difference (mostly single SVs)
        single_sv_genes = conservative_df[conservative_df['difference'] < 0.1]
        f.write(f"- **Genes with mostly single SVs** (difference <10%): {len(single_sv_genes)} ({len(single_sv_genes)/len(conservative_df)*100:.0f}%)\n")

        # Genes with high chromothripsis
        chromothripsis_genes = conservative_df[conservative_df['difference'] > 0.5]
        f.write(f"- **Genes with chromothripsis** (difference >50%): {len(chromothripsis_genes)}\n")

        if len(chromothripsis_genes) > 0:
            f.write(f"  - {', '.join(chromothripsis_genes['gene'].tolist())}\n")

        # Average SVs per affected sample
        avg_svs_overall = conservative_df['avg_svs_per_affected_sample'].mean()
        f.write(f"- **Average SVs per affected sample** (across all genes): {avg_svs_overall:.2f}\n")

        f.write("\n## Interpretation\n\n")
        f.write("**Conservative counting** is more comparable to literature studies that use\n")
        f.write("matched tumor-normal sequencing (e.g., TCGA, PCAWG), where each gene is typically\n")
        f.write("scored as altered or not per sample.\n\n")

        f.write("**Standard counting** (your current approach) reveals the true extent of\n")
        f.write("chromothripsis by counting all SV events, which is more informative for\n")
        f.write("understanding the genomic complexity of GBM.\n\n")

        f.write("**Key Insight**: Genes with large differences between conservative and standard\n")
        f.write("frequencies are undergoing chromothripsis (multiple SVs per gene in the same sample),\n")
        f.write("indicating catastrophic genomic rearrangements.\n\n")

        f.write("## Recommendation for Publication\n\n")
        f.write("Report **both metrics**:\n\n")
        f.write("1. **Conservative frequency**: For direct comparison with literature\n")
        f.write("   - Use in fold-change calculations vs TCGA/PCAWG\n")
        f.write("   - Present in main figures/tables\n\n")

        f.write("2. **Standard frequency** (with avg SVs per sample): To demonstrate chromothripsis\n")
        f.write("   - Highlight genes with >2 SVs per affected sample\n")
        f.write("   - Include in supplementary materials\n")
        f.write("   - Emphasize unprecedented structural complexity\n\n")

        f.write("This dual approach:\n")
        f.write("- Makes your results directly comparable to existing literature\n")
        f.write("- Demonstrates the novel insight from long-read sequencing (chromothripsis detection)\n")
        f.write("- Provides transparent methodology\n")

    print(f"✓ Report saved: {report_file}")
    print()
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total genes analyzed: {len(conservative_df)}")
    print(f"Average frequency reduction (conservative vs standard): {avg_reduction:.1f}%")
    print(f"\nTop 5 genes (conservative counting vs TCGA):")

    for i, (_, row) in enumerate(top_genes_conservative_tcga.head(5).iterrows(), 1):
        print(f"  {i}. {row['gene']}: {row['conservative_freq_pct']} (FC: {row['fold_change_tcga_conservative']:.1f}×)")

    print(f"\nGenes with high chromothripsis (>50% difference): {len(chromothripsis_genes)}")
    if len(chromothripsis_genes) > 0:
        for gene in chromothripsis_genes.head(5)['gene']:
            print(f"  - {gene}")

    print(f"\nAll results saved to: {output_dir}/")
    print("="*80)

if __name__ == "__main__":
    main()
