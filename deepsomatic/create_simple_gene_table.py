#!/usr/bin/env python3
"""
Create simple gene comparison table CSV
"""

import pandas as pd
import sys
from pathlib import Path

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
        return set()

    genes = set()
    for gene_str in df['Gene.refGene'].dropna():
        genes.update([g.strip() for g in str(gene_str).split(',')])
    return genes

def get_vaf_info(df, caller_name=''):
    """Get VAF range information"""
    if df is None or len(df) == 0:
        return "N/A"

    vafs = []

    # Check if VAF column exists
    if 'VAF' in df.columns:
        for vaf in df['VAF']:
            try:
                if isinstance(vaf, str):
                    vaf = float(vaf)
                vafs.append(vaf)
            except:
                continue
    # Extract from Otherinfo10 if VAF column doesn't exist
    elif 'Otherinfo10' in df.columns:
        for info in df['Otherinfo10']:
            try:
                # Format: GT:GQ:DP:AD:VAF:PL for DeepSomatic
                # Format: GT:GQ:DP:AF:AD:... for ClairS-TO (AF is at position 3)
                parts = str(info).split(':')
                if caller_name == 'ClairS-TO' and len(parts) > 3:
                    vaf = float(parts[3])
                elif caller_name == 'DeepSomatic' and len(parts) > 4:
                    vaf = float(parts[4])
                else:
                    continue
                vafs.append(vaf)
            except:
                continue

    if not vafs:
        return "N/A"

    min_vaf = min(vafs)
    max_vaf = max(vafs)

    if min_vaf == max_vaf:
        return f"{min_vaf:.1%}"
    else:
        return f"{min_vaf:.1%}-{max_vaf:.1%}"

def create_simple_gene_table():
    benchmark_dir = Path("/home/chbope/extension/script/deepsomatic/benchmark")
    ds_dir = benchmark_dir / "deepsomatic"
    clair_dir = benchmark_dir / "clairsto"
    output_dir = benchmark_dir / "comparison_results"

    # Load summary data to get sample list
    summary_df = pd.read_csv(output_dir / "per_sample_comparison.csv")

    results = []

    for _, row in summary_df.iterrows():
        sample_id = row['sample_id']
        print(f"Processing {sample_id}...")

        ds_df, clair_df = load_sample_data(sample_id, ds_dir, clair_dir)

        # Extract genes
        ds_genes = extract_genes(ds_df)
        clair_genes = extract_genes(clair_df)

        # Find common and unique genes
        common_genes = ds_genes & clair_genes
        ds_only_genes = ds_genes - clair_genes
        clair_only_genes = clair_genes - ds_genes

        results.append({
            'sample_id': sample_id,
            'deepsomatic_genes': ', '.join(sorted(ds_genes)) if ds_genes else 'None',
            'clairsto_genes': ', '.join(sorted(clair_genes)) if clair_genes else 'None',
            'common_genes': ', '.join(sorted(common_genes)) if common_genes else 'None',
            'missing_in_clairsto': ', '.join(sorted(ds_only_genes)) if ds_only_genes else 'None',
            'missing_in_deepsomatic': ', '.join(sorted(clair_only_genes)) if clair_only_genes else 'None'
        })

    # Create DataFrame
    result_df = pd.DataFrame(results)

    # Save to CSV
    output_file = output_dir / "gene_comparison_simple.csv"
    result_df.to_csv(output_file, index=False)

    print(f"\n✓ Simple gene comparison table saved to: {output_file}")
    print(f"  Total samples: {len(result_df)}")

    # Also print to screen
    print("\n" + "="*120)
    print("GENE COMPARISON TABLE")
    print("="*120)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', 80)
    pd.set_option('display.width', 200)
    print(result_df.to_string(index=False))

    return 0

if __name__ == "__main__":
    sys.exit(create_simple_gene_table())
