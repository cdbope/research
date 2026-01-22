#!/usr/bin/env python3
"""
Variant Caller Comparison Tool
Compare variants called by DeepSomatic vs Clair3/ClairS-TO
"""

import pandas as pd
import sys
from pathlib import Path

def load_deepsomatic(file_path):
    """Load DeepSomatic annotated variants"""
    print(f"Loading DeepSomatic variants from: {file_path}")
    df = pd.read_csv(file_path, sep='\t')

    # Create unique variant ID
    df['Variant_ID'] = df['Chr'].astype(str) + ':' + \
                       df['Start'].astype(str) + ':' + \
                       df['Ref'] + '>' + df['Alt']

    # Extract VAF from Otherinfo10 (format: GT:GQ:DP:AD:VAF:PL)
    if 'Otherinfo10' in df.columns:
        df['VAF_DS'] = df['Otherinfo10'].str.split(':').str[4].astype(float)
        df['Depth_DS'] = df['Otherinfo10'].str.split(':').str[2].astype(int)

    df['Caller'] = 'DeepSomatic'
    print(f"  Found {len(df)} variants")
    return df

def load_clair(file_path):
    """Load Clair3/ClairS-TO annotated variants"""
    print(f"Loading Clair3/ClairS-TO variants from: {file_path}")
    df = pd.read_csv(file_path, sep='\t')

    # Create unique variant ID
    df['Variant_ID'] = df['Chr'].astype(str) + ':' + \
                       df['Start'].astype(str) + ':' + \
                       df['Ref'] + '>' + df['Alt']

    # Rename columns for consistency
    df['VAF_Clair'] = df['AF'] if 'AF' in df.columns else None
    df['Depth_Clair'] = df['Depth'] if 'Depth' in df.columns else None
    df['Caller'] = 'Clair3/ClairS-TO'
    print(f"  Found {len(df)} variants")
    return df

def compare_variants(ds_df, clair_df, output_dir):
    """Compare variants between two callers"""
    print("\n" + "="*70)
    print("VARIANT COMPARISON ANALYSIS")
    print("="*70)

    # Get variant IDs
    ds_variants = set(ds_df['Variant_ID'])
    clair_variants = set(clair_df['Variant_ID'])

    # Calculate overlaps
    shared_variants = ds_variants & clair_variants
    ds_only = ds_variants - clair_variants
    clair_only = clair_variants - ds_variants

    total_unique = len(ds_variants | clair_variants)

    print(f"\n1. OVERALL STATISTICS")
    print(f"   {'DeepSomatic variants:':<30} {len(ds_variants):>6}")
    print(f"   {'Clair3/ClairS-TO variants:':<30} {len(clair_variants):>6}")
    print(f"   {'Total unique variants:':<30} {total_unique:>6}")
    print(f"   {'Shared variants:':<30} {len(shared_variants):>6} ({len(shared_variants)/total_unique*100:.1f}%)")
    print(f"   {'DeepSomatic only:':<30} {len(ds_only):>6} ({len(ds_only)/total_unique*100:.1f}%)")
    print(f"   {'Clair3/ClairS-TO only:':<30} {len(clair_only):>6} ({len(clair_only)/total_unique*100:.1f}%)")

    # Venn diagram data
    print(f"\n2. VENN DIAGRAM DATA")
    print(f"   DeepSomatic unique: {len(ds_only)}")
    print(f"   Shared: {len(shared_variants)}")
    print(f"   Clair3/ClairS-TO unique: {len(clair_only)}")

    # Concordance rate
    concordance = len(shared_variants) / total_unique * 100 if total_unique > 0 else 0
    print(f"\n3. CONCORDANCE RATE: {concordance:.2f}%")

    # Gene-level comparison
    print(f"\n4. GENE-LEVEL COMPARISON")
    ds_genes = set(ds_df['Gene.refGene'].dropna())
    clair_genes = set(clair_df['Gene.refGene'].dropna())
    shared_genes = ds_genes & clair_genes

    print(f"   {'Genes with variants (DeepSomatic):':<40} {len(ds_genes)}")
    print(f"   {'Genes with variants (Clair3/ClairS-TO):':<40} {len(clair_genes)}")
    print(f"   {'Shared genes:':<40} {len(shared_genes)}")

    # Functional impact comparison
    print(f"\n5. FUNCTIONAL IMPACT COMPARISON")
    for func_type in ['exonic', 'upstream', 'intronic']:
        ds_count = (ds_df['Func.refGene'] == func_type).sum()
        clair_count = (clair_df['Func.refGene'] == func_type).sum()
        print(f"   {func_type.capitalize():<20} DS: {ds_count:>4}   Clair: {clair_count:>4}")

    # Clinical significance comparison
    print(f"\n6. CLINICAL SIGNIFICANCE")
    if 'CLNSIG' in ds_df.columns and 'CLNSIG' in clair_df.columns:
        for sig in ['Pathogenic', 'Likely_pathogenic', 'Uncertain_significance']:
            ds_count = ds_df['CLNSIG'].astype(str).str.contains(sig, na=False).sum()
            clair_count = clair_df['CLNSIG'].astype(str).str.contains(sig, na=False).sum()
            print(f"   {sig.replace('_', ' '):<25} DS: {ds_count:>4}   Clair: {clair_count:>4}")

    # Create detailed comparison DataFrames
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Shared variants
    shared_df_ds = ds_df[ds_df['Variant_ID'].isin(shared_variants)].copy()
    shared_df_clair = clair_df[clair_df['Variant_ID'].isin(shared_variants)].copy()

    # Merge shared variants
    shared_merged = pd.merge(
        shared_df_ds[['Variant_ID', 'Chr', 'Start', 'Ref', 'Alt', 'Gene.refGene',
                       'ExonicFunc.refGene', 'AAChange.refGene', 'CLNSIG', 'COSMIC100',
                       'VAF_DS', 'Depth_DS']],
        shared_df_clair[['Variant_ID', 'VAF_Clair', 'Depth_Clair', 'Variant_caller']],
        on='Variant_ID',
        how='inner'
    )

    # Calculate VAF difference
    if 'VAF_DS' in shared_merged.columns and 'VAF_Clair' in shared_merged.columns:
        shared_merged['VAF_Diff'] = (shared_merged['VAF_DS'] - shared_merged['VAF_Clair']).abs()
        shared_merged['VAF_Diff_Pct'] = shared_merged['VAF_Diff'] * 100

    # Sort by VAF difference
    if 'VAF_Diff' in shared_merged.columns:
        shared_merged = shared_merged.sort_values('VAF_Diff', ascending=False)

    # DeepSomatic only variants
    ds_only_df = ds_df[ds_df['Variant_ID'].isin(ds_only)].copy()
    ds_only_df = ds_only_df.sort_values('VAF_DS', ascending=False)

    # Clair only variants
    clair_only_df = clair_df[clair_df['Variant_ID'].isin(clair_only)].copy()
    clair_only_df = clair_only_df.sort_values('VAF_Clair', ascending=False)

    # Save results
    shared_output = output_dir / "shared_variants.tsv"
    ds_only_output = output_dir / "deepsomatic_only.tsv"
    clair_only_output = output_dir / "clair_only.tsv"
    summary_output = output_dir / "comparison_summary.txt"

    shared_merged.to_csv(shared_output, sep='\t', index=False)
    ds_only_df.to_csv(ds_only_output, sep='\t', index=False)
    clair_only_df.to_csv(clair_only_output, sep='\t', index=False)

    print(f"\n7. OUTPUT FILES")
    print(f"   Shared variants: {shared_output}")
    print(f"   DeepSomatic only: {ds_only_output}")
    print(f"   Clair3/ClairS-TO only: {clair_only_output}")

    # Detailed analysis of shared variants
    if len(shared_merged) > 0:
        print(f"\n8. SHARED VARIANTS VAF ANALYSIS")
        if 'VAF_Diff' in shared_merged.columns:
            print(f"   Mean VAF difference: {shared_merged['VAF_Diff'].mean():.4f} ({shared_merged['VAF_Diff_Pct'].mean():.2f}%)")
            print(f"   Median VAF difference: {shared_merged['VAF_Diff'].median():.4f} ({shared_merged['VAF_Diff_Pct'].median():.2f}%)")
            print(f"   Max VAF difference: {shared_merged['VAF_Diff'].max():.4f} ({shared_merged['VAF_Diff_Pct'].max():.2f}%)")

            # VAF concordance (within 10%)
            concordant_vaf = (shared_merged['VAF_Diff'] < 0.10).sum()
            print(f"   VAF concordant (±10%): {concordant_vaf}/{len(shared_merged)} ({concordant_vaf/len(shared_merged)*100:.1f}%)")

    # Top discordant variants
    print(f"\n9. TOP 5 VARIANTS WITH LARGEST VAF DIFFERENCES")
    if len(shared_merged) > 0 and 'VAF_Diff' in shared_merged.columns:
        top_diff = shared_merged.head(5)
        print(f"   {'Variant':<30} {'Gene':<15} {'DS_VAF':<10} {'Clair_VAF':<10} {'Diff':<10}")
        print(f"   {'-'*80}")
        for _, row in top_diff.iterrows():
            variant = f"{row['Chr']}:{row['Start']}:{row['Ref']}>{row['Alt']}"
            gene = str(row['Gene.refGene'])[:14]
            print(f"   {variant:<30} {gene:<15} {row['VAF_DS']:<10.3f} {row['VAF_Clair']:<10.3f} {row['VAF_Diff']:<10.3f}")

    # Save summary to file
    with open(summary_output, 'w') as f:
        f.write("VARIANT CALLER COMPARISON SUMMARY\n")
        f.write("="*70 + "\n\n")
        f.write(f"DeepSomatic variants: {len(ds_variants)}\n")
        f.write(f"Clair3/ClairS-TO variants: {len(clair_variants)}\n")
        f.write(f"Shared variants: {len(shared_variants)}\n")
        f.write(f"DeepSomatic only: {len(ds_only)}\n")
        f.write(f"Clair3/ClairS-TO only: {len(clair_only)}\n")
        f.write(f"Concordance rate: {concordance:.2f}%\n")

    print(f"   Summary: {summary_output}")

    print(f"\n{'='*70}")
    print("Analysis complete!")
    print(f"{'='*70}\n")

    return {
        'shared': len(shared_variants),
        'ds_only': len(ds_only),
        'clair_only': len(clair_only),
        'concordance': concordance
    }

def main():
    # File paths - can be overridden by command line arguments
    if len(sys.argv) == 4:
        deepsomatic_file = sys.argv[1]
        clair_file = sys.argv[2]
        output_dir = sys.argv[3]
    else:
        # Default paths
        deepsomatic_file = "/home/chbope/extension/script/deepsomatic/output/T25-152_annotateandfilter_deep_somatic.csv"
        clair_file = "/home/chbope/extension/data/t25-152/merge_annot_clair3andclairsto/T25-152_merge_annotation_filter_snvs_allcall.csv"
        output_dir = "/home/chbope/extension/script/deepsomatic/comparison_results"

    # Check if files exist
    if not Path(deepsomatic_file).exists():
        print(f"Error: DeepSomatic file not found: {deepsomatic_file}")
        sys.exit(1)

    if not Path(clair_file).exists():
        print(f"Error: Clair3/ClairS-TO file not found: {clair_file}")
        sys.exit(1)

    # Load data
    ds_df = load_deepsomatic(deepsomatic_file)
    clair_df = load_clair(clair_file)

    # Compare
    results = compare_variants(ds_df, clair_df, output_dir)

    return 0

if __name__ == "__main__":
    sys.exit(main())
