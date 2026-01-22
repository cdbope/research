#!/usr/bin/env python3
"""
Comprehensive Benchmark Comparison: ClairS-TO vs DeepSomatic
Analyzes multiple samples to determine overall performance
"""

import pandas as pd
import sys
from pathlib import Path
import glob
import numpy as np

def extract_sample_id(filename):
    """Extract sample ID from filename"""
    # Extract pattern like KM20-342, T18-020, etc.
    parts = Path(filename).stem.split('_')[0]
    return parts

def load_annotated_csv(file_path, caller_name):
    """Load annotated CSV file"""
    try:
        df = pd.read_csv(file_path, sep='\t')

        # Create unique variant ID
        df['Variant_ID'] = df['Chr'].astype(str) + ':' + \
                           df['Start'].astype(str) + ':' + \
                           df['Ref'] + '>' + df['Alt']

        # Extract VAF if available
        if caller_name == 'DeepSomatic' and 'Otherinfo10' in df.columns:
            # DeepSomatic format: GT:GQ:DP:AD:VAF:PL
            df['VAF'] = df['Otherinfo10'].str.split(':').str[4].astype(float)
            df['Depth'] = df['Otherinfo10'].str.split(':').str[2].astype(int)
        elif caller_name == 'ClairS-TO' and 'Otherinfo10' in df.columns:
            # ClairS-TO format: GT:GQ:DP:AF:AD:...
            df['VAF'] = df['Otherinfo10'].str.split(':').str[3].astype(float)
            df['Depth'] = df['Otherinfo10'].str.split(':').str[2].astype(int)
        elif 'AF' in df.columns:  # Alternative format
            df['VAF'] = df['AF']
            df['Depth'] = df['Depth'] if 'Depth' in df.columns else None

        df['Caller'] = caller_name
        return df
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def analyze_sample(sample_id, ds_file, clair_file):
    """Analyze a single sample"""

    # Load data
    ds_df = load_annotated_csv(ds_file, 'DeepSomatic')
    clair_df = load_annotated_csv(clair_file, 'ClairS-TO')

    if ds_df is None or clair_df is None:
        return None

    # Get variant sets
    ds_variants = set(ds_df['Variant_ID'])
    clair_variants = set(clair_df['Variant_ID'])

    # Calculate metrics
    shared = ds_variants & clair_variants
    ds_only = ds_variants - clair_variants
    clair_only = clair_variants - ds_variants
    total_unique = len(ds_variants | clair_variants)

    # Concordance
    concordance = len(shared) / total_unique * 100 if total_unique > 0 else 0

    # Quality metrics
    cancer_genes = {
        'TP53', 'KRAS', 'EGFR', 'PIK3CA', 'BRAF', 'PTEN', 'ALK', 'ROS1',
        'ERBB2', 'MET', 'FGFR1', 'FGFR2', 'FGFR3', 'PDGFRA', 'KIT', 'RET',
        'NRAS', 'HRAS', 'JAK2', 'IDH1', 'IDH2', 'NOTCH1', 'NOTCH2',
        'H3-3A', 'H3-3B', 'TERT', 'ATRX', 'DAXX', 'SETD2', 'ARID1A',
        'SMARCA4', 'PTEN', 'PPM1D', 'HDAC9', 'KDM4C'
    }

    # Count cancer genes
    ds_cancer_genes = set()
    clair_cancer_genes = set()

    if 'Gene.refGene' in ds_df.columns:
        for gene_str in ds_df['Gene.refGene'].dropna():
            genes = str(gene_str).split(',')
            ds_cancer_genes.update(g.strip() for g in genes if g.strip() in cancer_genes)

    if 'Gene.refGene' in clair_df.columns:
        for gene_str in clair_df['Gene.refGene'].dropna():
            genes = str(gene_str).split(',')
            clair_cancer_genes.update(g.strip() for g in genes if g.strip() in cancer_genes)

    # Clinical significance
    ds_pathogenic = 0
    clair_pathogenic = 0

    if 'CLNSIG' in ds_df.columns:
        ds_pathogenic = ds_df['CLNSIG'].astype(str).str.contains('athogenic', case=False, na=False).sum()

    if 'CLNSIG' in clair_df.columns:
        clair_pathogenic = clair_df['CLNSIG'].astype(str).str.contains('athogenic', case=False, na=False).sum()

    # COSMIC hits
    ds_cosmic = ds_df['COSMIC100'].notna().sum() if 'COSMIC100' in ds_df.columns else 0
    clair_cosmic = clair_df['COSMIC100'].notna().sum() if 'COSMIC100' in clair_df.columns else 0

    # VAF analysis
    ds_min_vaf = ds_df['VAF'].min() if 'VAF' in ds_df.columns and len(ds_df) > 0 else np.nan
    clair_min_vaf = clair_df['VAF'].min() if 'VAF' in clair_df.columns and len(clair_df) > 0 else np.nan

    ds_mean_vaf = ds_df['VAF'].mean() if 'VAF' in ds_df.columns and len(ds_df) > 0 else np.nan
    clair_mean_vaf = clair_df['VAF'].mean() if 'VAF' in clair_df.columns and len(clair_df) > 0 else np.nan

    # Low VAF variants
    ds_low_vaf = ((ds_df['VAF'] < 0.30).sum() if 'VAF' in ds_df.columns else 0)
    clair_low_vaf = ((clair_df['VAF'] < 0.30).sum() if 'VAF' in clair_df.columns else 0)

    return {
        'sample_id': sample_id,
        'ds_variants': len(ds_variants),
        'clair_variants': len(clair_variants),
        'shared_variants': len(shared),
        'ds_only': len(ds_only),
        'clair_only': len(clair_only),
        'total_unique': total_unique,
        'concordance': concordance,
        'ds_cancer_genes': len(ds_cancer_genes),
        'clair_cancer_genes': len(clair_cancer_genes),
        'ds_pathogenic': ds_pathogenic,
        'clair_pathogenic': clair_pathogenic,
        'ds_cosmic': ds_cosmic,
        'clair_cosmic': clair_cosmic,
        'ds_min_vaf': ds_min_vaf,
        'clair_min_vaf': clair_min_vaf,
        'ds_mean_vaf': ds_mean_vaf,
        'clair_mean_vaf': clair_mean_vaf,
        'ds_low_vaf': ds_low_vaf,
        'clair_low_vaf': clair_low_vaf
    }

def main():
    benchmark_dir = Path("/home/chbope/extension/script/deepsomatic/benchmark")
    ds_dir = benchmark_dir / "deepsomatic"
    clair_dir = benchmark_dir / "clairsto"
    output_dir = benchmark_dir / "comparison_results"

    output_dir.mkdir(exist_ok=True)

    print("="*80)
    print("BENCHMARK COMPARISON: ClairS-TO vs DeepSomatic")
    print("="*80)
    print()

    # Find all samples
    ds_files = sorted(glob.glob(str(ds_dir / "*_annotateandfilter_deep_somatic.csv")))

    if not ds_files:
        print("Error: No DeepSomatic files found!")
        return 1

    print(f"Found {len(ds_files)} samples to analyze")
    print()

    results = []

    # Analyze each sample
    for ds_file in ds_files:
        sample_id = extract_sample_id(ds_file)

        # Find corresponding ClairS-TO file
        clair_file = clair_dir / f"{sample_id}_annotateandfilter_clairsto.csv"

        if not clair_file.exists():
            print(f"Warning: No ClairS-TO file for {sample_id}, skipping...")
            continue

        print(f"Analyzing {sample_id}...", end=' ')

        result = analyze_sample(sample_id, ds_file, str(clair_file))

        if result:
            results.append(result)
            print(f"✓ DS:{result['ds_variants']} Clair:{result['clair_variants']} Concordance:{result['concordance']:.1f}%")
        else:
            print("✗ Failed")

    if not results:
        print("Error: No samples analyzed successfully")
        return 1

    # Create DataFrame
    df = pd.DataFrame(results)

    # Save detailed results
    df.to_csv(output_dir / "per_sample_comparison.csv", index=False)
    print(f"\nDetailed results saved to: {output_dir / 'per_sample_comparison.csv'}")

    # Calculate aggregate statistics
    print("\n" + "="*80)
    print("AGGREGATE STATISTICS ACROSS ALL SAMPLES")
    print("="*80)

    print(f"\nTotal samples analyzed: {len(df)}")
    print(f"\nTotal variants:")
    print(f"  DeepSomatic:   {df['ds_variants'].sum()}")
    print(f"  ClairS-TO:     {df['clair_variants'].sum()}")
    print(f"  Shared:        {df['shared_variants'].sum()}")
    print(f"  DS only:       {df['ds_only'].sum()}")
    print(f"  Clair only:    {df['clair_only'].sum()}")

    print(f"\nAverage per sample:")
    print(f"  DeepSomatic:   {df['ds_variants'].mean():.1f} ± {df['ds_variants'].std():.1f}")
    print(f"  ClairS-TO:     {df['clair_variants'].mean():.1f} ± {df['clair_variants'].std():.1f}")
    print(f"  Concordance:   {df['concordance'].mean():.1f}% ± {df['concordance'].std():.1f}%")

    print(f"\nConcordance range:")
    print(f"  Minimum:       {df['concordance'].min():.1f}%")
    print(f"  Maximum:       {df['concordance'].max():.1f}%")
    print(f"  Median:        {df['concordance'].median():.1f}%")

    print(f"\nCancer genes detected:")
    print(f"  DeepSomatic:   {df['ds_cancer_genes'].sum()} (avg: {df['ds_cancer_genes'].mean():.1f}/sample)")
    print(f"  ClairS-TO:     {df['clair_cancer_genes'].sum()} (avg: {df['clair_cancer_genes'].mean():.1f}/sample)")

    print(f"\nPathogenic variants:")
    print(f"  DeepSomatic:   {df['ds_pathogenic'].sum()} total")
    print(f"  ClairS-TO:     {df['clair_pathogenic'].sum()} total")

    print(f"\nCOSMIC variants:")
    print(f"  DeepSomatic:   {df['ds_cosmic'].sum()} total")
    print(f"  ClairS-TO:     {df['clair_cosmic'].sum()} total")

    print(f"\nVAF Analysis:")
    print(f"  DeepSomatic minimum VAF:   {df['ds_min_vaf'].min():.4f} (best), {df['ds_min_vaf'].mean():.4f} (avg)")
    print(f"  ClairS-TO minimum VAF:     {df['clair_min_vaf'].min():.4f} (best), {df['clair_min_vaf'].mean():.4f} (avg)")
    print(f"  DeepSomatic mean VAF:      {df['ds_mean_vaf'].mean():.4f}")
    print(f"  ClairS-TO mean VAF:        {df['clair_mean_vaf'].mean():.4f}")

    print(f"\nLow VAF variants (<30%):")
    print(f"  DeepSomatic:   {df['ds_low_vaf'].sum()} total ({df['ds_low_vaf'].mean():.1f}/sample)")
    print(f"  ClairS-TO:     {df['clair_low_vaf'].sum()} total ({df['clair_low_vaf'].mean():.1f}/sample)")

    # Winner determination
    print("\n" + "="*80)
    print("PERFORMANCE COMPARISON")
    print("="*80)

    metrics = {
        'Total variants detected': (df['ds_variants'].sum(), df['clair_variants'].sum()),
        'Avg variants per sample': (df['ds_variants'].mean(), df['clair_variants'].mean()),
        'Cancer genes detected': (df['ds_cancer_genes'].sum(), df['clair_cancer_genes'].sum()),
        'Pathogenic variants': (df['ds_pathogenic'].sum(), df['clair_pathogenic'].sum()),
        'COSMIC hits': (df['ds_cosmic'].sum(), df['clair_cosmic'].sum()),
        'Low VAF detection': (df['ds_low_vaf'].sum(), df['clair_low_vaf'].sum()),
        'Best minimum VAF': (df['ds_min_vaf'].min(), df['clair_min_vaf'].min())
    }

    ds_wins = 0
    clair_wins = 0

    print(f"\n{'Metric':<35} {'DeepSomatic':<15} {'ClairS-TO':<15} {'Winner'}")
    print("-"*80)

    for metric, (ds_val, clair_val) in metrics.items():
        if metric == 'Best minimum VAF':
            winner = 'ClairS-TO ✓' if clair_val < ds_val else 'DeepSomatic ✓'
            if clair_val < ds_val:
                clair_wins += 1
            else:
                ds_wins += 1
        else:
            winner = 'ClairS-TO ✓' if clair_val > ds_val else 'DeepSomatic ✓'
            if clair_val > ds_val:
                clair_wins += 1
            else:
                ds_wins += 1

        if isinstance(ds_val, float):
            print(f"{metric:<35} {ds_val:<15.2f} {clair_val:<15.2f} {winner}")
        else:
            print(f"{metric:<35} {ds_val:<15} {clair_val:<15} {winner}")

    print(f"\n{'='*80}")
    print(f"FINAL SCORE: ClairS-TO {clair_wins} - {ds_wins} DeepSomatic")
    print(f"{'='*80}")

    if clair_wins > ds_wins:
        winner_name = "ClairS-TO"
        winner_score = clair_wins
        loser_score = ds_wins
    else:
        winner_name = "DeepSomatic"
        winner_score = ds_wins
        loser_score = clair_wins

    print(f"\n🏆 WINNER: {winner_name} ({winner_score}/{len(metrics)} metrics)")
    print(f"\nConclusion:")
    print(f"  Across {len(df)} samples, {winner_name} demonstrates better performance")
    print(f"  in {winner_score} out of {len(metrics)} key metrics.")

    # Save summary
    with open(output_dir / "benchmark_summary.txt", 'w') as f:
        f.write("BENCHMARK COMPARISON SUMMARY\n")
        f.write("="*80 + "\n\n")
        f.write(f"Samples analyzed: {len(df)}\n")
        f.write(f"Winner: {winner_name}\n")
        f.write(f"Score: {winner_score}-{loser_score}\n\n")
        f.write(f"DeepSomatic total variants: {df['ds_variants'].sum()}\n")
        f.write(f"ClairS-TO total variants: {df['clair_variants'].sum()}\n")
        f.write(f"Average concordance: {df['concordance'].mean():.1f}%\n")

    print(f"\nSummary saved to: {output_dir / 'benchmark_summary.txt'}")
    print("="*80)

    return 0

if __name__ == "__main__":
    sys.exit(main())
