#!/usr/bin/env python3
"""
Create quality metrics comparison table for matching variants
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

def get_variant_key(row):
    """Create unique variant identifier"""
    return f"{row['Chr']}:{row['Start']}-{row['End']}:{row['Ref']}>{row['Alt']}"

def extract_quality_metrics(df, caller_name):
    """Extract quality metrics from Otherinfo10 column"""
    if df is None or len(df) == 0:
        return {}

    variants = {}

    for idx, row in df.iterrows():
        key = get_variant_key(row)
        gene = row.get('Gene.refGene', 'Unknown')

        try:
            info = str(row['Otherinfo10'])
            parts = info.split(':')

            gt = parts[0] if len(parts) > 0 else ''
            gq = parts[1] if len(parts) > 1 else ''
            depth = parts[2] if len(parts) > 2 else ''

            if caller_name == 'DeepSomatic':
                # GT:GQ:DP:AD:VAF:PL
                ad = parts[3] if len(parts) > 3 else ''
                af = parts[4] if len(parts) > 4 else ''
            else:  # ClairS-TO
                # GT:GQ:DP:AF:AD:...
                af = parts[3] if len(parts) > 3 else ''
                ad = parts[4] if len(parts) > 4 else ''

            variants[key] = {
                'gene': gene,
                'chr': row['Chr'],
                'pos': row['Start'],
                'ref': row['Ref'],
                'alt': row['Alt'],
                'gt': gt,
                'gq': gq,
                'depth': depth,
                'ad': ad,
                'af': af
            }
        except:
            continue

    return variants

def create_quality_comparison():
    benchmark_dir = Path("/home/chbope/extension/script/deepsomatic/benchmark")
    ds_dir = benchmark_dir / "deepsomatic"
    clair_dir = benchmark_dir / "clairsto"
    output_dir = benchmark_dir / "comparison_results"

    # Load summary data to get sample list
    summary_df = pd.read_csv(output_dir / "per_sample_comparison.csv")

    all_results = []

    for _, row in summary_df.iterrows():
        sample_id = row['sample_id']
        print(f"Processing {sample_id}...")

        ds_df, clair_df = load_sample_data(sample_id, ds_dir, clair_dir)

        # Extract quality metrics
        ds_variants = extract_quality_metrics(ds_df, 'DeepSomatic')
        clair_variants = extract_quality_metrics(clair_df, 'ClairS-TO')

        # Find common variants
        ds_keys = set(ds_variants.keys())
        clair_keys = set(clair_variants.keys())
        common_keys = ds_keys & clair_keys

        print(f"  Common variants: {len(common_keys)}")

        # Compare quality metrics for common variants
        for key in sorted(common_keys):
            ds_var = ds_variants[key]
            clair_var = clair_variants[key]

            # Calculate differences
            try:
                gq_diff = abs(int(ds_var['gq']) - int(clair_var['gq'])) if ds_var['gq'] and clair_var['gq'] else None
                depth_diff = abs(int(ds_var['depth']) - int(clair_var['depth'])) if ds_var['depth'] and clair_var['depth'] else None
                af_diff = abs(float(ds_var['af']) - float(clair_var['af'])) if ds_var['af'] and clair_var['af'] else None
            except:
                gq_diff = None
                depth_diff = None
                af_diff = None

            all_results.append({
                'sample_id': sample_id,
                'gene': ds_var['gene'],
                'variant': f"{ds_var['chr']}:{ds_var['pos']} {ds_var['ref']}>{ds_var['alt']}",
                'ds_gt': ds_var['gt'],
                'clair_gt': clair_var['gt'],
                'gt_match': ds_var['gt'] == clair_var['gt'],
                'ds_gq': ds_var['gq'],
                'clair_gq': clair_var['gq'],
                'gq_diff': gq_diff,
                'ds_depth': ds_var['depth'],
                'clair_depth': clair_var['depth'],
                'depth_diff': depth_diff,
                'ds_ad': ds_var['ad'],
                'clair_ad': clair_var['ad'],
                'ds_af': ds_var['af'],
                'clair_af': clair_var['af'],
                'af_diff': round(af_diff, 4) if af_diff is not None else None
            })

    # Create DataFrame
    result_df = pd.DataFrame(all_results)

    if len(result_df) > 0:
        # Save to CSV
        output_file = output_dir / "quality_metrics_comparison.csv"
        result_df.to_csv(output_file, index=False)

        print(f"\n✓ Quality metrics comparison saved to: {output_file}")
        print(f"  Total common variants compared: {len(result_df)}")

        # Calculate statistics
        print("\n" + "="*80)
        print("QUALITY METRICS COMPARISON SUMMARY")
        print("="*80)

        if 'gt_match' in result_df.columns:
            gt_match_pct = (result_df['gt_match'].sum() / len(result_df)) * 100
            print(f"\nGenotype (GT) Agreement: {gt_match_pct:.1f}% ({result_df['gt_match'].sum}/{len(result_df)})")

        if 'gq_diff' in result_df.columns and result_df['gq_diff'].notna().any():
            avg_gq_diff = result_df['gq_diff'].mean()
            max_gq_diff = result_df['gq_diff'].max()
            print(f"\nGenotype Quality (GQ) Difference:")
            print(f"  Average: {avg_gq_diff:.1f}")
            print(f"  Maximum: {max_gq_diff:.0f}")

        if 'depth_diff' in result_df.columns and result_df['depth_diff'].notna().any():
            avg_depth_diff = result_df['depth_diff'].mean()
            max_depth_diff = result_df['depth_diff'].max()
            print(f"\nDepth (DP) Difference:")
            print(f"  Average: {avg_depth_diff:.1f}")
            print(f"  Maximum: {max_depth_diff:.0f}")

        if 'af_diff' in result_df.columns and result_df['af_diff'].notna().any():
            avg_af_diff = result_df['af_diff'].mean()
            max_af_diff = result_df['af_diff'].max()
            print(f"\nAllele Frequency (AF/VAF) Difference:")
            print(f"  Average: {avg_af_diff:.4f} ({avg_af_diff*100:.2f}%)")
            print(f"  Maximum: {max_af_diff:.4f} ({max_af_diff*100:.2f}%)")

        # Show examples with largest differences
        print("\n" + "="*80)
        print("TOP 5 VARIANTS WITH LARGEST AF DIFFERENCES")
        print("="*80)

        if 'af_diff' in result_df.columns:
            top_af_diff = result_df.nlargest(5, 'af_diff')[
                ['sample_id', 'gene', 'variant', 'ds_af', 'clair_af', 'af_diff', 'ds_depth', 'clair_depth']
            ]
            print(top_af_diff.to_string(index=False))

    else:
        print("\nNo common variants found for comparison!")

    return 0

if __name__ == "__main__":
    sys.exit(create_quality_comparison())
