#!/usr/bin/env python3
"""
Create detailed variant tables for each sample
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

def extract_variant_info(df, variant_caller):
    """Extract key variant information from dataframe"""
    if df is None or len(df) == 0:
        return pd.DataFrame()

    # Select key columns
    columns_to_keep = [
        'Gene.refGene',
        'ExonicFunc.refGene',
        'AAChange.refGene',
        'Chr',
        'Start',
        'End',
        'Ref',
        'Alt',
        'Func.refGene',
        'CLNSIG',
        'COSMIC100',
        'Otherinfo10'
    ]

    # Create simplified dataframe
    result = pd.DataFrame()

    result['Gene'] = df['Gene.refGene']
    result['Nucleotide_alteration'] = df['Chr'].astype(str) + ':' + df['Start'].astype(str) + df['Ref'] + '>' + df['Alt']
    result['Ref'] = df['Ref']
    result['Alt'] = df['Alt']
    result['Func'] = df['Func.refGene']
    result['CLNSIG'] = df['CLNSIG']
    result['COSMIC100'] = df['COSMIC100']

    # Extract GT, GQ, Depth, AD, AF/VAF from Otherinfo10
    gt_list = []
    gq_list = []
    depth_list = []
    ad_list = []
    af_list = []

    for info in df['Otherinfo10']:
        try:
            parts = str(info).split(':')
            gt = parts[0] if len(parts) > 0 else ''
            gq = parts[1] if len(parts) > 1 else ''
            depth = parts[2] if len(parts) > 2 else ''

            # Different format for DeepSomatic vs ClairS-TO
            if variant_caller == 'DeepSomatic':
                # GT:GQ:DP:AD:VAF:PL
                ad = parts[3] if len(parts) > 3 else ''
                af = parts[4] if len(parts) > 4 else ''
            else:  # ClairS-TO
                # GT:GQ:DP:AF:AD:...
                af = parts[3] if len(parts) > 3 else ''
                ad = parts[4] if len(parts) > 4 else ''

            gt_list.append(gt)
            gq_list.append(gq)
            depth_list.append(depth)
            ad_list.append(ad)
            af_list.append(af)
        except:
            gt_list.append('')
            gq_list.append('')
            depth_list.append('')
            ad_list.append('')
            af_list.append('')

    result['GQ'] = gq_list
    result['Depth'] = depth_list
    result['AD'] = ad_list
    result['GT'] = gt_list
    result['AF'] = af_list
    result['Chr'] = df['Chr']
    result['Start'] = df['Start']
    result['End'] = df['End']
    result['Variant_caller'] = variant_caller

    return result

def create_per_sample_tables():
    benchmark_dir = Path("/home/chbope/extension/script/deepsomatic/benchmark")
    ds_dir = benchmark_dir / "deepsomatic"
    clair_dir = benchmark_dir / "clairsto"
    output_dir = benchmark_dir / "comparison_results" / "per_sample_tables"

    # Create output directory
    output_dir.mkdir(exist_ok=True, parents=True)

    # Load summary data to get sample list
    summary_df = pd.read_csv(benchmark_dir / "comparison_results" / "per_sample_comparison.csv")

    for _, row in summary_df.iterrows():
        sample_id = row['sample_id']
        print(f"Processing {sample_id}...")

        ds_df, clair_df = load_sample_data(sample_id, ds_dir, clair_dir)

        # Extract variant info
        ds_variants = extract_variant_info(ds_df, 'DeepSomatic')
        clair_variants = extract_variant_info(clair_df, 'ClairS-TO')

        # Combine both
        combined = pd.concat([ds_variants, clair_variants], ignore_index=True)

        if len(combined) > 0:
            # Sort by chromosome and position
            combined = combined.sort_values(['Chr', 'Start'])

            # Save to CSV
            output_file = output_dir / f"{sample_id}_variants_detailed.csv"
            combined.to_csv(output_file, index=False)

            print(f"  ✓ Saved {len(combined)} variants to {output_file.name}")

    print(f"\n✓ All per-sample variant tables saved to: {output_dir}")
    print(f"  Total samples processed: {len(summary_df)}")

    return 0

if __name__ == "__main__":
    sys.exit(create_per_sample_tables())
