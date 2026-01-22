#!/usr/bin/env python3
"""
Create CSV file with detailed per-sample variant comparison
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

def get_variant_description(row):
    """Create human-readable variant description"""
    gene = row.get('Gene.refGene', 'Unknown')
    aa_change = row.get('AAChange.refGene', '')

    # Extract the actual amino acid change if available
    if pd.notna(aa_change) and aa_change != '':
        # Format: GENE:NM_xxx:exonX:c.xxx:p.xxx
        # We want the p.xxx part
        parts = str(aa_change).split(':')
        if len(parts) >= 5:
            aa = parts[4]
        else:
            aa = ''
    else:
        aa = ''

    # Get VAF
    vaf = row.get('VAF', 0)
    if isinstance(vaf, str):
        try:
            vaf = float(vaf)
        except:
            vaf = 0

    pos = f"{row['Chr']}:{row['Start']}"
    ref_alt = f"{row['Ref']}>{row['Alt']}"

    if aa and aa != '.':
        return f"{gene}:{aa} ({pos} {ref_alt} VAF={vaf:.1%})"
    else:
        return f"{gene} ({pos} {ref_alt} VAF={vaf:.1%})"

def create_variant_comparison_csv():
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

        if ds_df is None and clair_df is None:
            continue

        # Get variant keys
        ds_variants = {}
        clair_variants = {}

        if ds_df is not None and len(ds_df) > 0:
            for idx, variant_row in ds_df.iterrows():
                key = get_variant_key(variant_row)
                desc = get_variant_description(variant_row)
                ds_variants[key] = desc

        if clair_df is not None and len(clair_df) > 0:
            for idx, variant_row in clair_df.iterrows():
                key = get_variant_key(variant_row)
                desc = get_variant_description(variant_row)
                clair_variants[key] = desc

        # Find common and unique variants
        ds_keys = set(ds_variants.keys())
        clair_keys = set(clair_variants.keys())

        common_keys = ds_keys & clair_keys
        ds_only_keys = ds_keys - clair_keys
        clair_only_keys = clair_keys - ds_keys

        # Format variant lists
        deepsomatic_variants = "; ".join([ds_variants[k] for k in sorted(ds_variants.keys())])
        clairsto_variants = "; ".join([clair_variants[k] for k in sorted(clair_variants.keys())])
        common_variants = "; ".join([ds_variants[k] for k in sorted(common_keys)])
        missing_in_clairsto = "; ".join([ds_variants[k] for k in sorted(ds_only_keys)])
        missing_in_deepsomatic = "; ".join([clair_variants[k] for k in sorted(clair_only_keys)])

        results.append({
            'sample_id': sample_id,
            'deepsomatic_count': len(ds_variants),
            'clairsto_count': len(clair_variants),
            'common_count': len(common_keys),
            'missing_in_clairsto_count': len(ds_only_keys),
            'missing_in_deepsomatic_count': len(clair_only_keys),
            'deepsomatic_variants': deepsomatic_variants if deepsomatic_variants else 'None',
            'clairsto_variants': clairsto_variants if clairsto_variants else 'None',
            'common_variants': common_variants if common_variants else 'None',
            'missing_in_clairsto': missing_in_clairsto if missing_in_clairsto else 'None',
            'missing_in_deepsomatic': missing_in_deepsomatic if missing_in_deepsomatic else 'None'
        })

    # Create DataFrame
    result_df = pd.DataFrame(results)

    # Reorder columns
    columns_order = [
        'sample_id',
        'deepsomatic_count',
        'clairsto_count',
        'common_count',
        'missing_in_clairsto_count',
        'missing_in_deepsomatic_count',
        'deepsomatic_variants',
        'clairsto_variants',
        'common_variants',
        'missing_in_clairsto',
        'missing_in_deepsomatic'
    ]

    result_df = result_df[columns_order]

    # Save to CSV
    output_file = output_dir / "variant_comparison_detailed.csv"
    result_df.to_csv(output_file, index=False)

    print(f"\n✓ Variant comparison CSV saved to: {output_file}")
    print(f"  Total samples: {len(result_df)}")
    print(f"  Total DeepSomatic variants: {result_df['deepsomatic_count'].sum()}")
    print(f"  Total ClairS-TO variants: {result_df['clairsto_count'].sum()}")
    print(f"  Total common variants: {result_df['common_count'].sum()}")

    return 0

if __name__ == "__main__":
    sys.exit(create_variant_comparison_csv())
