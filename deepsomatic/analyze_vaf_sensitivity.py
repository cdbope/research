#!/usr/bin/env python3
"""
VAF Sensitivity Analysis
Analyze which caller is better at detecting low VAF variants
"""

import pandas as pd
import sys
from pathlib import Path
import numpy as np

def analyze_vaf_distribution(df, caller_name, vaf_column):
    """Analyze VAF distribution and sensitivity"""

    if vaf_column not in df.columns:
        print(f"Warning: VAF column '{vaf_column}' not found for {caller_name}")
        return None

    vaf_values = df[vaf_column].dropna()

    if len(vaf_values) == 0:
        print(f"Warning: No VAF values found for {caller_name}")
        return None

    # VAF ranges
    very_low = vaf_values[vaf_values < 0.05]  # < 5%
    low = vaf_values[(vaf_values >= 0.05) & (vaf_values < 0.10)]  # 5-10%
    moderate = vaf_values[(vaf_values >= 0.10) & (vaf_values < 0.30)]  # 10-30%
    high = vaf_values[(vaf_values >= 0.30) & (vaf_values < 0.70)]  # 30-70%
    very_high = vaf_values[vaf_values >= 0.70]  # > 70%

    results = {
        'caller': caller_name,
        'total_variants': len(vaf_values),
        'very_low_vaf': len(very_low),  # < 5%
        'low_vaf': len(low),  # 5-10%
        'moderate_vaf': len(moderate),  # 10-30%
        'high_vaf': len(high),  # 30-70%
        'very_high_vaf': len(very_high),  # > 70%
        'min_vaf': vaf_values.min(),
        'max_vaf': vaf_values.max(),
        'mean_vaf': vaf_values.mean(),
        'median_vaf': vaf_values.median(),
        'vaf_std': vaf_values.std()
    }

    return results

def compare_low_vaf_detection(ds_df, clair_df):
    """Compare low VAF detection capabilities"""

    print("\n" + "="*80)
    print("LOW VAF VARIANT DETECTION ANALYSIS")
    print("="*80)

    # Extract VAF for DeepSomatic
    if 'Otherinfo10' in ds_df.columns:
        # Extract VAF from format: GT:GQ:DP:AD:VAF:PL
        ds_df['VAF'] = ds_df['Otherinfo10'].str.split(':').str[4].astype(float)
        ds_vaf_col = 'VAF'
    elif 'VAF_DS' in ds_df.columns:
        ds_vaf_col = 'VAF_DS'
    else:
        ds_vaf_col = None

    # Extract VAF for Clair
    if 'AF' in clair_df.columns:
        clair_vaf_col = 'AF'
    elif 'VAF_Clair' in clair_df.columns:
        clair_vaf_col = 'VAF_Clair'
    else:
        clair_vaf_col = None

    # Analyze distributions
    ds_results = analyze_vaf_distribution(ds_df, "DeepSomatic", ds_vaf_col) if ds_vaf_col else None
    clair_results = analyze_vaf_distribution(clair_df, "Clair3/ClairS-TO", clair_vaf_col) if clair_vaf_col else None

    if not ds_results or not clair_results:
        print("Error: Could not extract VAF information from one or both files")
        return None

    # Print comparison
    print("\n1. VAF DISTRIBUTION COMPARISON")
    print(f"\n   {'VAF Range':<20} {'DeepSomatic':<20} {'Clair3/ClairS-TO':<20} {'Better Detection'}")
    print(f"   {'-'*80}")

    ranges = [
        ('Very Low (< 5%)', 'very_low_vaf', 'Critical for early detection'),
        ('Low (5-10%)', 'low_vaf', 'Subclonal mutations'),
        ('Moderate (10-30%)', 'moderate_vaf', 'Heterozygous variants'),
        ('High (30-70%)', 'high_vaf', 'Dominant clones'),
        ('Very High (> 70%)', 'very_high_vaf', 'Homozygous/main clone')
    ]

    ds_low_total = ds_results['very_low_vaf'] + ds_results['low_vaf']
    clair_low_total = clair_results['very_low_vaf'] + clair_results['low_vaf']

    winner_counts = {'ds': 0, 'clair': 0, 'tie': 0}

    for range_name, key, description in ranges:
        ds_count = ds_results[key]
        clair_count = clair_results[key]

        if ds_count > clair_count:
            winner = "DeepSomatic ✓"
            winner_counts['ds'] += 1
        elif clair_count > ds_count:
            winner = "Clair3/ClairS-TO ✓"
            winner_counts['clair'] += 1
        else:
            winner = "Tie"
            winner_counts['tie'] += 1

        ds_pct = (ds_count / ds_results['total_variants'] * 100) if ds_results['total_variants'] > 0 else 0
        clair_pct = (clair_count / clair_results['total_variants'] * 100) if clair_results['total_variants'] > 0 else 0

        print(f"   {range_name:<20} {ds_count:>3} ({ds_pct:>5.1f}%)        {clair_count:>3} ({clair_pct:>5.1f}%)        {winner}")

    print(f"\n   {'TOTAL Low VAF':<20} {ds_low_total:>3}                {clair_low_total:>3}                ", end='')
    if ds_low_total > clair_low_total:
        print("DeepSomatic ✓")
        low_vaf_winner = "DeepSomatic"
    elif clair_low_total > ds_low_total:
        print("Clair3/ClairS-TO ✓")
        low_vaf_winner = "Clair3/ClairS-TO"
    else:
        print("Tie")
        low_vaf_winner = "Tie"

    # Statistical summary
    print(f"\n2. VAF STATISTICS")
    print(f"\n   {'Metric':<25} {'DeepSomatic':<20} {'Clair3/ClairS-TO':<20}")
    print(f"   {'-'*65}")
    print(f"   {'Minimum VAF':<25} {ds_results['min_vaf']:<20.4f} {clair_results['min_vaf']:<20.4f}")
    print(f"   {'Maximum VAF':<25} {ds_results['max_vaf']:<20.4f} {clair_results['max_vaf']:<20.4f}")
    print(f"   {'Mean VAF':<25} {ds_results['mean_vaf']:<20.4f} {clair_results['mean_vaf']:<20.4f}")
    print(f"   {'Median VAF':<25} {ds_results['median_vaf']:<20.4f} {clair_results['median_vaf']:<20.4f}")
    print(f"   {'Std Deviation':<25} {ds_results['vaf_std']:<20.4f} {clair_results['vaf_std']:<20.4f}")

    print(f"\n3. LOW VAF SENSITIVITY ANALYSIS")
    print(f"   {'-'*80}")

    # Minimum detectable VAF
    print(f"\n   Minimum Detectable VAF:")
    print(f"     DeepSomatic:      {ds_results['min_vaf']:.4f} ({ds_results['min_vaf']*100:.2f}%)")
    print(f"     Clair3/ClairS-TO: {clair_results['min_vaf']:.4f} ({clair_results['min_vaf']*100:.2f}%)")

    if ds_results['min_vaf'] < clair_results['min_vaf']:
        print(f"     → DeepSomatic detects lower VAF variants ✓")
        min_vaf_winner = "DeepSomatic"
    elif clair_results['min_vaf'] < ds_results['min_vaf']:
        print(f"     → Clair3/ClairS-TO detects lower VAF variants ✓")
        min_vaf_winner = "Clair3/ClairS-TO"
    else:
        print(f"     → Equal minimum VAF detection")
        min_vaf_winner = "Tie"

    # Low VAF variant count
    print(f"\n   Low VAF Variant Count (< 10%):")
    print(f"     DeepSomatic:      {ds_low_total} variants")
    print(f"     Clair3/ClairS-TO: {clair_low_total} variants")

    if ds_low_total > clair_low_total:
        print(f"     → DeepSomatic finds more low VAF variants ✓")
    elif clair_low_total > ds_low_total:
        print(f"     → Clair3/ClairS-TO finds more low VAF variants ✓")
    else:
        print(f"     → Equal low VAF variant detection")

    # Low VAF in cancer genes
    print(f"\n   Low VAF in Cancer Genes:")

    cancer_genes = {
        'TP53', 'KRAS', 'EGFR', 'PIK3CA', 'BRAF', 'PTEN', 'ALK', 'ROS1',
        'ERBB2', 'MET', 'FGFR1', 'FGFR2', 'FGFR3', 'PDGFRA', 'KIT', 'RET',
        'NRAS', 'HRAS', 'JAK2', 'IDH1', 'IDH2', 'NOTCH1', 'NOTCH2',
        'H3-3A', 'H3-3B', 'TERT', 'ATRX', 'DAXX', 'SETD2', 'ARID1A'
    }

    ds_low_cancer = 0
    clair_low_cancer = 0

    if 'Gene.refGene' in ds_df.columns and ds_vaf_col in ds_df.columns:
        ds_low_df = ds_df[ds_df[ds_vaf_col] < 0.10]
        for gene in cancer_genes:
            if ds_low_df['Gene.refGene'].astype(str).str.contains(gene, na=False).any():
                ds_low_cancer += 1

    if 'Gene.refGene' in clair_df.columns and clair_vaf_col in clair_df.columns:
        clair_low_df = clair_df[clair_df[clair_vaf_col] < 0.10]
        for gene in cancer_genes:
            if clair_low_df['Gene.refGene'].astype(str).str.contains(gene, na=False).any():
                clair_low_cancer += 1

    print(f"     DeepSomatic:      {ds_low_cancer} cancer genes with low VAF variants")
    print(f"     Clair3/ClairS-TO: {clair_low_cancer} cancer genes with low VAF variants")

    if ds_low_cancer > clair_low_cancer:
        print(f"     → DeepSomatic finds more low VAF variants in cancer genes ✓")
        cancer_low_winner = "DeepSomatic"
    elif clair_low_cancer > ds_low_cancer:
        print(f"     → Clair3/ClairS-TO finds more low VAF variants in cancer genes ✓")
        cancer_low_winner = "Clair3/ClairS-TO"
    else:
        print(f"     → Equal detection in cancer genes")
        cancer_low_winner = "Tie"

    print(f"\n4. LOW VAF VARIANT DETAILS")
    print(f"   {'-'*80}")

    # List low VAF variants from both callers
    if ds_vaf_col in ds_df.columns:
        ds_low = ds_df[ds_df[ds_vaf_col] < 0.10].copy()
        if len(ds_low) > 0:
            print(f"\n   DeepSomatic Low VAF Variants (< 10%):")
            print(f"   {'Variant':<30} {'Gene':<15} {'VAF':<10} {'Clinical':<20}")
            print(f"   {'-'*80}")
            for _, row in ds_low.iterrows():
                variant = f"{row['Chr']}:{row['Start']}:{row['Ref']}>{row['Alt']}"
                gene = str(row.get('Gene.refGene', 'N/A'))[:14]
                vaf = row[ds_vaf_col]
                clnsig = str(row.get('CLNSIG', 'N/A'))[:19]
                print(f"   {variant:<30} {gene:<15} {vaf:<10.4f} {clnsig:<20}")
        else:
            print(f"\n   DeepSomatic: No variants with VAF < 10%")

    if clair_vaf_col in clair_df.columns:
        clair_low = clair_df[clair_df[clair_vaf_col] < 0.10].copy()
        if len(clair_low) > 0:
            print(f"\n   Clair3/ClairS-TO Low VAF Variants (< 10%):")
            print(f"   {'Variant':<30} {'Gene':<15} {'VAF':<10} {'Clinical':<20}")
            print(f"   {'-'*80}")
            for _, row in clair_low.iterrows():
                variant = f"{row['Chr']}:{row['Start']}:{row['Ref']}>{row['Alt']}"
                gene = str(row.get('Gene.refGene', 'N/A'))[:14]
                vaf = row[clair_vaf_col]
                clnsig = str(row.get('CLNSIG', 'N/A'))[:19]
                print(f"   {variant:<30} {gene:<15} {vaf:<10.4f} {clnsig:<20}")
        else:
            print(f"\n   Clair3/ClairS-TO: No variants with VAF < 10%")

    print(f"\n5. FINAL RECOMMENDATION FOR LOW VAF DETECTION")
    print(f"   {'-'*80}")

    # Score each criterion
    scores = {
        'ds': 0,
        'clair': 0
    }

    if min_vaf_winner == "DeepSomatic":
        scores['ds'] += 2
    elif min_vaf_winner == "Clair3/ClairS-TO":
        scores['clair'] += 2

    if low_vaf_winner == "DeepSomatic":
        scores['ds'] += 2
    elif low_vaf_winner == "Clair3/ClairS-TO":
        scores['clair'] += 2

    if cancer_low_winner == "DeepSomatic":
        scores['ds'] += 3
    elif cancer_low_winner == "Clair3/ClairS-TO":
        scores['clair'] += 3

    print(f"\n   Scoring:")
    print(f"     - Minimum detectable VAF: {min_vaf_winner}")
    print(f"     - Total low VAF variants: {low_vaf_winner}")
    print(f"     - Low VAF in cancer genes: {cancer_low_winner}")
    print(f"")
    print(f"   Final Score:")
    print(f"     DeepSomatic:      {scores['ds']}/7 points")
    print(f"     Clair3/ClairS-TO: {scores['clair']}/7 points")
    print(f"")

    if scores['ds'] > scores['clair']:
        winner = "DeepSomatic"
        reason = "Better at detecting low VAF variants"
    elif scores['clair'] > scores['ds']:
        winner = "Clair3/ClairS-TO"
        reason = "Better at detecting low VAF variants"
    else:
        winner = "Tie"
        reason = "Similar performance for low VAF detection"

    print(f"   🏆 WINNER FOR LOW VAF DETECTION: {winner}")
    print(f"   Reason: {reason}")

    # Clinical implications
    print(f"\n6. CLINICAL IMPLICATIONS")
    print(f"   {'-'*80}")
    print(f"""
   Low VAF variants are critical for:
   ✓ Early cancer detection (< 5% VAF)
   ✓ Minimal residual disease monitoring (< 1% VAF)
   ✓ Subclonal population detection (5-15% VAF)
   ✓ Liquid biopsy applications (often < 1% VAF)

   For T25-152:
   - DeepSomatic detected down to: {ds_results['min_vaf']*100:.2f}% VAF
   - Clair3/ClairS-TO detected down to: {clair_results['min_vaf']*100:.2f}% VAF

   Recommendation:
   - For low VAF detection: Use {winner}
   - For comprehensive screening: Consider both callers
   - For clinical reporting: Validate all low VAF variants (< 10%)
   """)

    print(f"   {'-'*80}")

    return {
        'winner': winner,
        'ds_min_vaf': ds_results['min_vaf'],
        'clair_min_vaf': clair_results['min_vaf'],
        'ds_low_count': ds_low_total,
        'clair_low_count': clair_low_total,
        'ds_low_cancer': ds_low_cancer,
        'clair_low_cancer': clair_low_cancer
    }

def main():
    # Import comparison functions
    sys.path.insert(0, str(Path(__file__).parent))
    from compare_callers import load_deepsomatic, load_clair

    # File paths
    if len(sys.argv) == 4:
        deepsomatic_file = sys.argv[1]
        clair_file = sys.argv[2]
        output_dir = sys.argv[3]
    else:
        deepsomatic_file = "/home/chbope/extension/script/deepsomatic/output/T25-152_annotateandfilter_deep_somatic.csv"
        clair_file = "/home/chbope/extension/data/t25-152/merge_annot_clair3andclairsto/T25-152_merge_annotation_filter_snvs_allcall.csv"
        output_dir = "/home/chbope/extension/script/deepsomatic/comparison_results"

    # Check files
    if not Path(deepsomatic_file).exists():
        print(f"Error: DeepSomatic file not found: {deepsomatic_file}")
        sys.exit(1)

    if not Path(clair_file).exists():
        print(f"Error: Clair3/ClairS-TO file not found: {clair_file}")
        sys.exit(1)

    print("="*80)
    print("VAF SENSITIVITY ANALYSIS")
    print("="*80)
    print(f"\nAnalyzing: {Path(deepsomatic_file).name}")
    print(f"          vs {Path(clair_file).name}")

    # Load data
    ds_df = load_deepsomatic(deepsomatic_file)
    clair_df = load_clair(clair_file)

    # Analyze
    results = compare_low_vaf_detection(ds_df, clair_df)

    if results:
        # Save results
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        vaf_report = output_dir / "vaf_sensitivity_report.txt"
        with open(vaf_report, 'w') as f:
            f.write("LOW VAF SENSITIVITY REPORT\n")
            f.write("="*80 + "\n\n")
            f.write(f"Winner: {results['winner']}\n\n")
            f.write(f"DeepSomatic minimum VAF: {results['ds_min_vaf']*100:.2f}%\n")
            f.write(f"Clair3/ClairS-TO minimum VAF: {results['clair_min_vaf']*100:.2f}%\n\n")
            f.write(f"DeepSomatic low VAF variants (< 10%): {results['ds_low_count']}\n")
            f.write(f"Clair3/ClairS-TO low VAF variants (< 10%): {results['clair_low_count']}\n\n")
            f.write(f"DeepSomatic low VAF in cancer genes: {results['ds_low_cancer']}\n")
            f.write(f"Clair3/ClairS-TO low VAF in cancer genes: {results['clair_low_cancer']}\n")

        print(f"\nReport saved to: {vaf_report}")

    print("\n" + "="*80)

    return 0

if __name__ == "__main__":
    sys.exit(main())
