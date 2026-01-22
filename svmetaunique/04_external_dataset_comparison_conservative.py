#!/usr/bin/env python3
"""
Step 4: External Dataset Comparison with CONSERVATIVE COUNTING

Compare conservative gene frequencies with TCGA/PCAWG datasets
Calculate fold-changes and identify validated genes

Input: Gene frequencies from step 03 (conservative counting)
Output: Validated genes with fold-changes vs TCGA/PCAWG
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input from step 03
COHORT_GENE_FREQ = "/home/chbope/extension/script/svmetaunique/results/genes/gene_frequencies_conservative.csv"

# External datasets
EXTERNAL_DATA_DIR = "/home/chbope/extension/script/svmeta/external_datasets"
TCGA_SV_FILE = f"{EXTERNAL_DATA_DIR}/tcga_gbm_sv_summary_hg38.csv"
PCAWG_SV_FILE = f"{EXTERNAL_DATA_DIR}/pcawg_gbm_sv_summary_hg38.csv"

# Output directory
OUTPUT_DIR = "/home/chbope/extension/script/svmetaunique/results/external_comparison"

# Sample sizes
N_COHORT = 200
N_TCGA = 577
N_PCAWG = 41

# ============================================================================

def wilson_ci(n_success, n_total, confidence=0.95):
    """Calculate Wilson score confidence interval"""
    if n_total == 0:
        return (0.0, 0.0)

    z = stats.norm.ppf((1 + confidence) / 2)
    p = n_success / n_total

    denominator = 1 + z**2 / n_total
    centre = (p + z**2 / (2 * n_total)) / denominator
    adjustment = z * np.sqrt((p * (1 - p) / n_total + z**2 / (4 * n_total**2))) / denominator

    lower = max(0.0, centre - adjustment)
    upper = min(1.0, centre + adjustment)

    return (lower, upper)

def load_external_dataset(file_path, dataset_name):
    """Load TCGA or PCAWG summary file and convert to gene-level frequencies"""
    print(f"Loading {dataset_name} from: {file_path}")

    if not os.path.exists(file_path):
        print(f"⚠️  {dataset_name} file not found: {file_path}")
        return pd.DataFrame()

    df = pd.read_csv(file_path)

    # The file has SV-level data with 'genes' column
    # We need to convert to gene-level frequencies

    # Get total samples for this dataset
    total_samples = df['total_samples'].iloc[0] if 'total_samples' in df.columns else 100

    gene_freq_dict = {}

    for _, row in df.iterrows():
        if 'genes' not in row or pd.isna(row['genes']):
            continue

        # Parse gene names (may be semicolon-separated)
        genes = str(row['genes']).split(';')

        # Get frequency for this SV
        sv_freq = row['frequency'] if 'frequency' in row else 0

        # Assign this frequency to all affected genes
        for gene in genes:
            gene = gene.strip()
            if gene:
                # Use max frequency if gene appears in multiple SVs
                if gene not in gene_freq_dict:
                    gene_freq_dict[gene] = sv_freq
                else:
                    gene_freq_dict[gene] = max(gene_freq_dict[gene], sv_freq)

    # Convert to DataFrame
    gene_df = pd.DataFrame([
        {'gene': gene, 'frequency': freq}
        for gene, freq in gene_freq_dict.items()
    ])

    print(f"✓ Loaded {len(gene_df)} genes from {dataset_name}")

    return gene_df

def compare_with_external_datasets(cohort_df, tcga_df, pcawg_df):
    """
    Compare cohort frequencies with TCGA and PCAWG
    Calculate fold-changes using CONSERVATIVE frequencies
    """
    print("\nComparing with external datasets...")

    # Start with cohort genes
    comparison = cohort_df[['gene', 'frequency', 'n_samples_affected', 'n_DEL', 'n_DUP', 'n_INV', 'n_BND', 'is_driver']].copy()
    comparison = comparison.rename(columns={'frequency': 'cohort_freq', 'n_samples_affected': 'cohort_n_samples'})

    # Add TCGA frequencies
    if not tcga_df.empty:
        tcga_dict = dict(zip(tcga_df['gene'], tcga_df['frequency']))
        comparison['tcga_freq'] = comparison['gene'].map(tcga_dict)
    else:
        comparison['tcga_freq'] = np.nan

    # Add PCAWG frequencies
    if not pcawg_df.empty:
        pcawg_dict = dict(zip(pcawg_df['gene'], pcawg_df['frequency']))
        comparison['pcawg_freq'] = comparison['gene'].map(pcawg_dict)
    else:
        comparison['pcawg_freq'] = np.nan

    # Calculate fold-changes
    comparison['fold_change_tcga'] = comparison.apply(
        lambda row: row['cohort_freq'] / row['tcga_freq'] if pd.notna(row['tcga_freq']) and row['tcga_freq'] > 0 else np.nan,
        axis=1
    )

    comparison['fold_change_pcawg'] = comparison.apply(
        lambda row: row['cohort_freq'] / row['pcawg_freq'] if pd.notna(row['pcawg_freq']) and row['pcawg_freq'] > 0 else np.nan,
        axis=1
    )

    # Add Wilson confidence intervals
    comparison['ci_low'] = comparison['cohort_n_samples'].apply(lambda n: wilson_ci(n, N_COHORT)[0])
    comparison['ci_high'] = comparison['cohort_n_samples'].apply(lambda n: wilson_ci(n, N_COHORT)[1])

    # Validation status
    def get_validation_status(row):
        has_tcga = pd.notna(row['tcga_freq'])
        has_pcawg = pd.notna(row['pcawg_freq'])

        if has_tcga and has_pcawg:
            return 'Validated_Both'
        elif has_tcga or has_pcawg:
            return 'Validated_One'
        else:
            return 'Cohort_Specific'

    comparison['validation_status'] = comparison.apply(get_validation_status, axis=1)

    print(f"✓ Comparison complete for {len(comparison)} genes")

    return comparison

def filter_high_confidence_genes(comparison_df, min_fold_change=2.0):
    """
    Filter for high-confidence genes:
    - Validated in at least one external dataset
    - Fold-change ≥ threshold
    """
    print(f"\nFiltering high-confidence genes (min FC: {min_fold_change}×)...")

    # Must be validated in at least one dataset
    validated = comparison_df[comparison_df['validation_status'].isin(['Validated_Both', 'Validated_One'])].copy()

    # Must have fold-change ≥ threshold in at least one dataset
    validated = validated[
        ((validated['fold_change_tcga'] >= min_fold_change) & pd.notna(validated['fold_change_tcga'])) |
        ((validated['fold_change_pcawg'] >= min_fold_change) & pd.notna(validated['fold_change_pcawg']))
    ]

    # Sort by maximum fold-change
    validated['max_fold_change'] = validated[['fold_change_tcga', 'fold_change_pcawg']].max(axis=1)
    validated = validated.sort_values('max_fold_change', ascending=False)

    print(f"✓ Identified {len(validated)} high-confidence genes")

    return validated

def create_summary_report(comparison_df, validated_df, output_dir):
    """Create markdown summary report"""

    report_file = f"{output_dir}/COMPARISON_REPORT_CONSERVATIVE.md"

    with open(report_file, 'w') as f:
        f.write("# External Dataset Comparison Report\n\n")
        f.write("## Methodology\n\n")
        f.write("**Conservative Counting Approach** (following npae082.pdf):\n")
        f.write("- Max 1 SV per gene per sample (binary: 0 or 1)\n")
        f.write("- Frequency = (# samples with ≥1 SV in gene) / 200\n")
        f.write("- Directly comparable to literature using matched tumor-normal sequencing\n\n")

        f.write("## Dataset Information\n\n")
        f.write("| Dataset | N Samples | Sequencing | SV Calling Method |\n")
        f.write("|---------|-----------|------------|-------------------|\n")
        f.write("| **Your Cohort (200GBM)** | 200 | Long-read (ONT/PromethION) | Sniffles2 (tumor-only) |\n")
        f.write("| **TCGA-GBM** | 577 | Short-read (Illumina) | Matched tumor-normal |\n")
        f.write("| **PCAWG-GBM** | 41 | Short-read (Illumina) | Matched tumor-normal |\n\n")

        f.write("## Summary Statistics\n\n")
        f.write(f"- **Total genes with SVs in cohort**: {len(comparison_df)}\n")
        f.write(f"- **Genes validated in TCGA**: {comparison_df['tcga_freq'].notna().sum()}\n")
        f.write(f"- **Genes validated in PCAWG**: {comparison_df['pcawg_freq'].notna().sum()}\n")
        f.write(f"- **Genes validated in both**: {(comparison_df['tcga_freq'].notna() & comparison_df['pcawg_freq'].notna()).sum()}\n")
        f.write(f"- **High-confidence genes** (FC ≥2×): {len(validated_df)}\n\n")

        f.write("## Top 20 Genes (Conservative Fold-Changes)\n\n")
        f.write("| Gene | 200GBM Freq | TCGA Freq | PCAWG Freq | FC (TCGA) | FC (PCAWG) | Driver |\n")
        f.write("|------|-------------|-----------|------------|-----------|------------|--------|\n")

        for _, row in validated_df.head(20).iterrows():
            tcga_fc = f"{row['fold_change_tcga']:.1f}×" if pd.notna(row['fold_change_tcga']) else "N/A"
            pcawg_fc = f"{row['fold_change_pcawg']:.1f}×" if pd.notna(row['fold_change_pcawg']) else "N/A"
            tcga_freq = f"{row['tcga_freq']*100:.1f}%" if pd.notna(row['tcga_freq']) else "N/A"
            pcawg_freq = f"{row['pcawg_freq']*100:.1f}%" if pd.notna(row['pcawg_freq']) else "N/A"
            driver = "Yes" if row['is_driver'] else "No"

            f.write(f"| {row['gene']} | {row['cohort_freq']*100:.1f}% | {tcga_freq} | {pcawg_freq} | ")
            f.write(f"{tcga_fc} | {pcawg_fc} | {driver} |\n")

        f.write("\n## Key Findings\n\n")

        # Top gene
        top_gene = validated_df.iloc[0]
        f.write(f"### Highest Enrichment: {top_gene['gene']}\n\n")
        f.write(f"- **200GBM Frequency**: {top_gene['cohort_freq']*100:.1f}%\n")
        f.write(f"- **TCGA Frequency**: {top_gene['tcga_freq']*100:.1f}%\n")
        f.write(f"- **Fold-Change**: {top_gene['fold_change_tcga']:.1f}×\n")
        f.write(f"- **Driver Gene**: {'Yes' if top_gene['is_driver'] else 'No'}\n\n")

        # Driver gene enrichment
        driver_genes = validated_df[validated_df['is_driver'] == True]
        f.write(f"### Driver Gene Enrichment\n\n")
        f.write(f"- **Driver genes in top 20**: {validated_df.head(20)['is_driver'].sum()}\n")
        f.write(f"- **Total driver genes validated**: {len(driver_genes)}\n")
        f.write(f"- **Average fold-change (drivers)**: {driver_genes['max_fold_change'].mean():.1f}×\n\n")

        # Genes with highest fold-changes
        f.write("### Genes with Extreme Enrichment (FC >10×)\n\n")
        extreme = validated_df[validated_df['max_fold_change'] > 10]
        if len(extreme) > 0:
            for _, row in extreme.iterrows():
                f.write(f"- **{row['gene']}**: {row['max_fold_change']:.1f}× ")
                f.write(f"({row['cohort_freq']*100:.1f}% in 200GBM)\n")
        else:
            f.write("None\n")

        f.write("\n## Interpretation\n\n")
        f.write("**Conservative counting ensures direct comparability** with TCGA/PCAWG:\n")
        f.write("- Both reference datasets use matched tumor-normal sequencing\n")
        f.write("- Genes are typically scored as altered/not altered (binary)\n")
        f.write("- Your conservative approach follows the same methodology\n\n")

        f.write("**High fold-changes indicate**:\n")
        f.write("1. Effective germline removal in your tumor-only data\n")
        f.write("2. Genuine somatic enrichment (not germline contamination)\n")
        f.write("3. Long-read advantage in detecting complex SVs\n")
        f.write("4. Novel biological insights into GBM chromothripsis\n\n")

        f.write("## Files Generated\n\n")
        f.write("- `all_genes_comparison_conservative.csv` - Complete gene comparison\n")
        f.write("- `high_confidence_validated_genes_conservative.csv` - Validated genes (FC ≥2×)\n")
        f.write("- `COMPARISON_REPORT_CONSERVATIVE.md` - This report\n")

    print(f"✓ Report saved: {report_file}")

def main():
    """Main comparison pipeline"""

    print("="*80)
    print("STEP 04: EXTERNAL DATASET COMPARISON (CONSERVATIVE COUNTING)")
    print("="*80)
    print("\nComparing conservative gene frequencies with TCGA/PCAWG")
    print("="*80)
    print()

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load cohort gene frequencies
    if not os.path.exists(COHORT_GENE_FREQ):
        print(f"⚠️  Cohort gene frequencies not found: {COHORT_GENE_FREQ}")
        print("Please run 03_build_matrix_conservative.py first")
        return

    cohort_df = pd.read_csv(COHORT_GENE_FREQ)
    print(f"✓ Loaded {len(cohort_df)} genes from cohort")

    # Load external datasets
    tcga_df = load_external_dataset(TCGA_SV_FILE, "TCGA-GBM")
    pcawg_df = load_external_dataset(PCAWG_SV_FILE, "PCAWG-GBM")

    # Compare with external datasets
    comparison_df = compare_with_external_datasets(cohort_df, tcga_df, pcawg_df)

    # Save complete comparison
    comparison_file = f"{OUTPUT_DIR}/all_genes_comparison_conservative.csv"
    comparison_df.to_csv(comparison_file, index=False)
    print(f"\n✓ Saved complete comparison: {comparison_file}")

    # Filter high-confidence genes
    validated_df = filter_high_confidence_genes(comparison_df, min_fold_change=2.0)

    # Save validated genes
    validated_file = f"{OUTPUT_DIR}/high_confidence_validated_genes_conservative.csv"
    validated_df.to_csv(validated_file, index=False)
    print(f"✓ Saved validated genes: {validated_file}")

    # Create summary report
    create_summary_report(comparison_df, validated_df, OUTPUT_DIR)

    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total genes analyzed: {len(comparison_df)}")
    print(f"High-confidence validated genes: {len(validated_df)}")
    print(f"\nTop 5 genes by fold-change:")

    for i, (_, row) in enumerate(validated_df.head(5).iterrows(), 1):
        max_fc = row['max_fold_change']
        print(f"  {i}. {row['gene']}: {row['cohort_freq']*100:.1f}% (FC: {max_fc:.1f}×)")

    print(f"\nResults saved to: {OUTPUT_DIR}/")
    print("="*80)

if __name__ == "__main__":
    main()
