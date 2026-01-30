#!/usr/bin/env python3
"""
Version 4: Nanopore Neural Glioblastoma Classification Pipeline with Normalization

Key changes from v3:
- Adds quantile normalization to match Illumina training distribution
- Uses imputer statistics as the reference distribution for normalization
- This addresses the systematic beta value differences between Nanopore and Illumina

The normalization step transforms Nanopore beta values to match the distribution
of the Illumina array data used to train the classifier.

Usage:
    python run_nanopore_classification_v4.py <bedmethyl_file> [--output_dir <dir>]

Example:
    python run_nanopore_classification_v4.py /path/to/T001.wf_mods.bedmethyl.gz --output_dir results/

Author: Generated for Drexler et al. pipeline adaptation
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import joblib

# Add parent directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from nanopore_to_drexler_v2 import load_epic_probe_coordinates, parse_bedmethyl, map_probes_to_nanopore


def load_classifier_data(base_path):
    """Load all required classifier data."""
    print("\nLoading classifier data...")

    # Load CpG probe lists (same order as original notebook)
    low_cgs = pd.read_csv(
        os.path.join(base_path, "code/neural_group_classification/low_manifest_all.csv"),
        index_col=0
    ).index.tolist()

    high_cgs = pd.read_csv(
        os.path.join(base_path, "code/neural_group_classification/high_manifest_all.csv"),
        index_col=0
    ).index.tolist()

    # Combined list in the correct order (low + high, as in original notebook)
    cgs = low_cgs + high_cgs
    print(f"  Required CpG probes: {len(cgs)} (low: {len(low_cgs)}, high: {len(high_cgs)})")

    # Load imputer (contains training data statistics)
    imputer = joblib.load(os.path.join(base_path, "code/neural_group_classification/imputer.pkl"))
    print(f"  Loaded imputer (training data reference)")

    # Load classifier (logistic regression trained on all 1289 CpGs)
    clf = joblib.load(os.path.join(base_path, "code/neural_group_classification/logregCV_allCpG.pkl"))
    print(f"  Loaded classifier (LogisticRegressionCV, features: {clf.coef_.shape[1]})")

    return cgs, low_cgs, high_cgs, imputer, clf


def get_training_reference(imputer, cgs):
    """
    Extract training data distribution from imputer for normalization reference.

    The imputer contains the mean values for each probe from the Illumina training data.
    We use these as the reference distribution for quantile normalization.
    """
    training_means = pd.Series(imputer.statistics_, index=imputer.feature_names_in_)
    return training_means.reindex(cgs)


def quantile_normalize(sample_values, reference_values):
    """
    Quantile normalize sample values to match the distribution of reference values.

    This transforms the Nanopore beta values to have the same distribution
    as the Illumina training data.

    Args:
        sample_values: pd.Series of Nanopore beta values
        reference_values: pd.Series of Illumina training means

    Returns:
        pd.Series of normalized beta values
    """
    # Get valid (non-NA) values for both
    valid_mask = sample_values.notna() & reference_values.notna()

    if valid_mask.sum() == 0:
        return sample_values

    sample_valid = sample_values[valid_mask].copy()
    ref_valid = reference_values[valid_mask].copy()

    n = len(sample_valid)

    # Rank the sample values (1 to n)
    sample_ranks = sample_valid.rank(method='average')

    # Sort reference values to get the target distribution
    ref_sorted = np.sort(ref_valid.values)

    # Map ranks to reference distribution
    # Convert ranks to 0-based indices
    indices = ((sample_ranks - 1) / (n - 1) * (n - 1)).round().astype(int).clip(0, n-1)

    # Create normalized values
    normalized = pd.Series(ref_sorted[indices.values], index=sample_valid.index)

    # Build result with NAs preserved
    result = sample_values.copy()
    result[valid_mask] = normalized

    return result


def zscore_normalize(sample_values, reference_values):
    """
    Z-score normalize sample values to match reference distribution.

    Transforms sample values to have the same mean and std as reference.

    Args:
        sample_values: pd.Series of Nanopore beta values
        reference_values: pd.Series of Illumina training means

    Returns:
        pd.Series of normalized beta values
    """
    valid_mask = sample_values.notna() & reference_values.notna()

    sample_valid = sample_values[valid_mask]
    ref_valid = reference_values[valid_mask]

    # Calculate statistics
    sample_mean = sample_valid.mean()
    sample_std = sample_valid.std()
    ref_mean = ref_valid.mean()
    ref_std = ref_valid.std()

    # Z-score transform then rescale
    normalized = (sample_values - sample_mean) / sample_std * ref_std + ref_mean

    # Clip to valid beta range [0, 1]
    normalized = normalized.clip(0, 1)

    return normalized


def impute_missing_values(df, imputer, cgs):
    """
    Impute missing values using the training set mean.
    """
    print("\nImputing missing values...")

    # Get fill values from imputer (training set means)
    df_fill_values = pd.DataFrame(
        imputer.statistics_,
        columns=["fill_value"],
        index=imputer.feature_names_in_
    )

    # Find columns with missing values
    missing_cols = df.columns[df.isna().any()].tolist()
    print(f"  Columns with missing values: {len(missing_cols)}")

    # Impute
    for idx in df.index:
        for col in missing_cols:
            if pd.isnull(df.loc[idx, col]):
                if col in df_fill_values.index:
                    df.loc[idx, col] = df_fill_values.loc[col, "fill_value"]

    # Check remaining missing
    still_missing = df.isna().sum().sum()
    if still_missing > 0:
        print(f"  WARNING: {still_missing} values still missing after imputation")

    return df


def run_classification(df, clf):
    """
    Run the logistic regression classifier.
    """
    print("\nRunning logistic regression classification...")

    preds = clf.predict(np.array(df))
    preds_proba = clf.predict_proba(np.array(df))
    preds_score = preds_proba.max(axis=1)

    results = pd.DataFrame({
        'Prediction': preds.astype(int),
        'Prediction_score': preds_score,
        'Neural_Classification': pd.Series(preds).map({0: 'Low-Neural', 1: 'High-Neural'}).values
    }, index=df.index)

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Classify Nanopore methylation data using the Drexler classifier (v4 - with normalization)'
    )
    parser.add_argument('bedmethyl', help='Input bedMethyl file (from modkit)')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    parser.add_argument('--sample_name', help='Sample name (default: derived from filename)')
    parser.add_argument('--min_coverage', type=int, default=5,
                        help='Minimum read coverage per CpG (default: 5)')
    parser.add_argument('--epic_bed',
                        default='/home/chbope/extension/nWGS_manuscript_data/data/reference/EPIC_sites_NEW.bed',
                        help='EPIC probe coordinates BED file (hg38)')
    parser.add_argument('--normalization', choices=['none', 'quantile', 'zscore'], default='quantile',
                        help='Normalization method (default: quantile)')

    args = parser.parse_args()

    # Setup paths
    base_path = os.path.dirname(script_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    # Determine sample name
    if args.sample_name:
        sample_name = args.sample_name
    else:
        sample_name = os.path.basename(args.bedmethyl)
        for suffix in ['.wf_mods.bedmethyl.gz', '.bedmethyl.gz', '.bedmethyl', '.gz']:
            sample_name = sample_name.replace(suffix, '')

    print("=" * 70)
    print("Nanopore Neural Glioblastoma Classification (v4)")
    print("With Normalization to Illumina Training Distribution")
    print("=" * 70)
    print(f"Sample: {sample_name}")
    print(f"Input: {args.bedmethyl}")
    print(f"Output directory: {args.output_dir}")
    print(f"Minimum coverage: {args.min_coverage}x")
    print(f"Normalization: {args.normalization}")

    # Load classifier data
    cgs, low_cgs, high_cgs, imputer, clf = load_classifier_data(base_path)

    # Get training reference distribution
    training_ref = get_training_reference(imputer, cgs)

    # Load EPIC probe coordinates
    print("\n" + "=" * 70)
    print("Step 1: Loading EPIC probe coordinates")
    print("=" * 70)
    probe_coords = load_epic_probe_coordinates(args.epic_bed)
    probe_coords_required = probe_coords[probe_coords.index.isin(cgs)]
    print(f"  Found coordinates for {len(probe_coords_required)}/{len(cgs)} required probes")

    # Parse bedMethyl file
    print("\n" + "=" * 70)
    print("Step 2: Parsing Nanopore methylation data")
    print("=" * 70)
    methylation_data = parse_bedmethyl(args.bedmethyl)

    # Map to probes
    print("\n" + "=" * 70)
    print("Step 3: Mapping to EPIC probes")
    print("=" * 70)
    mapped_values = map_probes_to_nanopore(
        probe_coords_required, methylation_data,
        min_coverage=args.min_coverage
    )

    # Create DataFrame
    df = pd.DataFrame(mapped_values, index=[sample_name]).T
    df = df.reindex(cgs)

    # Get raw values as Series for normalization
    raw_values = df[sample_name]

    # Report raw statistics
    n_valid = raw_values.notna().sum()
    coverage_pct = 100.0 * n_valid / len(cgs)
    print(f"\nProbe coverage: {n_valid}/{len(cgs)} ({coverage_pct:.1f}%)")

    print(f"\nRaw Nanopore beta statistics:")
    print(f"  Mean: {raw_values.mean():.4f}")
    print(f"  Std:  {raw_values.std():.4f}")
    print(f"  Min:  {raw_values.min():.4f}")
    print(f"  Max:  {raw_values.max():.4f}")

    print(f"\nIllumina training reference statistics:")
    print(f"  Mean: {training_ref.mean():.4f}")
    print(f"  Std:  {training_ref.std():.4f}")

    # Apply normalization
    print("\n" + "=" * 70)
    print(f"Step 4: Applying {args.normalization} normalization")
    print("=" * 70)

    if args.normalization == 'quantile':
        normalized_values = quantile_normalize(raw_values, training_ref)
    elif args.normalization == 'zscore':
        normalized_values = zscore_normalize(raw_values, training_ref)
    else:
        normalized_values = raw_values
        print("  Skipping normalization (using raw values)")

    if args.normalization != 'none':
        print(f"\nNormalized beta statistics:")
        print(f"  Mean: {normalized_values.mean():.4f}")
        print(f"  Std:  {normalized_values.std():.4f}")
        print(f"  Min:  {normalized_values.min():.4f}")
        print(f"  Max:  {normalized_values.max():.4f}")

    # Update DataFrame with normalized values
    df[sample_name] = normalized_values
    df = df.T  # Transpose so samples are rows

    # Save beta values (normalized)
    betas_file = os.path.join(args.output_dir, f"{sample_name}_betas_normalized.csv")
    df.T.to_csv(betas_file)
    print(f"\nSaved normalized beta values to: {betas_file}")

    # Impute missing values
    print("\n" + "=" * 70)
    print("Step 5: Imputing missing values")
    print("=" * 70)
    df = impute_missing_values(df, imputer, cgs)

    # Run classification
    print("\n" + "=" * 70)
    print("Step 6: Running Logistic Regression classification")
    print("=" * 70)
    results = run_classification(df, clf)

    # Save results
    results_file = os.path.join(args.output_dir, f"{sample_name}_prediction_normalized.csv")
    results.to_csv(results_file)
    print(f"\nSaved predictions to: {results_file}")

    # Print results
    print("\n" + "=" * 70)
    print("CLASSIFICATION RESULTS")
    print("=" * 70)
    for idx, row in results.iterrows():
        print(f"\nSample: {idx}")
        print(f"  Prediction: {int(row['Prediction'])} ({row['Neural_Classification']})")
        print(f"  Prediction score: {row['Prediction_score']:.4f}")

    print("\n" + "=" * 70)
    print("Pipeline completed successfully!")
    print("=" * 70)

    return results


if __name__ == '__main__':
    main()
