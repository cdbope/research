#!/usr/bin/env python3
"""
Complete pipeline to classify Nanopore methylation data using the Drexler neural glioblastoma classifier.

This script:
1. Converts bedMethyl data to beta values at Illumina 450K probe positions
2. Applies imputation for missing values
3. Runs the logistic regression classifier
4. Outputs the neural classification prediction

Usage:
    python run_nanopore_classification.py <bedmethyl_file> [--output_dir <dir>]

Example:
    python run_nanopore_classification.py /path/to/T001.wf_mods.bedmethyl.gz --output_dir results/

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

from nanopore_to_drexler import load_probe_coordinates, parse_bedmethyl, map_probes_to_nanopore


def load_classifier_data(base_path):
    """Load all required classifier data."""
    print("\nLoading classifier data...")

    # Load CpG probe lists
    low_cgs = pd.read_csv(
        os.path.join(base_path, "code/neural_group_classification/low_manifest_all.csv"),
        index_col=0
    ).index.tolist()

    high_cgs = pd.read_csv(
        os.path.join(base_path, "code/neural_group_classification/high_manifest_all.csv"),
        index_col=0
    ).index.tolist()

    cgs = low_cgs + high_cgs
    print(f"  Required CpG probes: {len(cgs)} (low: {len(low_cgs)}, high: {len(high_cgs)})")

    # Load imputer
    imputer = joblib.load(os.path.join(base_path, "code/neural_group_classification/imputer.pkl"))
    print(f"  Loaded imputer")

    # Load classifier
    clf = joblib.load(os.path.join(base_path, "code/neural_group_classification/logregCV_allCpG.pkl"))
    print(f"  Loaded classifier (features: {clf.coef_.shape[1]})")

    return cgs, low_cgs, high_cgs, imputer, clf


def impute_missing_values(df, imputer, cgs):
    """Impute missing values using the training set mean."""
    print("\nImputing missing values...")

    # Get fill values from imputer
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


def calculate_neural_signature_score(df, high_cgs):
    """
    Calculate the neural signature score as the mean beta value of high-neural CpG probes.

    This continuous score (typically 0.1-0.8) corresponds to the Y-axis in the paper's
    neural signature figures. Higher values indicate more neural-like methylation patterns.
    """
    # Get the high-neural probes that are present in the data
    available_high_cgs = [cg for cg in high_cgs if cg in df.columns]

    # Calculate mean beta value across high-neural probes for each sample
    neural_scores = df[available_high_cgs].mean(axis=1)

    return neural_scores


# Paper's neural signature threshold (Drexler et al., Nature Medicine 2024)
NEURAL_SIGNATURE_THRESHOLD = 0.41


def run_classification(df, clf, high_cgs):
    """Run the classifier and return predictions."""
    print("\nRunning classification...")

    # Calculate neural signature score (continuous, like in the paper's figures)
    neural_signature = calculate_neural_signature_score(df, high_cgs)

    # Classify based on paper's threshold of 0.41
    # Score >= 0.41 -> High-Neural, Score < 0.41 -> Low-Neural
    preds = (neural_signature >= NEURAL_SIGNATURE_THRESHOLD).astype(int)

    # Also compute logistic regression probabilities for reference
    preds_proba = clf.predict_proba(np.array(df))
    classes = clf.classes_

    results = pd.DataFrame({
        'Neural_Signature_Score': neural_signature,
        'Prediction': preds,
        'Neural_Classification': preds.map({0: 'Low-Neural', 1: 'High-Neural'}),
        'LR_Prob_Class_0': preds_proba[:, 0],
        'LR_Prob_Class_1': preds_proba[:, 1] if len(classes) > 1 else np.nan,
    }, index=df.index)

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Classify Nanopore methylation data using the Drexler neural glioblastoma classifier'
    )
    parser.add_argument('bedmethyl', help='Input bedMethyl file (from modkit)')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    parser.add_argument('--sample_name', help='Sample name (default: derived from filename)')
    parser.add_argument('--min_coverage', type=int, default=5,
                        help='Minimum read coverage per CpG (default: 5)')
    parser.add_argument('--genome', choices=['hg19', 'hg38'], default='hg38',
                        help='Reference genome build used for bedMethyl file (default: hg38)')
    parser.add_argument('--probe_coords', help='Pre-computed probe coordinates CSV')

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

    print("=" * 60)
    print("Nanopore Neural Glioblastoma Classification")
    print("=" * 60)
    print(f"Sample: {sample_name}")
    print(f"Input: {args.bedmethyl}")
    print(f"Output directory: {args.output_dir}")
    print(f"Minimum coverage: {args.min_coverage}x")
    print(f"Reference genome: {args.genome}")

    # Load classifier data
    cgs, low_cgs, high_cgs, imputer, clf = load_classifier_data(base_path)

    # Load probe coordinates
    print("\nLoading probe coordinates...")
    probe_coords = load_probe_coordinates(args.probe_coords, genome=args.genome)

    # Filter to required probes
    probe_coords_required = probe_coords[probe_coords.index.isin(cgs)]
    print(f"  Found coordinates for {len(probe_coords_required)}/{len(cgs)} required probes")

    # Parse bedMethyl file
    print("\n" + "=" * 60)
    print("Step 1: Parsing Nanopore methylation data")
    print("=" * 60)
    methylation_data = parse_bedmethyl(args.bedmethyl)

    # Map to probes
    print("\n" + "=" * 60)
    print("Step 2: Mapping to Illumina 450K probes")
    print("=" * 60)
    mapped_values = map_probes_to_nanopore(
        probe_coords_required, methylation_data,
        min_coverage=args.min_coverage
    )

    # Create DataFrame
    df = pd.DataFrame(mapped_values, index=[sample_name]).T

    # Ensure correct probe order
    df = df.reindex(cgs)
    df = df.T  # Transpose so samples are rows

    # Report coverage
    n_valid = df.loc[sample_name].notna().sum()
    coverage_pct = 100.0 * n_valid / len(cgs)
    print(f"\nProbe coverage: {n_valid}/{len(cgs)} ({coverage_pct:.1f}%)")

    # Save beta values
    betas_file = os.path.join(args.output_dir, f"{sample_name}_betas.csv")
    df.T.to_csv(betas_file)
    print(f"Saved beta values to: {betas_file}")

    # Impute missing values
    print("\n" + "=" * 60)
    print("Step 3: Imputing missing values")
    print("=" * 60)
    df = impute_missing_values(df, imputer, cgs)

    # Run classification
    print("\n" + "=" * 60)
    print("Step 4: Running classification")
    print("=" * 60)
    results = run_classification(df, clf, high_cgs)

    # Save results
    results_file = os.path.join(args.output_dir, f"{sample_name}_prediction.csv")
    results.to_csv(results_file)
    print(f"\nSaved predictions to: {results_file}")

    # Print results
    print("\n" + "=" * 60)
    print("CLASSIFICATION RESULTS")
    print("=" * 60)
    print(f"  (Using paper's threshold: {NEURAL_SIGNATURE_THRESHOLD})")
    for idx, row in results.iterrows():
        print(f"\nSample: {idx}")
        print(f"  Neural Signature Score: {row['Neural_Signature_Score']:.3f}")
        print(f"  Classification: {row['Neural_Classification']}")
        print(f"  Logistic Regression Probabilities (for reference):")
        print(f"    Low-Neural: {row['LR_Prob_Class_0']:.3f}")
        print(f"    High-Neural: {row['LR_Prob_Class_1']:.3f}")

    print("\n" + "=" * 60)
    print("Pipeline completed successfully!")
    print("=" * 60)

    return results


if __name__ == '__main__':
    main()
