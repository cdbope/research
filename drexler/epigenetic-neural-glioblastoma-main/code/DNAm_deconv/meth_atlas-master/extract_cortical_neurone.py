#!/usr/bin/env python3
"""
Extract Cortical_neurons row from deconvolution output and transpose to have
samples as rows and Cortical_neurons values as a column.

Usage:
    python extract_cortical_neurone.py <input_file> [--output <output_file>]

Example:
    python extract_cortical_neurone.py combined_deconv_betas_deconv_output.csv
"""

import pandas as pd
import argparse
import os


def extract_cortical_neurons(input_file, output_file=None):
    """Extract Cortical_neurons row and transpose."""

    # Read the deconvolution output
    df = pd.read_csv(input_file, index_col=0)

    # Extract the Cortical_neurons row
    if 'Cortical_neurons' not in df.index:
        raise ValueError("Cortical_neurons not found in the input file")

    cortical_values = df.loc['Cortical_neurons']

    # Create output DataFrame with samples as rows
    result = pd.DataFrame({
        'Cortical_neurons': cortical_values
    })
    result.index.name = 'Sample_id'

    # Set output file name
    if output_file is None:
        output_file = 'GBM_Moss_cortical_neurone_only_signatures.csv'

    # Save to CSV
    result.to_csv(output_file)

    print(f"Extracted Cortical_neurons for {len(result)} samples")
    print(f"Saved to: {output_file}")
    print(f"\nPreview:")
    print(result.head(10))

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Extract Cortical_neurons from deconvolution output'
    )
    parser.add_argument('input_file', help='Input deconvolution output CSV file')
    parser.add_argument('--output', '-o', default='GBM_Moss_signatures.csv',
                        help='Output file name (default: GBM_Moss_signatures.csv)')

    args = parser.parse_args()

    extract_cortical_neurons(args.input_file, args.output)


if __name__ == '__main__':
    main()
