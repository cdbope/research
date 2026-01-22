#!/usr/bin/env python3
"""
Generate beta values CSV file from Nanopore bedMethyl data for deconvolution.
Version 2: Supports batch processing of multiple samples.

This script:
1. Reads bedMethyl files (Nanopore methylation data) for multiple samples
2. Reads EPIC probe coordinates (probe ID -> chr:position)
3. Reads reference atlas to get required probe IDs
4. Maps Nanopore methylation to EPIC probes
5. Outputs beta values CSV files compatible with deconvolve.py

Usage:
    # Single sample
    python generate_beta_for_deconv_v2.py --bedmethyl /path/to/sample.wf_mods.bedmethyl.gz --output_dir results/

    # Batch processing
    python generate_beta_for_deconv_v2.py --input_dir /path/to/bedmethyl_files --sample_ids samples.txt --output_dir results/

Example:
    python generate_beta_for_deconv_v2.py \
        --input_dir /media/chbope/Expansion/200gbms_bedmethyl \
        --sample_ids /media/chbope/Expansion/200gbms_bedmethyl/samplesids.txt \
        --output_dir /home/chbope/extension/script/drexler/results/200gbm/deconv

Author: Generated for Drexler pipeline
"""

import argparse
import gzip
import os
import sys
import pandas as pd
from collections import defaultdict


def parse_bedmethyl(bedmethyl_path):
    """
    Parse bedMethyl file and return methylation data.

    Returns dict: {(chr, pos): (percent_methylated, coverage)}

    Note: bedMethyl files from modkit use '.' for strand and may have multiple
    entries per position (one for each strand/modification type). We aggregate
    by taking the entry with highest coverage.
    """
    print(f"  Parsing bedMethyl file: {os.path.basename(bedmethyl_path)}")

    methylation_data = {}

    open_func = gzip.open if bedmethyl_path.endswith('.gz') else open
    mode = 'rt' if bedmethyl_path.endswith('.gz') else 'r'

    with open_func(bedmethyl_path, mode) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            chrom = fields[0]
            start = int(fields[1])  # 0-based
            coverage = int(fields[9])
            percent_meth = float(fields[10])

            # Store by (chr, pos) only - ignore strand since bedMethyl uses '.'
            key = (chrom, start)

            # Keep the entry with highest coverage if multiple entries exist
            if key not in methylation_data or coverage > methylation_data[key][1]:
                methylation_data[key] = (percent_meth, coverage)

    print(f"    Total unique CpG positions: {len(methylation_data):,}")
    return methylation_data


def load_epic_coordinates(epic_bed_path):
    """
    Load EPIC probe coordinates from BED file.

    Returns dict: {probe_id: (chr, start, end, strand)}
    """
    print(f"Loading EPIC probe coordinates: {epic_bed_path}")

    probe_coords = {}

    with open(epic_bed_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            probe_id = fields[4]
            strand = fields[5]

            # Skip non-cg probes (like MGMT region)
            if probe_id.startswith('cg'):
                probe_coords[probe_id] = (chrom, start, end, strand)

    print(f"  Loaded {len(probe_coords):,} EPIC probe coordinates")
    return probe_coords


def load_atlas_probes(atlas_path):
    """
    Load reference atlas to get required probe IDs.

    Returns set of probe IDs
    """
    print(f"Loading reference atlas probes: {atlas_path}")

    df = pd.read_csv(atlas_path, index_col=0)

    probes = set(df.index)
    print(f"  Atlas contains {len(probes):,} probes")
    return probes


def map_probes_to_nanopore(probe_coords, methylation_data, atlas_probes, min_coverage=1):
    """
    Map EPIC probes to Nanopore methylation data.

    Returns dict: {probe_id: beta_value}
    """
    mapped_values = {}
    matched = 0
    low_coverage = 0
    not_in_epic = 0
    not_in_nanopore = 0

    for probe_id in atlas_probes:
        if probe_id not in probe_coords:
            not_in_epic += 1
            continue

        chrom, start, end, strand = probe_coords[probe_id]

        # Try to find matching position in methylation data
        # BED is 0-based, half-open [start, end)
        key = (chrom, start)

        if key in methylation_data:
            percent_meth, coverage = methylation_data[key]
            if coverage >= min_coverage:
                beta = percent_meth / 100.0
                mapped_values[probe_id] = beta
                matched += 1
            else:
                low_coverage += 1
        else:
            not_in_nanopore += 1

    print(f"    Matched: {matched:,}, Low cov: {low_coverage:,}, Not found: {not_in_nanopore:,} ({100.0 * matched / len(atlas_probes):.1f}%)")

    return mapped_values


def process_single_sample(bedmethyl_path, sample_name, probe_coords, atlas_probes, min_coverage):
    """Process a single sample and return beta values."""
    methylation_data = parse_bedmethyl(bedmethyl_path)
    mapped_values = map_probes_to_nanopore(probe_coords, methylation_data, atlas_probes, min_coverage)
    return mapped_values


def main():
    parser = argparse.ArgumentParser(
        description='Generate beta values CSV from Nanopore bedMethyl for deconvolution (v2 - batch processing)'
    )

    # Input options (mutually exclusive groups)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--bedmethyl', help='Single input bedMethyl file')
    input_group.add_argument('--input_dir', help='Directory containing bedMethyl files (use with --sample_ids)')

    parser.add_argument('--sample_ids', help='Text file with sample IDs (one per line), required with --input_dir')
    parser.add_argument('--epic_bed',
                        default='/home/chbope/extension/nWGS_manuscript_data/data/reference/EPIC_sites_NEW.bed',
                        help='EPIC probe coordinates BED file')
    parser.add_argument('--atlas',
                        default='/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/reference_atlas.csv',
                        help='Reference atlas CSV file')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    parser.add_argument('--sample_name', help='Sample name (for single sample mode)')
    parser.add_argument('--min_coverage', type=int, default=1,
                        help='Minimum read coverage per CpG (default: 1)')
    parser.add_argument('--combine', action='store_true',
                        help='Combine all samples into a single CSV file')

    args = parser.parse_args()

    # Validate arguments
    if args.input_dir and not args.sample_ids:
        parser.error("--sample_ids is required when using --input_dir")

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("Generate Beta Values for Deconvolution (v2)")
    print("=" * 60)

    # Load reference data
    print("\nStep 1: Loading reference data")
    print("-" * 40)
    epic_coords = load_epic_coordinates(args.epic_bed)
    atlas_probes = load_atlas_probes(args.atlas)

    # Determine samples to process
    if args.bedmethyl:
        # Single sample mode
        if args.sample_name:
            sample_name = args.sample_name
        else:
            sample_name = os.path.basename(args.bedmethyl)
            for suffix in ['.wf_mods.bedmethyl.gz', '.bedmethyl.gz', '.bedmethyl', '.gz']:
                sample_name = sample_name.replace(suffix, '')

        samples = [(sample_name, args.bedmethyl)]
    else:
        # Batch mode
        with open(args.sample_ids, 'r') as f:
            sample_ids = [line.strip() for line in f if line.strip()]
        # Remove duplicates while preserving order
        sample_ids = list(dict.fromkeys(sample_ids))

        samples = []
        for sample_id in sample_ids:
            bedmethyl_path = os.path.join(args.input_dir, f"{sample_id}.wf_mods.bedmethyl.gz")
            if os.path.exists(bedmethyl_path):
                samples.append((sample_id, bedmethyl_path))
            else:
                print(f"  WARNING: File not found for {sample_id}: {bedmethyl_path}")

    print(f"\nSamples to process: {len(samples)}")

    # Process samples
    print("\nStep 2: Processing samples")
    print("-" * 40)

    all_results = {}
    for i, (sample_name, bedmethyl_path) in enumerate(samples, 1):
        print(f"\n[{i}/{len(samples)}] {sample_name}")
        mapped_values = process_single_sample(
            bedmethyl_path, sample_name, epic_coords, atlas_probes, args.min_coverage
        )
        all_results[sample_name] = mapped_values

        # Save individual file
        df = pd.DataFrame({sample_name: mapped_values})
        df.index.name = ''

        # Save full file
        output_file = os.path.join(args.output_dir, f"{sample_name}_deconv_betas.csv")
        df.to_csv(output_file)

        # Save non-zero file
        df_nonzero = df[df[sample_name] > 0]
        output_file_nonzero = os.path.join(args.output_dir, f"{sample_name}_deconv_betas_nonzero.csv")
        df_nonzero.to_csv(output_file_nonzero)

        print(f"    Saved: {sample_name}_deconv_betas.csv ({len(mapped_values)} probes, {len(df_nonzero)} non-zero)")

    # Combine all samples if requested
    if args.combine and len(samples) > 1:
        print("\nStep 3: Combining all samples")
        print("-" * 40)

        combined_df = pd.DataFrame(all_results)
        combined_df.index.name = 'CpGs'

        # Save combined file (format: CpGs,sample1,sample2,...)
        combined_file = os.path.join(args.output_dir, "combined_deconv_betas.csv")
        combined_df.to_csv(combined_file)
        print(f"Saved combined file: {combined_file}")
        print(f"  Samples: {len(samples)}")
        print(f"  Probes: {len(combined_df)}")

        # Save combined non-zero file (probes with at least one non-zero value)
        combined_nonzero = combined_df[(combined_df > 0).any(axis=1)]
        combined_nonzero_file = os.path.join(args.output_dir, "combined_deconv_betas_nonzero.csv")
        combined_nonzero.to_csv(combined_nonzero_file)
        print(f"Saved combined non-zero file: {combined_nonzero_file}")
        print(f"  Probes with non-zero values: {len(combined_nonzero)}")

    print("\n" + "=" * 60)
    print(f"Done! Processed {len(samples)} samples")
    print("=" * 60)


if __name__ == '__main__':
    main()
