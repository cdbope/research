#!/usr/bin/env python3
"""
Generate beta values CSV file from Nanopore bedMethyl data for deconvolution.

This script:
1. Reads a bedMethyl file (Nanopore methylation data)
2. Reads EPIC probe coordinates (probe ID -> chr:position)
3. Reads reference atlas to get required probe IDs
4. Maps Nanopore methylation to EPIC probes
5. Outputs a beta values CSV file compatible with deconvolve.py

Usage:
    python generate_beta_for_deconv.py <bedmethyl_file> [options]

Example:
    python generate_beta_for_deconv.py /path/to/sample.wf_mods.bedmethyl.gz \
        --epic_bed /path/to/EPIC_sites_NEW.bed \
        --atlas /path/to/full_atlas.csv.gz \
        --output_dir results/ \
        --sample_name T22-116

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
    print(f"Parsing bedMethyl file: {bedmethyl_path}")

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

            if line_num % 1000000 == 0:
                print(f"  Processed {line_num:,} lines...")

    print(f"  Total unique CpG positions: {len(methylation_data):,}")
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

    df = pd.read_csv(atlas_path, index_col=0, nrows=0)  # Just read header to get index
    df = pd.read_csv(atlas_path, index_col=0)

    probes = set(df.index)
    print(f"  Atlas contains {len(probes):,} probes")
    return probes


def map_probes_to_nanopore(probe_coords, methylation_data, atlas_probes, min_coverage=1):
    """
    Map EPIC probes to Nanopore methylation data.

    Returns dict: {probe_id: beta_value}
    """
    print(f"Mapping probes to Nanopore data (min coverage: {min_coverage}x)...")

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

    print(f"  Matched: {matched:,}")
    print(f"  Low coverage (< {min_coverage}x): {low_coverage:,}")
    print(f"  Not in EPIC BED: {not_in_epic:,}")
    print(f"  Not in Nanopore data: {not_in_nanopore:,}")
    print(f"  Coverage: {100.0 * matched / len(atlas_probes):.1f}%")

    return mapped_values


def main():
    parser = argparse.ArgumentParser(
        description='Generate beta values CSV from Nanopore bedMethyl for deconvolution'
    )
    parser.add_argument('bedmethyl', help='Input bedMethyl file (from modkit)')
    parser.add_argument('--epic_bed',
                        default='/home/chbope/extension/nWGS_manuscript_data/data/reference/EPIC_sites_NEW.bed',
                        help='EPIC probe coordinates BED file')
    parser.add_argument('--atlas',
                        default='/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master/reference_atlas.csv',
                        help='Reference atlas CSV file')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    parser.add_argument('--sample_name', help='Sample name (default: derived from filename)')
    parser.add_argument('--min_coverage', type=int, default=1,
                        help='Minimum read coverage per CpG (default: 1)')

    args = parser.parse_args()

    # Determine sample name
    if args.sample_name:
        sample_name = args.sample_name
    else:
        sample_name = os.path.basename(args.bedmethyl)
        for suffix in ['.wf_mods.bedmethyl.gz', '.bedmethyl.gz', '.bedmethyl', '.gz']:
            sample_name = sample_name.replace(suffix, '')

    print("=" * 60)
    print("Generate Beta Values for Deconvolution")
    print("=" * 60)
    print(f"Sample: {sample_name}")
    print(f"Input: {args.bedmethyl}")
    print(f"EPIC BED: {args.epic_bed}")
    print(f"Atlas: {args.atlas}")
    print(f"Output directory: {args.output_dir}")
    print(f"Minimum coverage: {args.min_coverage}x")

    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    print("\n" + "=" * 60)
    print("Step 1: Loading reference data")
    print("=" * 60)
    epic_coords = load_epic_coordinates(args.epic_bed)
    atlas_probes = load_atlas_probes(args.atlas)

    # Parse bedMethyl
    print("\n" + "=" * 60)
    print("Step 2: Parsing Nanopore methylation data")
    print("=" * 60)
    methylation_data = parse_bedmethyl(args.bedmethyl)

    # Map probes
    print("\n" + "=" * 60)
    print("Step 3: Mapping probes to Nanopore data")
    print("=" * 60)
    mapped_values = map_probes_to_nanopore(
        epic_coords, methylation_data, atlas_probes,
        min_coverage=args.min_coverage
    )

    # Create output DataFrame
    print("\n" + "=" * 60)
    print("Step 4: Saving beta values")
    print("=" * 60)

    df = pd.DataFrame({sample_name: mapped_values})
    df.index.name = ''

    # Save full file (including zeros)
    output_file = os.path.join(args.output_dir, f"{sample_name}_deconv_betas.csv")
    df.to_csv(output_file)
    print(f"Saved: {output_file}")
    print(f"Probes with values: {len(mapped_values):,} / {len(atlas_probes):,}")

    # Save filtered file (excluding zeros)
    df_nonzero = df[df[sample_name] > 0]
    output_file_nonzero = os.path.join(args.output_dir, f"{sample_name}_deconv_betas_nonzero.csv")
    df_nonzero.to_csv(output_file_nonzero)
    print(f"Saved (non-zero only): {output_file_nonzero}")
    print(f"Probes with non-zero values: {len(df_nonzero):,}")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)

    return df


if __name__ == '__main__':
    main()
