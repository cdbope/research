#!/usr/bin/env python3
"""
Generate beta values CSV file from Nanopore bedMethyl data for deconvolution.
5mC-only version: Extracts only 5mC (m) modification, ignoring 5hmC (h).

This script:
1. Reads reference atlas to get required probe IDs
2. Reads EPIC probe coordinates for those probes
3. Reads bedMethyl files, ONLY keeping 5mC (m) entries at positions matching atlas probes
4. Outputs beta values CSV files compatible with deconvolve.py

Key difference from merged versions:
- Only processes modification code 'm' (5mC methylation)
- Ignores modification code 'h' (5hmC hydroxymethylation)
- Use this when you need pure 5mC signal without 5hmC contribution

Usage:
    # Single sample
    python generate_beta_for_deconv_5mc.py --bedmethyl /path/to/sample.wf_mods.bedmethyl.gz --output_dir results/

    # Batch processing with 8 CPUs
    python generate_beta_for_deconv_5mc.py --input_dir /path/to/bedmethyl_files --sample_ids samples.txt --output_dir results/ --cpus 8

Author: Generated for Drexler pipeline
"""

import argparse
import gzip
import os
import sys
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool, cpu_count


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


def build_target_positions(atlas_probes, probe_coords):
    """
    Build a set of target positions (chr, start) from atlas probes.

    Returns:
        - target_positions: set of (chr, start) tuples to look for in bedMethyl
        - position_to_probe: dict mapping (chr, start) -> probe_id
    """
    target_positions = set()
    position_to_probe = {}

    for probe_id in atlas_probes:
        if probe_id in probe_coords:
            chrom, start, end, strand = probe_coords[probe_id]
            key = (chrom, start)
            target_positions.add(key)
            position_to_probe[key] = probe_id

    print(f"  Target positions to extract: {len(target_positions):,}")
    return target_positions, position_to_probe


def extract_5mc_only_targeted(bedmethyl_path, target_positions, position_to_probe, min_coverage=1):
    """
    Parse bedMethyl file and extract ONLY 5mC (m) entries at positions matching target probes.

    This extracts pure 5mC methylation signal, ignoring 5hmC hydroxymethylation.

    Args:
        bedmethyl_path: Input bedMethyl file
        target_positions: Set of (chr, start) positions to extract
        position_to_probe: Dict mapping (chr, start) -> probe_id
        min_coverage: Minimum coverage threshold

    Returns dict: {probe_id: beta_value}
    """
    print(f"  Extracting 5mC-only from: {os.path.basename(bedmethyl_path)}")

    # Collect only 5mC (m) entries for target positions
    position_data = {}

    open_func = gzip.open if bedmethyl_path.endswith('.gz') else open
    mode = 'rt' if bedmethyl_path.endswith('.gz') else 'r'

    lines_read = 0
    m_entries_found = 0
    h_entries_skipped = 0

    with open_func(bedmethyl_path, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue

            lines_read += 1

            fields = line.strip().split('\t')
            if len(fields) < 18:
                continue

            chrom = fields[0]
            start = int(fields[1])
            mod_code = fields[3]

            key = (chrom, start)

            # Only process if this position is in our target set
            if key in target_positions:
                if mod_code == 'm':
                    # Only keep 5mC entries
                    position_data[key] = fields
                    m_entries_found += 1
                elif mod_code == 'h':
                    # Skip 5hmC entries
                    h_entries_skipped += 1

    print(f"    Read {lines_read:,} lines")
    print(f"    5mC (m) entries found: {m_entries_found:,}")
    print(f"    5hmC (h) entries skipped: {h_entries_skipped:,}")

    # Calculate beta values from 5mC entries only
    mapped_values = {}
    low_coverage_count = 0

    for key, fields in position_data.items():
        probe_id = position_to_probe[key]

        n_valid = int(fields[9])
        percent_meth = float(fields[10])

        if n_valid >= min_coverage:
            beta = percent_meth / 100.0
            mapped_values[probe_id] = beta
        else:
            low_coverage_count += 1

    print(f"    Mapped probes: {len(mapped_values):,}, Low coverage skipped: {low_coverage_count:,}")

    return mapped_values


def process_single_sample(sample_info, target_positions, position_to_probe, min_coverage, output_dir):
    """
    Process a single sample: extract 5mC-only from targeted positions.

    This function is designed to be called in parallel.
    """
    sample_name, bedmethyl_path = sample_info

    print(f"\nProcessing: {sample_name}")

    # Extract 5mC-only from targeted positions
    mapped_values = extract_5mc_only_targeted(
        bedmethyl_path, target_positions, position_to_probe, min_coverage
    )

    # Save individual file
    df = pd.DataFrame({sample_name: mapped_values})
    df.index.name = ''

    output_file = os.path.join(output_dir, f"{sample_name}_deconv_betas_5mc.csv")
    df.to_csv(output_file)

    print(f"    Saved: {sample_name}_deconv_betas_5mc.csv ({len(mapped_values)} probes)")

    return sample_name, mapped_values


def process_sample_wrapper(args):
    """Wrapper function for multiprocessing Pool.map()"""
    sample_info, target_positions, position_to_probe, min_coverage, output_dir = args
    return process_single_sample(sample_info, target_positions, position_to_probe, min_coverage, output_dir)


def main():
    parser = argparse.ArgumentParser(
        description='Generate beta values CSV from Nanopore bedMethyl for deconvolution (5mC only, no 5hmC)'
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
    parser.add_argument('--cpus', type=int, default=1,
                        help=f'Number of CPUs to use for parallel processing (default: 1, max available: {cpu_count()})')

    args = parser.parse_args()

    # Validate arguments
    if args.input_dir and not args.sample_ids:
        parser.error("--sample_ids is required when using --input_dir")

    # Validate CPU count
    max_cpus = cpu_count()
    if args.cpus < 1:
        args.cpus = 1
    elif args.cpus > max_cpus:
        print(f"Warning: Requested {args.cpus} CPUs but only {max_cpus} available. Using {max_cpus}.")
        args.cpus = max_cpus

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("Generate Beta Values for Deconvolution (5mC ONLY)")
    print("=" * 60)
    print(f"Using {args.cpus} CPU(s)")
    print("NOTE: This extracts only 5mC (m) entries, skipping 5hmC (h)")

    # Load reference data FIRST (this is the optimization)
    print("\nStep 1: Loading reference data and building target positions")
    print("-" * 40)
    atlas_probes = load_atlas_probes(args.atlas)
    epic_coords = load_epic_coordinates(args.epic_bed)

    # Build set of target positions to look for
    target_positions, position_to_probe = build_target_positions(atlas_probes, epic_coords)

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
    print("\nStep 2: Processing samples (5mC-only extraction)")
    print("-" * 40)

    all_results = {}

    if args.cpus == 1 or len(samples) == 1:
        # Sequential processing
        for i, sample_info in enumerate(samples, 1):
            print(f"\n[{i}/{len(samples)}] {sample_info[0]}")
            sample_name, mapped_values = process_single_sample(
                sample_info, target_positions, position_to_probe, args.min_coverage, args.output_dir
            )
            all_results[sample_name] = mapped_values
    else:
        # Parallel processing
        print(f"\nUsing {args.cpus} parallel workers")

        # Prepare arguments for each sample
        process_args = [
            (sample_info, target_positions, position_to_probe, args.min_coverage, args.output_dir)
            for sample_info in samples
        ]

        with Pool(processes=args.cpus) as pool:
            results = pool.map(process_sample_wrapper, process_args)

        for sample_name, mapped_values in results:
            all_results[sample_name] = mapped_values

    # Combine all samples if requested
    if args.combine and len(samples) > 1:
        print("\nStep 3: Combining all samples")
        print("-" * 40)

        combined_df = pd.DataFrame(all_results)
        combined_df.index.name = 'CpGs'

        # Save combined file (format: CpGs,sample1,sample2,...)
        combined_file = os.path.join(args.output_dir, "combined_deconv_betas_5mc.csv")
        combined_df.to_csv(combined_file)
        print(f"Saved combined file: {combined_file}")
        print(f"  Samples: {len(samples)}")
        print(f"  Probes: {len(combined_df)}")

    print("\n" + "=" * 60)
    print(f"Done! Processed {len(samples)} samples (5mC only)")
    print("=" * 60)


if __name__ == '__main__':
    main()
