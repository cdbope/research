#!/usr/bin/env python3
"""
Generate beta values CSV file from Nanopore bedMethyl data for deconvolution.
Version 4: Merges 5mC (m) and 5hmC (h) modifications with multi-CPU support.

This script:
1. Reads bedMethyl files (Nanopore methylation data) for multiple samples
2. MERGES 5mC and 5hmC entries at the same position (sums N_modified counts)
3. Saves reformatted bedMethyl files with merged modifications
4. Reads EPIC probe coordinates (probe ID -> chr:position)
5. Reads reference atlas to get required probe IDs
6. Maps Nanopore methylation to EPIC probes
7. Outputs beta values CSV files compatible with deconvolve.py

Changes from v3:
- Added --cpus parameter for parallel processing of multiple samples
- Removed non-zero file generation

The merge logic:
- For each position, sum N_modified from both 5mC (m) and 5hmC (h) entries
- Recalculate percent_methylated as: (N_mod_5mC + N_mod_5hmC) / N_valid * 100
- Output as 'm' modification type

Example:
    Input:
    chr1    13692   13693   h       34      .       13692   13693   255,0,0 34      2.94    1       24      9       0       2       0       2
    chr1    13692   13693   m       34      .       13692   13693   255,0,0 34      26.47   9       24      1       0       2       0       2

    Merged output:
    chr1    13692   13693   m       34      .       13692   13693   255,0,0 34      29.41   10      24      10      0       2       0       2

Usage:
    # Single sample
    python generate_beta_for_deconv_v4.py --bedmethyl /path/to/sample.wf_mods.bedmethyl.gz --output_dir results/

    # Batch processing with 8 CPUs
    python generate_beta_for_deconv_v4.py --input_dir /path/to/bedmethyl_files --sample_ids samples.txt --output_dir results/ --cpus 8

Author: Generated for Drexler pipeline
"""

import argparse
import gzip
import os
import sys
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from functools import partial


def merge_5mc_5hmc_bedmethyl(bedmethyl_path, output_path):
    """
    Parse bedMethyl file and merge 5mC (m) and 5hmC (h) entries at the same position.

    bedMethyl format (from modkit):
    col 0: chrom
    col 1: start (0-based)
    col 2: end
    col 3: modification code (m=5mC, h=5hmC)
    col 4: score
    col 5: strand
    col 6-8: thick start/end, color
    col 9: N_valid (coverage)
    col 10: percent methylated
    col 11: N_modified
    col 12: N_canonical
    col 13: N_other_mod
    col 14: N_delete
    col 15: N_fail
    col 16: N_diff
    col 17: N_nocall

    Returns dict: {(chr, pos): (percent_methylated, coverage)}
    """
    print(f"  Parsing and merging 5mC+5hmC from: {os.path.basename(bedmethyl_path)}")

    # First pass: collect all entries by position
    # Key: (chrom, start, strand)
    # Value: dict with 'm' and 'h' entries
    position_data = defaultdict(lambda: {'m': None, 'h': None, 'other': []})

    open_func = gzip.open if bedmethyl_path.endswith('.gz') else open
    mode = 'rt' if bedmethyl_path.endswith('.gz') else 'r'

    with open_func(bedmethyl_path, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 18:
                continue

            chrom = fields[0]
            start = int(fields[1])
            mod_code = fields[3]
            strand = fields[5]

            key = (chrom, start, strand)

            if mod_code == 'm':
                position_data[key]['m'] = fields
            elif mod_code == 'h':
                position_data[key]['h'] = fields
            else:
                position_data[key]['other'].append(fields)

    # Second pass: merge and write output
    print(f"    Found {len(position_data):,} unique positions")

    merged_count = 0
    m_only_count = 0
    h_only_count = 0

    # Determine output compression
    if output_path.endswith('.gz'):
        out_func = gzip.open
        out_mode = 'wt'
    else:
        out_func = open
        out_mode = 'w'

    methylation_data = {}

    with out_func(output_path, out_mode) as out_f:
        for key in sorted(position_data.keys()):
            chrom, start, strand = key
            data = position_data[key]

            m_entry = data['m']
            h_entry = data['h']

            if m_entry and h_entry:
                # Merge 5mC and 5hmC
                merged_count += 1

                # Extract values from both entries
                n_valid = int(m_entry[9])  # Should be same for both
                n_mod_m = int(m_entry[11])
                n_mod_h = int(h_entry[11])
                n_canonical_m = int(m_entry[12])
                n_other_mod_m = int(m_entry[13])

                # Sum modifications
                n_mod_total = n_mod_m + n_mod_h

                # Recalculate percent methylated
                if n_valid > 0:
                    percent_meth = (n_mod_total / n_valid) * 100.0
                else:
                    percent_meth = 0.0

                # Create merged entry (use 'm' fields as template)
                merged_fields = m_entry.copy()
                merged_fields[10] = f"{percent_meth:.2f}"
                merged_fields[11] = str(n_mod_total)
                # N_other_mod becomes the 5hmC count (for reference)
                merged_fields[13] = str(n_mod_h)

                out_f.write('\t'.join(merged_fields) + '\n')

                # Store for beta extraction
                methylation_data[(chrom, start)] = (percent_meth, n_valid)

            elif m_entry:
                # Only 5mC entry exists
                m_only_count += 1
                out_f.write('\t'.join(m_entry) + '\n')

                n_valid = int(m_entry[9])
                percent_meth = float(m_entry[10])
                methylation_data[(chrom, start)] = (percent_meth, n_valid)

            elif h_entry:
                # Only 5hmC entry exists - treat as methylation
                h_only_count += 1

                # Convert h entry to m entry
                converted_fields = h_entry.copy()
                converted_fields[3] = 'm'  # Change mod code to 'm'

                out_f.write('\t'.join(converted_fields) + '\n')

                n_valid = int(h_entry[9])
                percent_meth = float(h_entry[10])
                methylation_data[(chrom, start)] = (percent_meth, n_valid)

            # Write other modification types as-is
            for other_entry in data['other']:
                out_f.write('\t'.join(other_entry) + '\n')

    print(f"    Merged (5mC+5hmC): {merged_count:,}")
    print(f"    5mC only: {m_only_count:,}")
    print(f"    5hmC only: {h_only_count:,}")
    print(f"    Saved merged bedMethyl to: {os.path.basename(output_path)}")

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


def process_single_sample(sample_info, probe_coords, atlas_probes, min_coverage, output_dir, keep_merged):
    """
    Process a single sample: merge 5mC+5hmC and return beta values.

    This function is designed to be called in parallel.
    """
    sample_name, bedmethyl_path = sample_info

    print(f"\nProcessing: {sample_name}")

    # Create merged bedMethyl file
    merged_bedmethyl_path = os.path.join(output_dir, f"{sample_name}_merged.bedmethyl.gz")
    methylation_data = merge_5mc_5hmc_bedmethyl(bedmethyl_path, merged_bedmethyl_path)

    # Map to probes
    mapped_values = map_probes_to_nanopore(probe_coords, methylation_data, atlas_probes, min_coverage)

    # Save individual file
    df = pd.DataFrame({sample_name: mapped_values})
    df.index.name = ''

    # Save full file
    output_file = os.path.join(output_dir, f"{sample_name}_deconv_betas_merged.csv")
    df.to_csv(output_file)

    print(f"    Saved: {sample_name}_deconv_betas_merged.csv ({len(mapped_values)} probes)")

    # Cleanup merged bedMethyl file if not keeping
    if not keep_merged and os.path.exists(merged_bedmethyl_path):
        os.remove(merged_bedmethyl_path)

    return sample_name, mapped_values


def process_sample_wrapper(args):
    """Wrapper function for multiprocessing Pool.map()"""
    sample_info, probe_coords, atlas_probes, min_coverage, output_dir, keep_merged = args
    return process_single_sample(sample_info, probe_coords, atlas_probes, min_coverage, output_dir, keep_merged)


def main():
    parser = argparse.ArgumentParser(
        description='Generate beta values CSV from Nanopore bedMethyl for deconvolution (v4 - merges 5mC+5hmC, multi-CPU)'
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
    parser.add_argument('--keep_merged', action='store_true',
                        help='Keep the merged bedMethyl files (default: delete after processing)')
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
    print("Generate Beta Values for Deconvolution (v4 - 5mC+5hmC merged)")
    print("=" * 60)
    print(f"Using {args.cpus} CPU(s)")

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
    print("\nStep 2: Processing samples (merging 5mC+5hmC)")
    print("-" * 40)

    all_results = {}

    if args.cpus == 1 or len(samples) == 1:
        # Sequential processing
        for i, sample_info in enumerate(samples, 1):
            print(f"\n[{i}/{len(samples)}] {sample_info[0]}")
            sample_name, mapped_values = process_single_sample(
                sample_info, epic_coords, atlas_probes, args.min_coverage, args.output_dir, args.keep_merged
            )
            all_results[sample_name] = mapped_values
    else:
        # Parallel processing
        print(f"\nUsing {args.cpus} parallel workers")

        # Prepare arguments for each sample
        process_args = [
            (sample_info, epic_coords, atlas_probes, args.min_coverage, args.output_dir, args.keep_merged)
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
        combined_file = os.path.join(args.output_dir, "combined_deconv_betas_merged.csv")
        combined_df.to_csv(combined_file)
        print(f"Saved combined file: {combined_file}")
        print(f"  Samples: {len(samples)}")
        print(f"  Probes: {len(combined_df)}")

    print("\n" + "=" * 60)
    print(f"Done! Processed {len(samples)} samples (5mC+5hmC merged)")
    print("=" * 60)


if __name__ == '__main__':
    main()
