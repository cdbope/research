#!/usr/bin/env python3
"""
Convert Nanopore bedMethyl data to Drexler pipeline format.

This script maps Nanopore methylation calls (from modkit) to Illumina 450K CpG probe positions,
allowing you to use Nanopore data with the neural glioblastoma classification pipeline.

Usage:
    python nanopore_to_drexler.py <bedmethyl_file> <output_csv> [--probe_coords <probe_coords.csv>]

Example:
    python nanopore_to_drexler.py T001.wf_mods.bedmethyl.gz betas.csv

Author: Generated for Drexler et al. pipeline adaptation
"""

import argparse
import gzip
import pandas as pd
import numpy as np
import os
import sys

def load_probe_coordinates(probe_coords_file=None, genome='hg38'):
    """
    Load Illumina 450K probe coordinates.
    If probe_coords_file is not provided, generate from R/Bioconductor.

    Args:
        probe_coords_file: Pre-computed probe coordinates CSV file
        genome: Reference genome build ('hg19' or 'hg38')
    """
    # Check for cached file based on genome build
    cache_file = f"probe_coordinates_450k_{genome}.csv"

    if probe_coords_file and os.path.exists(probe_coords_file):
        print(f"Loading probe coordinates from {probe_coords_file}")
        df = pd.read_csv(probe_coords_file, index_col=0)
        return df

    if os.path.exists(cache_file):
        print(f"Loading cached probe coordinates ({genome}): {cache_file}")
        df = pd.read_csv(cache_file, index_col=0)
        return df

    # Generate probe coordinates using R
    print(f"Generating probe coordinates for {genome}...")

    if genome == 'hg19':
        r_script = '''
        library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        data(Locations)

        # Convert to data frame
        locs <- as.data.frame(Locations)
        locs$probe_id <- rownames(locs)

        # Write to CSV
        write.csv(locs, "probe_coordinates_450k_hg19.csv", row.names=FALSE)
        cat("Saved probe coordinates (hg19)\\n")
        '''
    elif genome == 'hg38':
        r_script = '''
        library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        library(rtracklayer)
        library(GenomicRanges)

        cat("Loading Illumina 450K probe coordinates (hg19)...\\n")
        data(Locations)

        # Convert to GRanges
        locs <- as.data.frame(Locations)
        locs$probe_id <- rownames(locs)

        gr_hg19 <- GRanges(
            seqnames = locs$chr,
            ranges = IRanges(start = locs$pos, width = 1),
            strand = locs$strand,
            probe_id = locs$probe_id
        )

        cat(sprintf("Loaded %d probes\\n", length(gr_hg19)))

        # Download and decompress chain file if not present
        chain_file_gz <- "hg19ToHg38.over.chain.gz"
        chain_file <- "hg19ToHg38.over.chain"

        if (!file.exists(chain_file)) {
            if (!file.exists(chain_file_gz)) {
                cat("Downloading hg19ToHg38 chain file...\\n")
                download.file(
                    "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                    chain_file_gz,
                    mode = "wb"
                )
            }
            cat("Decompressing chain file...\\n")
            system(paste("gunzip -k", chain_file_gz))
        }

        # Load chain file
        cat("Loading chain file...\\n")
        chain <- import.chain(chain_file)

        # Perform liftOver
        cat("Performing liftOver hg19 -> hg38...\\n")
        gr_hg38_list <- liftOver(gr_hg19, chain)

        # Convert to data frame (some probes may not lift over)
        lifted_idx <- which(lengths(gr_hg38_list) == 1)
        gr_hg38 <- unlist(gr_hg38_list[lifted_idx])

        cat(sprintf("Successfully lifted %d/%d probes (%.1f%%)\\n",
            length(gr_hg38), length(gr_hg19),
            100 * length(gr_hg38) / length(gr_hg19)))

        # Create output data frame
        result <- data.frame(
            probe_id = gr_hg38$probe_id,
            chr = as.character(seqnames(gr_hg38)),
            pos = start(gr_hg38),
            strand = as.character(strand(gr_hg38)),
            stringsAsFactors = FALSE
        )

        # Save to CSV
        write.csv(result, "probe_coordinates_450k_hg38.csv", row.names = FALSE)
        cat("Saved probe coordinates (hg38)\\n")
        '''
    else:
        print(f"Error: Unknown genome build '{genome}'. Use 'hg19' or 'hg38'.")
        sys.exit(1)

    import subprocess
    result = subprocess.run(['Rscript', '-e', r_script], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error generating probe coordinates: {result.stderr}")
        sys.exit(1)

    print(result.stdout)

    df = pd.read_csv(cache_file)
    df = df.set_index('probe_id')
    return df


def load_required_probes(base_path):
    """
    Load the list of CpG probes required for the Drexler classifier.
    """
    low_cgs = pd.read_csv(os.path.join(base_path, "code/neural_group_classification/low_manifest_all.csv"),
                          index_col=0).index.tolist()
    high_cgs = pd.read_csv(os.path.join(base_path, "code/neural_group_classification/high_manifest_all.csv"),
                           index_col=0).index.tolist()
    return low_cgs + high_cgs


def parse_bedmethyl(bedmethyl_file):
    """
    Parse modkit bedMethyl file and extract 5mC methylation values.

    bedMethyl format (from modkit):
    col 0: chrom
    col 1: start (0-based)
    col 2: end
    col 3: modification code (m=5mC, h=5hmC)
    col 4: score
    col 5: strand
    col 10: percent methylated
    col 11: N_modified
    col 12: N_canonical
    col 13: N_other_mod
    col 14: N_delete
    col 15: N_fail
    col 16: N_diff
    col 17: N_nocall
    """
    print(f"Parsing bedMethyl file: {bedmethyl_file}")

    # Determine if gzipped
    open_func = gzip.open if bedmethyl_file.endswith('.gz') else open

    methylation_data = {}

    with open_func(bedmethyl_file, 'rt') as f:
        for line_num, line in enumerate(f):
            if line_num % 5000000 == 0 and line_num > 0:
                print(f"  Processed {line_num:,} lines...")

            fields = line.strip().split('\t')

            # Only use 5mC (m) calls, skip 5hmC (h) and others
            mod_code = fields[3]
            if mod_code != 'm':
                continue

            chrom = fields[0]
            # Ensure chr prefix
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom

            pos = int(fields[1]) + 1  # Convert 0-based to 1-based
            strand = fields[5]
            percent_meth = float(fields[10])

            # Coverage calculation
            n_mod = int(fields[11])
            n_canonical = int(fields[12])
            coverage = n_mod + n_canonical

            # Create key: chr:pos:strand
            key = f"{chrom}:{pos}:{strand}"

            # Store methylation percentage (as beta value: 0-1 scale)
            methylation_data[key] = {
                'beta': percent_meth / 100.0,
                'coverage': coverage
            }

    print(f"  Loaded {len(methylation_data):,} 5mC methylation calls")
    return methylation_data


def map_probes_to_nanopore(probe_coords, methylation_data, min_coverage=5):
    """
    Map Illumina probes to Nanopore methylation values.

    For each probe, look up the corresponding position in the Nanopore data.
    """
    print("Mapping probes to Nanopore methylation values...")

    mapped_values = {}
    matched = 0
    unmatched = 0
    low_coverage = 0

    for probe_id, row in probe_coords.iterrows():
        chrom = row['chr']
        pos = row['pos']
        strand = row['strand']

        # Try multiple strand options since:
        # - CpG sites are typically measured on both strands
        # - Some bedMethyl files use '.' for unstranded data
        key_plus = f"{chrom}:{pos}:+"
        key_minus = f"{chrom}:{pos}:-"
        key_unstranded = f"{chrom}:{pos}:."

        value = None
        coverage = 0

        # Check the specified strand first
        if strand == '+' and key_plus in methylation_data:
            value = methylation_data[key_plus]['beta']
            coverage = methylation_data[key_plus]['coverage']
        elif strand == '-' and key_minus in methylation_data:
            value = methylation_data[key_minus]['beta']
            coverage = methylation_data[key_minus]['coverage']
        # Then try unstranded data (common in Nanopore output)
        elif key_unstranded in methylation_data:
            value = methylation_data[key_unstranded]['beta']
            coverage = methylation_data[key_unstranded]['coverage']
        # Then try the opposite strand
        elif key_plus in methylation_data:
            value = methylation_data[key_plus]['beta']
            coverage = methylation_data[key_plus]['coverage']
        elif key_minus in methylation_data:
            value = methylation_data[key_minus]['beta']
            coverage = methylation_data[key_minus]['coverage']

        if value is not None:
            if coverage >= min_coverage:
                mapped_values[probe_id] = value
                matched += 1
            else:
                mapped_values[probe_id] = np.nan
                low_coverage += 1
        else:
            mapped_values[probe_id] = np.nan
            unmatched += 1

    print(f"  Matched: {matched:,} probes")
    print(f"  Low coverage (<{min_coverage}x): {low_coverage:,} probes")
    print(f"  Not found: {unmatched:,} probes")

    return mapped_values


def main():
    parser = argparse.ArgumentParser(
        description='Convert Nanopore bedMethyl to Drexler pipeline format'
    )
    parser.add_argument('bedmethyl', help='Input bedMethyl file (gzipped or plain)')
    parser.add_argument('output', help='Output CSV file (betas format)')
    parser.add_argument('--probe_coords', help='Pre-computed probe coordinates CSV')
    parser.add_argument('--min_coverage', type=int, default=5,
                        help='Minimum coverage for a probe (default: 5)')
    parser.add_argument('--sample_name', help='Sample name (default: derived from filename)')
    parser.add_argument('--base_path',
                        default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        help='Base path to the Drexler code directory')

    args = parser.parse_args()

    # Determine sample name
    if args.sample_name:
        sample_name = args.sample_name
    else:
        sample_name = os.path.basename(args.bedmethyl).replace('.wf_mods.bedmethyl.gz', '').replace('.bedmethyl.gz', '').replace('.bedmethyl', '')

    print(f"Sample name: {sample_name}")
    print(f"Base path: {args.base_path}")

    # Load probe coordinates
    probe_coords = load_probe_coordinates(args.probe_coords)
    print(f"Loaded {len(probe_coords):,} probe coordinates")

    # Load required probes for classifier
    required_probes = load_required_probes(args.base_path)
    print(f"Classifier requires {len(required_probes):,} probes")

    # Filter to required probes
    probe_coords_required = probe_coords[probe_coords.index.isin(required_probes)]
    print(f"Found coordinates for {len(probe_coords_required):,} required probes")

    # Parse bedMethyl
    methylation_data = parse_bedmethyl(args.bedmethyl)

    # Map probes
    mapped_values = map_probes_to_nanopore(probe_coords_required, methylation_data,
                                            min_coverage=args.min_coverage)

    # Create output DataFrame in the format expected by run_example.ipynb
    # The betas.csv format has CpG IDs as rows and samples as columns
    df_out = pd.DataFrame(mapped_values, index=[sample_name]).T
    df_out.index.name = None

    # Save
    df_out.to_csv(args.output)
    print(f"\nSaved to {args.output}")

    # Report coverage of required probes
    n_valid = df_out[sample_name].notna().sum()
    n_missing = df_out[sample_name].isna().sum()
    coverage_pct = 100.0 * n_valid / len(df_out)

    print(f"\nProbe coverage summary:")
    print(f"  Valid values: {n_valid:,} ({coverage_pct:.1f}%)")
    print(f"  Missing: {n_missing:,} ({100-coverage_pct:.1f}%)")

    if coverage_pct < 80:
        print("\nWARNING: Low probe coverage may affect classification accuracy.")
        print("Consider using the imputer from the Drexler pipeline.")


if __name__ == '__main__':
    main()
