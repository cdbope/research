#!/usr/bin/env python3
"""
Script to merge GFF/GTF annotation files by adding NM_ID from MANE Select transcripts.
"""

import re
import sys
from collections import defaultdict


def parse_attributes(attr_string):
    """Parse the attributes column (column 9) into a dictionary."""
    attributes = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        if '=' in item:
            key, value = item.split('=', 1)
            attributes[key] = value
    return attributes


def extract_nm_id(transcript_id):
    """Extract clean NM_ID from transcript_id like 'NM_002168.4' -> 'NM_002168'"""
    if transcript_id:
        # Remove version number after the dot
        match = re.match(r'(NM_\d+)', transcript_id)
        if match:
            return match.group(1)
    return None


def read_mane_select_file(filename):
    """
    Read the first file (BestRefSeq) and extract MANE Select NM_IDs mapped to gene names.
    Returns a dictionary: {gene_name: nm_id}
    """
    gene_to_nm = {}

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 9:
                    continue

                # Check if it's an mRNA entry
                feature_type = fields[2]
                if feature_type != 'mRNA':
                    continue

                attributes = parse_attributes(fields[8])

                # Check for MANE Select tag
                tag = attributes.get('tag', '')
                if 'MANE Select' not in tag:
                    continue

                # Extract gene name and transcript ID
                gene_name = attributes.get('gene', '')
                transcript_id = attributes.get('transcript_id', '')

                if gene_name and transcript_id:
                    nm_id = extract_nm_id(transcript_id)
                    if nm_id:
                        gene_to_nm[gene_name] = nm_id
                        print(f"Found MANE Select: {gene_name} -> {nm_id}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found", file=sys.stderr)
        sys.exit(1)

    return gene_to_nm


def process_gene_file(input_filename, output_filename, gene_to_nm):
    """
    Read the second file (gene annotations) and add NM_ID column.
    """
    processed_count = 0
    matched_count = 0

    try:
        with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
            for line in infile:
                line = line.strip()
                if not line or line.startswith('#'):
                    outfile.write(line + '\n')
                    continue

                fields = line.split('\t')
                if len(fields) < 9:
                    # If line doesn't have enough fields, write as is
                    outfile.write(line + '\n')
                    continue

                processed_count += 1

                # Parse attributes to find gene_name (field 9 in this format)
                # Handle files with attributes in different columns
                attr_column = 9 if len(fields) > 9 else 8
                attributes = parse_attributes(fields[attr_column])
                gene_name = attributes.get('gene_name', '')

                # Look up NM_ID
                nm_id = gene_to_nm.get(gene_name, 'NA')
                if nm_id != 'NA':
                    matched_count += 1

                # Replace the last column if it exists, otherwise add new column
                if len(fields) > 10:
                    fields[-1] = nm_id
                    output_line = '\t'.join(fields)
                else:
                    output_line = '\t'.join(fields) + '\t' + nm_id
                outfile.write(output_line + '\n')

        print(f"\nProcessing complete:", file=sys.stderr)
        print(f"  Total lines processed: {processed_count}", file=sys.stderr)
        print(f"  Lines with matched NM_ID: {matched_count}", file=sys.stderr)
        print(f"  Output written to: {output_filename}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found", file=sys.stderr)
        sys.exit(1)


def main():
    if len(sys.argv) != 4:
        print("Usage: python merge_annotations.py <mane_select_file> <gene_file> <output_file>")
        print("\nExample:")
        print("  python merge_annotations.py refseq_mane.gff genes.gtf output.tsv")
        sys.exit(1)

    mane_select_file = sys.argv[1]
    gene_file = sys.argv[2]
    output_file = sys.argv[3]

    print("Step 1: Reading MANE Select transcripts...", file=sys.stderr)
    gene_to_nm = read_mane_select_file(mane_select_file)
    print(f"Found {len(gene_to_nm)} MANE Select genes\n", file=sys.stderr)

    print("Step 2: Processing gene file and adding NM_IDs...", file=sys.stderr)
    process_gene_file(gene_file, output_file, gene_to_nm)


if __name__ == "__main__":
    main()
