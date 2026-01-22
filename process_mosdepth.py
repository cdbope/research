#!/usr/bin/env python3

import os
import sys
import glob
from pathlib import Path

def read_sample_ids(sample_id_file):
    """Read sample IDs from file."""
    with open(sample_id_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def get_coverage_from_summary(file_path):
    """Extract coverage value from the last line of mosdepth summary file."""
    try:
        with open(file_path, 'r') as f:
            # Read all lines and get the last one
            lines = f.readlines()
            if not lines:
                return None
            
            last_line = lines[-1].strip()
            # Split the line and get the 4th column (0-based index)
            columns = last_line.split('\t')
            if len(columns) >= 4:
                return float(columns[3])
    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}", file=sys.stderr)
    return None

def main():
    if len(sys.argv) != 4:
        print("Usage: python process_mosdepth.py <sample_id_file> <mosdepth_dir> <output_file>")
        sys.exit(1)

    sample_id_file = sys.argv[1]
    mosdepth_dir = sys.argv[2]
    output_file = sys.argv[3]

    # Read sample IDs
    sample_ids = read_sample_ids(sample_id_file)
    print(f"Found {len(sample_ids)} sample IDs")

    # Process each sample and collect results
    results = []
    for sample_id in sample_ids:
        # Look for the mosdepth summary file
        pattern = os.path.join(mosdepth_dir, f"**/{sample_id}*.mosdepth.summary.txt")
        matching_files = glob.glob(pattern, recursive=True)
        
        if not matching_files:
            print(f"Warning: No mosdepth summary file found for {sample_id}", file=sys.stderr)
            continue
            
        if len(matching_files) > 1:
            print(f"Warning: Multiple mosdepth files found for {sample_id}, using first one", file=sys.stderr)
        
        # Get coverage from the first matching file
        coverage = get_coverage_from_summary(matching_files[0])
        if coverage is not None:
            results.append((sample_id, coverage))
        else:
            print(f"Warning: Could not extract coverage for {sample_id}", file=sys.stderr)

    # Write results to output file
    with open(output_file, 'w') as f:
        for sample_id, coverage in results:
            f.write(f"{sample_id}\t{coverage:.3f}\n")

    print(f"Processed {len(results)} samples. Results written to {output_file}")

if __name__ == "__main__":
    main()
