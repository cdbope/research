import csv
import sys
from collections import defaultdict

def filter_breakpoints(input_file, output_file, paired_output_file=None):
    with open(input_file, newline='') as fin:
        reader = csv.DictReader(fin, delimiter='\t')
        rows = list(reader)
        header = reader.fieldnames

    # Group by (ID, Genes)
    groups = defaultdict(list)
    for row in rows:
        key = (row['ID'], row['Genes'])
        groups[key].append(row)

    # Find keys with both start and end
    to_remove = set()
    for key, group in groups.items():
        breakings = set(row['breaking'] for row in group)
        if 'start' in breakings and 'end' in breakings:
            to_remove.add(key)

    # Write output (rows for (ID, Genes) without both start and end)
    with open(output_file, 'w', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter='\t')
        writer.writeheader()
        for row in rows:
            key = (row['ID'], row['Genes'])
            if key not in to_remove:
                writer.writerow(row)

    # If paired_output_file is provided, filter the output file further
    if paired_output_file:
        with open(output_file, newline='') as fin:
            reader = csv.DictReader(fin, delimiter='\t')
            out_rows = list(reader)
            out_header = reader.fieldnames
        # Group by ID
        id_groups = defaultdict(list)
        for row in out_rows:
            id_groups[row['ID']].append(row)
        # Find IDs with both start and end
        paired_ids = set()
        for key, group in id_groups.items():
            breakings = set(row['breaking'] for row in group)
            if 'start' in breakings and 'end' in breakings:
                paired_ids.add(key)
        # Write only unique rows for IDs with both start and end
        seen = set()
        with open(paired_output_file, 'w', newline='') as fout2:
            writer2 = csv.DictWriter(fout2, fieldnames=out_header, delimiter='\t')
            writer2.writeheader()
            for row in out_rows:
                if row['ID'] in paired_ids:
                    row_key = (row['chr'], row['star'], row['end'], row['ID'], row['svtype'], row['breaking'], row['Genes'])
                    if row_key not in seen:
                        writer2.writerow(row)
                        seen.add(row_key)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Filter breakpoints: remove pairs with both start and end for the same (ID, Genes). Optionally, further filter output to keep only IDs with both start and end, and remove duplicate breakpoints.")
    parser.add_argument('--in', dest='input_file', required=True, help='Input file')
    parser.add_argument('--out', dest='output_file', required=True, help='Output file')
    parser.add_argument('--paired', dest='paired_output_file', required=False, help='Paired output file (optional)')
    args = parser.parse_args()
    filter_breakpoints(args.input_file, args.output_file, args.paired_output_file)
