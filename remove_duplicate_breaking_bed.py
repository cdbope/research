import csv
from collections import defaultdict

def read_gene_list(gene_file):
    """Read gene list from a one-column file."""
    genes = set()
    with open(gene_file, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene:  # Skip empty lines
                genes.add(gene)
    return genes

def filter_breakpoints(input_file, output_file, paired_output_file=None, gene_file=None, filtered_output_file=None):
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
        # Group by ID and breaking
        id_breaking_groups = defaultdict(lambda: {'start': [], 'end': []})
        for row in out_rows:
            id_breaking_groups[row['ID']][row['breaking']].append(row)
        # For each ID, pair as many start and end as possible
        paired_rows = []
        for id_, breaks in id_breaking_groups.items():
            n_start = len(breaks['start'])
            n_end = len(breaks['end'])
            n_pair = min(n_start, n_end)
            # Keep only up to n_pair of each
            paired_rows.extend(breaks['start'][:n_pair])
            paired_rows.extend(breaks['end'][:n_pair])
        # Remove duplicates
        seen = set()
        with open(paired_output_file, 'w', newline='') as fout2:
            writer2 = csv.DictWriter(fout2, fieldnames=out_header, delimiter='\t')
            writer2.writeheader()
            for row in paired_rows:
                row_key = (row['chr'], row['star'], row['end'], row['ID'], row['svtype'], row['breaking'], row['Genes'])
                if row_key not in seen:
                    writer2.writerow(row)
                    seen.add(row_key)
        
        # If gene_file is provided, filter paired output by genes
        if gene_file and filtered_output_file:
            # Read gene list
            gene_list = read_gene_list(gene_file)
            print(f"Loaded {len(gene_list)} genes from {gene_file}")
            
            # Filter paired output to keep only rows with genes in the gene list
            with open(paired_output_file, newline='') as fin:
                reader = csv.DictReader(fin, delimiter='\t')
                paired_rows = list(reader)
                paired_header = reader.fieldnames
            
            # Filter rows where Genes column contains a gene from the gene list AND SVTYPE is not "INV"
            filtered_rows = []
            for row in paired_rows:
                gene = row['Genes']
                svtype = row['svtype']
                if gene in gene_list and svtype != "INV":
                    filtered_rows.append(row)
            
            # Write filtered results
            with open(filtered_output_file, 'w', newline='') as fout3:
                writer3 = csv.DictWriter(fout3, fieldnames=paired_header, delimiter='\t')
                writer3.writeheader()
                for row in filtered_rows:
                    writer3.writerow(row)
            
            print(f"Kept {len(filtered_rows)} rows with genes from the gene list (excluding INV variants)")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Filter breakpoints: remove pairs with both start and end for the same (ID, Genes). Optionally, further filter output to keep only as many start and end as can be paired, and remove duplicate breakpoints. Can also filter by gene list.")
    parser.add_argument('--in', dest='input_file', required=True, help='Input file')
    parser.add_argument('--out', dest='output_file', required=True, help='Output file')
    parser.add_argument('--paired', dest='paired_output_file', required=False, help='Paired output file (optional)')
    parser.add_argument('--gene-list', dest='gene_file', required=False, help='Gene list file (one gene per line)')
    parser.add_argument('--filtered', dest='filtered_output_file', required=False, help='Filtered output file (genes from gene list only)')
    args = parser.parse_args()
    filter_breakpoints(args.input_file, args.output_file, args.paired_output_file, args.gene_file, args.filtered_output_file)