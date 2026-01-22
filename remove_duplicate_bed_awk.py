import csv
import sys
from collections import defaultdict

def extract_gene_name(attributes):
    """Extract gene_name from the attributes string (equivalent to the first awk command)."""
    gene = "NA"
    if attributes:
        fields = attributes.split(';')
        for field in fields:
            if field.startswith('gene_name='):
                parts = field.split('=', 1)
                if len(parts) == 2:
                    gene = parts[1]
                break
    return gene

def process_intersectbed_output(input_file, output_file):
    """Process intersectBed output to extract gene names and format for filter_breakpoints."""
    with open(input_file, 'r') as fin, open(output_file, 'w', newline='') as fout:
        writer = csv.writer(fout, delimiter='\t')
        writer.writerow(['chr', 'star', 'end', 'ID', 'svtype', 'breaking', 'Genes'])
        
        for line in fin:
            if line.strip() == "":
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            
            # Extract fields from intersectBed output
            chr_col = fields[0]
            start_col = fields[1]
            end_col = fields[2]
            name_col = fields[3]
            attributes = fields[-1]  # Last column contains attributes
            
            # Extract gene name from attributes (equivalent to first awk)
            gene = extract_gene_name(attributes)
            
            # Parse the name field to extract ID, svtype, breaking (equivalent to second awk)
            if '|' in name_col:
                parts = name_col.split('|')
                if len(parts) >= 3:
                    id_part = parts[0]
                    svtype_part = parts[1]
                    breaking_part = parts[2]
                    
                    # Write formatted output
                    writer.writerow([chr_col, start_col, end_col, id_part, svtype_part, breaking_part, gene])

def read_gene_list_old(gene_file):
    """Read gene list from a one-column file."""
    genes = set()
    with open(gene_file, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene:  # Skip empty lines
                genes.add(gene)
    return genes

def read_gene_list(gene_file):
    genes = set()
    with open(gene_file, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            gene = line.strip().replace('\r', '').upper()
            if i == 0 and gene == "GENES":
                continue  # skip header
            if gene:
                genes.add(gene)
    return genes

def filter_breakpoints(input_file, output_file, paired_output_file=None, gene_file=None, filtered_output_file=None):
    """Filter breakpoints function with improved gene name handling and encoding robustness."""
    import csv
    from collections import defaultdict

    # Step 1: Initial filtering based on start/end pairs
    with open(input_file, newline='', encoding='utf-8') as fin:
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

    # Write filtered output
    with open(output_file, 'w', newline='', encoding='utf-8') as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter='\t')
        writer.writeheader()
        for row in rows:
            key = (row['ID'], row['Genes'])
            if key not in to_remove:
                writer.writerow(row)

    # Step 2: Pair start and end
    if paired_output_file:
        with open(output_file, newline='', encoding='utf-8') as fin:
            reader = csv.DictReader(fin, delimiter='\t')
            out_rows = list(reader)
            out_header = reader.fieldnames

        id_breaking_groups = defaultdict(lambda: {'start': [], 'end': []})
        for row in out_rows:
            id_breaking_groups[row['ID']][row['breaking']].append(row)

        paired_rows = []
        for id_, breaks in id_breaking_groups.items():
            n_pair = min(len(breaks['start']), len(breaks['end']))
            paired_rows.extend(breaks['start'][:n_pair])
            paired_rows.extend(breaks['end'][:n_pair])

        # Remove duplicates
        seen = set()
        with open(paired_output_file, 'w', newline='', encoding='utf-8') as fout2:
            writer2 = csv.DictWriter(fout2, fieldnames=out_header, delimiter='\t')
            writer2.writeheader()
            for row in paired_rows:
                row_key = tuple(row[col] for col in ['chr', 'star', 'end', 'ID', 'svtype', 'breaking', 'Genes'])
                if row_key not in seen:
                    writer2.writerow(row)
                    seen.add(row_key)

        # Step 3: Final filtering by gene list
        if gene_file and filtered_output_file:
            gene_list = set()
            with open(gene_file, 'r', encoding='utf-8') as f:
                for line in f:
                    gene = line.strip().replace('\r', '').upper()
                    if gene:
                        gene_list.add(gene)

            print(f"Loaded {len(gene_list)} genes from {gene_file}")

            with open(paired_output_file, newline='', encoding='utf-8') as fin:
                reader = csv.DictReader(fin, delimiter='\t')
                paired_rows = list(reader)
                paired_header = reader.fieldnames

            filtered_rows = []
            for row in paired_rows:
                gene = row['Genes'].strip().replace('\r', '').upper()
                svtype = row['svtype'].strip().upper()
                if gene in gene_list and svtype != "INV":
                    filtered_rows.append(row)

            with open(filtered_output_file, 'w', newline='', encoding='utf-8') as fout3:
                writer3 = csv.DictWriter(fout3, fieldnames=paired_header, delimiter='\t')
                writer3.writeheader()
                for row in filtered_rows:
                    writer3.writerow(row)

            print(f"Kept {len(filtered_rows)} rows with genes from the gene list (excluding INV variants)")


def filter_breakpoints_old(input_file, output_file, paired_output_file=None, gene_file=None, filtered_output_file=None):
    """Filter breakpoints function (same as in remove_duplicate_breaking_bed.py)."""
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
        with open(output_file, newline='',encoding='utf-8') as fin:
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
    parser = argparse.ArgumentParser(description="Process intersectBed output and filter breakpoints: extract gene names, format data, and apply filtering.")
    parser.add_argument('--intersectbed', dest='intersectbed_file', required=True, help='Input intersectBed output file')
    parser.add_argument('--formatted', dest='formatted_file', required=True, help='Formatted output file (ready for filtering)')
    parser.add_argument('--out', dest='output_file', required=True, help='Filtered output file')
    parser.add_argument('--paired', dest='paired_output_file', required=False, help='Paired output file (optional)')
    parser.add_argument('--gene-list', dest='gene_file', required=False, help='Gene list file (one gene per line)')
    parser.add_argument('--filtered', dest='filtered_output_file', required=False, help='Filtered output file (genes from gene list only)')
    args = parser.parse_args()
    
    # Step 1: Process intersectBed output to extract gene names and format data
    print("Processing intersectBed output...")
    process_intersectbed_output(args.intersectbed_file, args.formatted_file)
    print(f"Formatted data saved to {args.formatted_file}")
    
    # Step 2: Apply filtering
    print("Applying breakpoint filtering...")
    filter_breakpoints(args.formatted_file, args.output_file, args.paired_output_file, args.gene_file, args.filtered_output_file)
    print("Processing complete!") 
