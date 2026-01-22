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

def process_intersectbed_output(input_file):
    """Process intersectBed output to extract gene names and format for filter_breakpoints. Returns processed data."""
    processed_data = []
    
    try:
        with open(input_file, 'r', encoding='utf-8') as fin:
            line_count = 0
            for line in fin:
                line_count += 1
                if line.strip() == "":
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    print(f"Warning: Line {line_count} has fewer than 10 fields: {len(fields)}")
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
                        
                        # Create row dictionary
                        row = {
                            'chr': chr_col,
                            'star': start_col,
                            'end': end_col,
                            'ID': id_part,
                            'svtype': svtype_part,
                            'breaking': breaking_part,
                            'Genes': gene
                        }
                        processed_data.append(row)
                    else:
                        print(f"Warning: Line {line_count} has invalid name format: {name_col}")
                else:
                    print(f"Warning: Line {line_count} has no '|' in name field: {name_col}")
        
        print(f"Successfully processed {len(processed_data)} rows from {input_file}")
        return processed_data
    
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")
        return []

def read_gene_list(gene_file):
    """Read gene list from a one-column file."""
    genes = set()
    try:
        with open(gene_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                gene = line.strip()
                if gene:  # Skip empty lines
                    genes.add(gene)
        print(f"Successfully loaded {len(genes)} genes from {gene_file}")
        return genes
    except Exception as e:
        print(f"Error reading gene file {gene_file}: {str(e)}")
        return set()

def filter_breakpoints_from_data(processed_data, output_file, paired_output_file=None, gene_file=None, filtered_output_file=None):
    """Filter breakpoints function that works with data in memory instead of reading from file."""
    if not processed_data:
        print("Warning: No data to process!")
        return
    
    rows = processed_data
    header = ['chr', 'star', 'end', 'ID', 'svtype', 'breaking', 'Genes']
    
    print(f"Starting filtering with {len(rows)} input rows")

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

    print(f"Found {len(to_remove)} (ID, Genes) pairs with both start and end")

    # Write output (rows for (ID, Genes) without both start and end)
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as fout:
            writer = csv.DictWriter(fout, fieldnames=header, delimiter='\t')
            writer.writeheader()
            written_count = 0
            for row in rows:
                key = (row['ID'], row['Genes'])
                if key not in to_remove:
                    writer.writerow(row)
                    written_count += 1
        print(f"Wrote {written_count} rows to {output_file}")
    except Exception as e:
        print(f"Error writing to {output_file}: {str(e)}")

    # If paired_output_file is provided, filter the output file further
    if paired_output_file:
        # Use the filtered data directly instead of reading from file
        out_rows = [row for row in rows if (row['ID'], row['Genes']) not in to_remove]
        out_header = header
        print(f"Processing {len(out_rows)} rows for pairing")
        
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
        
        print(f"Created {len(paired_rows)} paired rows")
        
        # Remove duplicates
        seen = set()
        try:
            with open(paired_output_file, 'w', newline='', encoding='utf-8') as fout2:
                writer2 = csv.DictWriter(fout2, fieldnames=out_header, delimiter='\t')
                writer2.writeheader()
                unique_count = 0
                for row in paired_rows:
                    row_key = (row['chr'], row['star'], row['end'], row['ID'], row['svtype'], row['breaking'], row['Genes'])
                    if row_key not in seen:
                        writer2.writerow(row)
                        seen.add(row_key)
                        unique_count += 1
            print(f"Wrote {unique_count} unique paired rows to {paired_output_file}")
        except Exception as e:
            print(f"Error writing to {paired_output_file}: {str(e)}")
        
        # If gene_file is provided, filter paired output by genes
        if gene_file and filtered_output_file:
            # Read gene list
            gene_list = read_gene_list(gene_file)
            print(f"Loaded {len(gene_list)} genes from {gene_file}")
            
            # Filter rows where Genes column contains a gene from the gene list AND SVTYPE is not "INV"
            filtered_rows = []
            for row in paired_rows:
                row_key = (row['chr'], row['star'], row['end'], row['ID'], row['svtype'], row['breaking'], row['Genes'])
                if row_key not in seen:
                    gene = row['Genes']
                    svtype = row['svtype']
                    if gene in gene_list and svtype != "INV":
                        filtered_rows.append(row)
                    seen.add(row_key)
            
            # Write filtered results
            try:
                with open(filtered_output_file, 'w', newline='', encoding='utf-8') as fout3:
                    writer3 = csv.DictWriter(fout3, fieldnames=out_header, delimiter='\t')
                    writer3.writeheader()
                    for row in filtered_rows:
                        writer3.writerow(row)
                print(f"Wrote {len(filtered_rows)} filtered rows to {filtered_output_file}")
            except Exception as e:
                print(f"Error writing to {filtered_output_file}: {str(e)}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process intersectBed output and filter breakpoints: extract gene names, format data, and apply filtering.")
    parser.add_argument('--intersectbed', dest='intersectbed_file', required=True, help='Input intersectBed output file')
    parser.add_argument('--out', dest='output_file', required=True, help='Filtered output file')
    parser.add_argument('--paired', dest='paired_output_file', required=False, help='Paired output file (optional)')
    parser.add_argument('--gene-list', dest='gene_file', required=False, help='Gene list file (one gene per line)')
    parser.add_argument('--filtered', dest='filtered_output_file', required=False, help='Filtered output file (genes from gene list only)')
    args = parser.parse_args()
    
    print(f"Python version: {sys.version}")
    print(f"Input file: {args.intersectbed_file}")
    
    # Step 1: Process intersectBed output to extract gene names and format data
    print("Processing intersectBed output...")
    processed_data = process_intersectbed_output(args.intersectbed_file)
    print(f"Processed {len(processed_data)} rows from intersectBed output")
    
    if not processed_data:
        print("ERROR: No data was processed. Check your input file format.")
        sys.exit(1)
    
    # Step 2: Apply filtering
    print("Applying breakpoint filtering...")
    filter_breakpoints_from_data(processed_data, args.output_file, args.paired_output_file, args.gene_file, args.filtered_output_file)
    print("Processing complete!") 
