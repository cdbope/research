#!/usr/bin/env python3
"""
Create gene BED file from GENCODE GFF3 annotation
Output format: chr start end gene_name
"""

import gzip
import sys

def extract_gene_name(attributes):
    """Extract gene_name from GFF3 attributes field"""
    for attr in attributes.split(';'):
        if attr.startswith('gene_name='):
            return attr.split('=')[1]
    return None

def main():
    input_gff3 = "/home/chbope/extension/nWGS_manuscript_data/data/reference/gencode.v48.annotation.gff3"
    output_bed = "/home/chbope/extension/script/svmeta/external_datasets/refseq_genes_hg38.bed"

    print(f"Reading GFF3: {input_gff3}")
    print("Extracting gene features...")

    genes = []
    gene_count = 0

    with open(input_gff3, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]

            # Extract gene name
            gene_name = extract_gene_name(attributes)

            if gene_name:
                genes.append((chrom, start, end, gene_name))
                gene_count += 1

                if gene_count % 5000 == 0:
                    print(f"  Processed {gene_count} genes...")

    print(f"✓ Extracted {len(genes)} genes")

    # Write BED file
    print(f"\nWriting BED file: {output_bed}")

    with open(output_bed, 'w') as f:
        for chrom, start, end, gene_name in sorted(genes):
            f.write(f"{chrom}\t{start}\t{end}\t{gene_name}\n")

    print(f"✓ Created BED file with {len(genes)} genes")
    print(f"\nFile saved: {output_bed}")

    # Show first few lines
    print("\nFirst 10 genes:")
    with open(output_bed, 'r') as f:
        for i, line in enumerate(f):
            if i >= 10:
                break
            print(f"  {line.strip()}")

if __name__ == "__main__":
    main()
