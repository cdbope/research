#!/bin/bash
###############################################################################
# Step 1: Normalize VCF files with bcftools
#
# This script normalizes all VCF files to ensure consistent representation
# of variants across samples.
#
# Usage: ./01_normalize_vcfs.sh
###############################################################################

set -e  # Exit on error

# Configuration
REF_FASTA="/home/chbope/extension/nWGS_manuscript_data/data/reference/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa"
INPUT_VCF_DIR="/home/chbope/extension/script/svanalysis/vcf"
OUTPUT_VCF_DIR="/home/chbope/extension/script/svmeta/results/normalized_vcfs"

# Create output directory
mkdir -p "$OUTPUT_VCF_DIR"

echo "============================================================================"
echo "STEP 1: VCF NORMALIZATION"
echo "============================================================================"
echo "Reference: $REF_FASTA"
echo "Input directory: $INPUT_VCF_DIR"
echo "Output directory: $OUTPUT_VCF_DIR"
echo "============================================================================"
echo

# Process each VCF file
for vcf in "$INPUT_VCF_DIR"/*.vcf.gz; do
    filename=$(basename "$vcf")
    sample=$(basename "$vcf" .vcf.gz)
    output="$OUTPUT_VCF_DIR/${sample}.norm.vcf.gz"

    echo "Processing: $filename"

    # Normalize: left-align indels, sort, compress
    bcftools norm -f "$REF_FASTA" "$vcf" 2>/dev/null | \
        bcftools sort -Oz -o "$output"

    # Index
    bcftools index "$output"

    echo "  ✓ Normalized: $output"
done

echo
echo "============================================================================"
echo "VCF normalization complete!"
echo "Normalized VCFs: $OUTPUT_VCF_DIR"
echo "============================================================================"
