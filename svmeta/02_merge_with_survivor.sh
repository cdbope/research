#!/bin/bash
###############################################################################
# Step 2: Merge VCFs across samples using SURVIVOR
#
# After merging:
#  - Sort merged VCF
#  - Compress with bgzip
#  - Index with tabix
#
# Output: merged_SV.vcf.gz + merged_SV.vcf.gz.tbi
###############################################################################

set -euo pipefail

# Configuration
PREPARED_VCF_DIR="/home/chbope/extension/script/svmeta/results/prepared_vcfs/filter_vcf"
OUTPUT_DIR="/home/chbope/extension/script/svmeta/results/merged"
SURVIVOR_BIN="/home/chbope/extension/script/svmeta/tools/SURVIVOR/Debug/SURVIVOR"

# SURVIVOR parameters
MAX_DIST=1000          # Maximum distance between SVs to merge (bp)
MIN_SUPPORT=1          # Minimum number of samples to support an event
TYPE_MATCH=1           # Require same SV type (1=yes)
STRAND_MATCH=0         # Ignore strand (0=no)
USE_SIZE=1             # Use size similarity
SIZE_THRESHOLD=0.5     # Size similarity threshold (30%)

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create VCF list for input
VCF_LIST="$OUTPUT_DIR/vcf_list.txt"
ls "$PREPARED_VCF_DIR"/*.vcf > "$VCF_LIST"

n_vcfs=$(wc -l < "$VCF_LIST")
echo "Found $n_vcfs VCF files for merging."

echo
echo "==========================="
echo "Running SURVIVOR merge..."
echo "==========================="

# Run SURVIVOR merge → produces merged_SV.vcf
$SURVIVOR_BIN merge \
    "$VCF_LIST" \
    "$MAX_DIST" \
    "$MIN_SUPPORT" \
    "$TYPE_MATCH" \
    "$STRAND_MATCH" \
    "$USE_SIZE" \
    "$SIZE_THRESHOLD" \
    "$OUTPUT_DIR/merged_SV.vcf"

echo "✓ SURVIVOR merge complete."
echo

# ---------------------------------------------------------------------------
# Sort, compress, and index merged VCF
# ---------------------------------------------------------------------------

echo "Sorting merged VCF..."
bcftools sort "$OUTPUT_DIR/merged_SV.vcf" -Oz -o "$OUTPUT_DIR/merged_SV.vcf.gz"

echo "Indexing merged VCF with tabix..."
tabix -p vcf "$OUTPUT_DIR/merged_SV.vcf.gz"

echo "✓ Sorting + bgzip + tabix indexing complete."
echo

# Optional: show quick stats
echo "Quick statistics:"
bcftools stats "$OUTPUT_DIR/merged_SV.vcf.gz" | grep '^SN' | grep 'number of'

echo
echo "==========================================================="
echo "SURVIVOR merge pipeline complete!"
echo "Output: $OUTPUT_DIR/merged_SV.vcf.gz"
echo "Index:  $OUTPUT_DIR/merged_SV.vcf.gz.tbi"
echo "==========================================================="
