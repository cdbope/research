#!/bin/bash
###############################################################################
# Conservative Counting Pipeline - Master Runner
#
# This script runs the complete conservative counting analysis:
# - Step 03: Build gene matrix with conservative counting
# - Step 04: External dataset comparison
#
# ALL outputs go to: /home/chbope/extension/script/svmetaunique/results/
# NO files are modified in svmeta/
#
# Prerequisites:
# - Steps 01-02 must be completed in svmeta/ folder first
# - Merged VCF must exist: svmeta/results/merged/merged_SV.vcf.gz
###############################################################################

set -euo pipefail

SCRIPT_DIR="/home/chbope/extension/script/svmetaunique"
OUTPUT_DIR="${SCRIPT_DIR}/results"

echo "============================================================================"
echo "CONSERVATIVE COUNTING PIPELINE"
echo "============================================================================"
echo "Following npae082.pdf methodology:"
echo "Max 1 SV per gene per sample (binary: 0 or 1)"
echo "============================================================================"
echo ""

# Check prerequisites
echo "Checking prerequisites..."

MERGED_VCF="/home/chbope/extension/script/svmeta/results/merged/merged_SV.vcf.gz"
if [ ! -f "$MERGED_VCF" ]; then
    echo "⚠️  ERROR: Merged VCF not found: $MERGED_VCF"
    echo "Please run steps 01-02 in svmeta/ folder first:"
    echo "  cd /home/chbope/extension/script/svmeta"
    echo "  bash 01_prepare_vcfs.sh"
    echo "  bash 02_merge_with_survivor.sh"
    exit 1
fi

echo "✓ Merged VCF found: $MERGED_VCF"

GENE_BED="/home/chbope/extension/script/svmeta/external_datasets/refseq_genes_hg38.bed"
if [ ! -f "$GENE_BED" ]; then
    echo "⚠️  WARNING: Gene BED not found: $GENE_BED"
    echo "Will attempt to create from GFF3..."
fi

echo ""
echo "============================================================================"
echo "STEP 03: BUILD GENE MATRIX (CONSERVATIVE COUNTING)"
echo "============================================================================"
echo ""

cd "$SCRIPT_DIR"
python3 03_build_matrix_conservative.py

if [ $? -ne 0 ]; then
    echo "⚠️  ERROR: Step 03 failed"
    exit 1
fi

echo ""
echo "============================================================================"
echo "STEP 04: EXTERNAL DATASET COMPARISON"
echo "============================================================================"
echo ""

python3 04_external_dataset_comparison_conservative.py

if [ $? -ne 0 ]; then
    echo "⚠️  ERROR: Step 04 failed"
    exit 1
fi

echo ""
echo "============================================================================"
echo "PIPELINE COMPLETE!"
echo "============================================================================"
echo "All results saved to: $OUTPUT_DIR/"
echo ""
echo "Key output files:"
echo "  • ${OUTPUT_DIR}/matrices/gene_sample_matrix_conservative.csv"
echo "  • ${OUTPUT_DIR}/genes/gene_frequencies_conservative.csv"
echo "  • ${OUTPUT_DIR}/external_comparison/high_confidence_validated_genes_conservative.csv"
echo "  • ${OUTPUT_DIR}/external_comparison/COMPARISON_REPORT_CONSERVATIVE.md"
echo ""
echo "Next steps:"
echo "  1. Review the comparison report"
echo "  2. Compare with standard counting results (optional):"
echo "     python3 05_conservative_gene_counting.py"
echo "============================================================================"
