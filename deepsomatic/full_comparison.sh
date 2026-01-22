#!/bin/bash

# Complete Variant Caller Comparison with Recommendation
# Runs both comparison and evaluation

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default paths
DEEPSOMATIC_FILE="/home/chbope/extension/script/deepsomatic/output/T25-152_annotateandfilter_deep_somatic.csv"
CLAIR_FILE="/home/chbope/extension/data/t25-152/merge_annot_clair3andclairsto/T25-152_merge_annotation_filter_snvs_allcall.csv"
OUTPUT_DIR="/home/chbope/extension/script/deepsomatic/comparison_results"

# Usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Complete variant caller comparison and recommendation

OPTIONS:
    -d FILE     DeepSomatic annotated CSV file
    -c FILE     Clair3/ClairS-TO annotated CSV file
    -o DIR      Output directory
    -h          Show this help message

EXAMPLES:
    # Use default paths (T25-152)
    $0

    # Custom files
    $0 -d deepsomatic.csv -c clair.csv -o output

EOF
}

# Parse arguments
while getopts "d:c:o:h" opt; do
    case $opt in
        d) DEEPSOMATIC_FILE="$OPTARG" ;;
        c) CLAIR_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

echo "======================================================================"
echo "COMPLETE VARIANT CALLER COMPARISON AND RECOMMENDATION"
echo "======================================================================"
echo ""
echo "Input files:"
echo "  DeepSomatic:      $DEEPSOMATIC_FILE"
echo "  Clair3/ClairS-TO: $CLAIR_FILE"
echo "  Output:           $OUTPUT_DIR"
echo ""

# Check files
if [ ! -f "$DEEPSOMATIC_FILE" ]; then
    echo "Error: DeepSomatic file not found: $DEEPSOMATIC_FILE"
    exit 1
fi

if [ ! -f "$CLAIR_FILE" ]; then
    echo "Error: Clair3/ClairS-TO file not found: $CLAIR_FILE"
    exit 1
fi

# Check Python and pandas
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 not found"
    exit 1
fi

if ! python3 -c "import pandas" 2>/dev/null; then
    echo "Installing pandas..."
    pip3 install pandas
fi

mkdir -p "$OUTPUT_DIR"

# Step 1: Run comparison
echo ""
echo "Step 1: Running variant comparison..."
echo "----------------------------------------------------------------------"
python3 "$SCRIPT_DIR/compare_callers.py" \
    "$DEEPSOMATIC_FILE" \
    "$CLAIR_FILE" \
    "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo "Error: Comparison failed"
    exit 1
fi

# Step 2: Run evaluation
echo ""
echo "Step 2: Running caller evaluation and generating recommendation..."
echo "----------------------------------------------------------------------"
python3 "$SCRIPT_DIR/evaluate_callers.py" \
    "$DEEPSOMATIC_FILE" \
    "$CLAIR_FILE" \
    "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo "Error: Evaluation failed"
    exit 1
fi

# Summary
echo ""
echo "======================================================================"
echo "ANALYSIS COMPLETE!"
echo "======================================================================"
echo ""
echo "Generated files in: $OUTPUT_DIR"
echo ""
ls -lh "$OUTPUT_DIR"
echo ""
echo "Key files:"
echo "  - shared_variants.tsv           : High-confidence variants (both callers)"
echo "  - deepsomatic_only.tsv          : DeepSomatic unique variants"
echo "  - clair_only.tsv                : Clair3/ClairS-TO unique variants"
echo "  - comparison_summary.txt        : Quick statistics"
echo "  - caller_recommendation.txt     : Which caller is better"
echo ""
echo "Quick summary:"
echo "----------------------------------------------------------------------"
cat "$OUTPUT_DIR/caller_recommendation.txt" | head -20
echo ""
echo "For full details, see: $OUTPUT_DIR/caller_recommendation.txt"
echo "======================================================================"
