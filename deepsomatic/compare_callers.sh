#!/bin/bash

# Variant Caller Comparison Wrapper Script
# Compare DeepSomatic vs Clair3/ClairS-TO variant calls

# Default paths
DEEPSOMATIC_FILE="/home/chbope/extension/script/deepsomatic/output/T25-152_annotateandfilter_deep_somatic.csv"
CLAIR_FILE="/home/chbope/extension/data/t25-152/merge_annot_clair3andclairsto/T25-152_merge_annotation_filter_snvs_allcall.csv"
OUTPUT_DIR="/home/chbope/extension/script/deepsomatic/comparison_results"

# Usage information
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Compare variant calls between DeepSomatic and Clair3/ClairS-TO

OPTIONS:
    -d FILE     DeepSomatic annotated CSV file
    -c FILE     Clair3/ClairS-TO annotated CSV file
    -o DIR      Output directory (default: comparison_results)
    -h          Show this help message

EXAMPLES:
    # Use default paths
    $0

    # Specify custom files
    $0 -d deepsomatic.csv -c clair.csv -o output_dir

    # For different sample
    $0 -d T001_annotateandfilter_deep_somatic.csv -c T001_clair.csv

EOF
}

# Parse command line arguments
while getopts "d:c:o:h" opt; do
    case $opt in
        d) DEEPSOMATIC_FILE="$OPTARG" ;;
        c) CLAIR_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

echo "=========================================="
echo "Variant Caller Comparison Tool"
echo "=========================================="
echo ""
echo "Input files:"
echo "  DeepSomatic: $DEEPSOMATIC_FILE"
echo "  Clair3/ClairS-TO: $CLAIR_FILE"
echo "  Output: $OUTPUT_DIR"
echo ""

# Check if files exist
if [ ! -f "$DEEPSOMATIC_FILE" ]; then
    echo "Error: DeepSomatic file not found: $DEEPSOMATIC_FILE"
    exit 1
fi

if [ ! -f "$CLAIR_FILE" ]; then
    echo "Error: Clair3/ClairS-TO file not found: $CLAIR_FILE"
    exit 1
fi

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 not found. Please install Python 3."
    exit 1
fi

# Check if pandas is available
if ! python3 -c "import pandas" 2>/dev/null; then
    echo "Error: pandas not installed. Installing..."
    pip3 install pandas
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run Python comparison script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$SCRIPT_DIR/compare_callers.py" \
    "$DEEPSOMATIC_FILE" \
    "$CLAIR_FILE" \
    "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
    echo ""
    echo "Comparison complete! Check results in: $OUTPUT_DIR"
    echo ""
    echo "Generated files:"
    ls -lh "$OUTPUT_DIR"
else
    echo "Error: Comparison failed"
    exit 1
fi
