#!/bin/bash
#
# Structural Variant Analysis Runner
# Quick script to run the complete SV analysis pipeline
#

echo "========================================"
echo "SV Analysis for GBM Samples"
echo "========================================"
echo ""

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo "Error: Python not found. Please install Python 3.x"
    exit 1
fi

# Check if required packages are installed
echo "Checking dependencies..."
python -c "import pandas, numpy, matplotlib, seaborn, scipy, sklearn" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip install -r requirements.txt
    if [ $? -ne 0 ]; then
        echo "Error: Failed to install dependencies"
        exit 1
    fi
fi

echo "Dependencies OK"
echo ""

# Run analysis
echo "Running SV analysis pipeline..."
echo ""
python sv_analysis.py

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "Analysis completed successfully!"
    echo "========================================"
    echo ""
    echo "Results saved to: ./results/"
    echo ""
    echo "View key files:"
    echo "  - results/sv_summary_display.csv (main table)"
    echo "  - results/clustering_dendrogram.png"
    echo "  - results/pca_clustering.png"
    echo "  - results/vaf_analysis.png"
    echo "  - results/analysis_report.txt"
    echo ""
else
    echo ""
    echo "Error: Analysis failed. Please check error messages above."
    exit 1
fi
