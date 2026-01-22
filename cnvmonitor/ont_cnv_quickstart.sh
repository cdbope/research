#!/bin/bash

################################################################################
# ONT Time-Series CNV Analysis - Quick Start Script
#
# This script helps you quickly set up and run the CNV analysis pipeline
################################################################################

set -euo pipefail

# Colors
BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  ONT Time-Series CNV Analysis Pipeline - Quick Start      ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

show_menu() {
    echo -e "${GREEN}What would you like to do?${NC}"
    echo ""
    echo "1) Setup - Check/install dependencies"
    echo "2) Test - Run a test with example data"
    echo "3) Run - Start continuous monitoring"
    echo "4) Process - Process existing BAM files once"
    echo "5) View - Check analysis results"
    echo "6) Help - Show detailed documentation"
    echo "7) Exit"
    echo ""
}

check_dependencies() {
    echo -e "${BLUE}[INFO] Checking dependencies...${NC}"
    echo ""

    # Check conda
    if command -v conda &> /dev/null; then
        echo -e "${GREEN}✓${NC} Conda is installed"
    else
        echo -e "${RED}✗${NC} Conda is not installed"
        echo "Please install Miniconda or Anaconda first"
        exit 1
    fi

    # Check if cnvgen environment exists
    if conda env list | grep -q "^cnvgen "; then
        echo -e "${GREEN}✓${NC} cnvgen environment exists"
    else
        echo -e "${YELLOW}⚠${NC} cnvgen environment not found"
        echo ""
        echo "Creating cnvgen environment..."
        if [[ -f "${SCRIPT_DIR}/cnvgen_environment.yml" ]]; then
            conda env create -f "${SCRIPT_DIR}/cnvgen_environment.yml"
            echo -e "${GREEN}✓${NC} cnvgen environment created"
        else
            echo -e "${RED}✗${NC} cnvgen_environment.yml not found"
            exit 1
        fi
    fi

    # Activate and check key tools
    echo ""
    echo "Checking tools in cnvgen environment..."

    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate cnvgen

    if command -v cnvpytor &> /dev/null; then
        echo -e "${GREEN}✓${NC} CNVpytor: $(cnvpytor --version 2>&1 | head -1)"
    else
        echo -e "${RED}✗${NC} CNVpytor not found"
    fi

    if command -v samtools &> /dev/null; then
        echo -e "${GREEN}✓${NC} Samtools: $(samtools --version | head -1)"
    else
        echo -e "${RED}✗${NC} Samtools not found"
    fi

    if Rscript -e "library(ggplot2)" 2>&1 | grep -q "Error"; then
        echo -e "${RED}✗${NC} R packages not properly installed"
    else
        echo -e "${GREEN}✓${NC} R packages installed"
    fi

    echo ""
    echo -e "${GREEN}✓ Dependencies check complete!${NC}"
}

show_example_commands() {
    echo -e "${BLUE}[INFO] Example Commands${NC}"
    echo ""

    echo -e "${YELLOW}Continuous monitoring (real-time analysis):${NC}"
    echo "./ont_timeseries_cnv_monitor.sh \\"
    echo "  -i /path/to/bam_directory \\"
    echo "  -o /path/to/output \\"
    echo "  -r /path/to/reference.fa \\"
    echo "  -s SAMPLE_ID \\"
    echo "  -c 300 \\"
    echo "  -t 4"
    echo ""

    echo -e "${YELLOW}One-time processing (existing BAMs):${NC}"
    echo "./ont_timeseries_cnv_monitor.sh \\"
    echo "  -i /path/to/bam_directory \\"
    echo "  -o /path/to/output \\"
    echo "  -r /path/to/reference.fa \\"
    echo "  -s SAMPLE_ID \\"
    echo "  --no-continuous \\"
    echo "  -t 8"
    echo ""

    echo -e "${YELLOW}Process single BAM:${NC}"
    echo "./ont_process_individual_bam_cnv.sh \\"
    echo "  -b sample.bam \\"
    echo "  -o ./output \\"
    echo "  -r reference.fa \\"
    echo "  -s SAMPLE_ID"
    echo ""
}

run_continuous() {
    echo -e "${BLUE}[INFO] Setting up continuous monitoring...${NC}"
    echo ""

    read -p "BAM input directory: " INPUT_DIR
    read -p "Output directory: " OUTPUT_DIR
    read -p "Reference genome FASTA: " REFERENCE
    read -p "Sample ID: " SAMPLE_ID
    read -p "Check interval (seconds) [300]: " CHECK_INTERVAL
    CHECK_INTERVAL=${CHECK_INTERVAL:-300}
    read -p "Number of threads [4]: " THREADS
    THREADS=${THREADS:-4}

    echo ""
    echo -e "${YELLOW}Starting continuous monitoring...${NC}"
    echo "Press Ctrl+C to stop"
    echo ""

    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate cnvgen

    "${SCRIPT_DIR}/ont_timeseries_cnv_monitor.sh" \
        -i "$INPUT_DIR" \
        -o "$OUTPUT_DIR" \
        -r "$REFERENCE" \
        -s "$SAMPLE_ID" \
        -c "$CHECK_INTERVAL" \
        -t "$THREADS"
}

run_batch() {
    echo -e "${BLUE}[INFO] Setting up batch processing...${NC}"
    echo ""

    read -p "BAM input directory: " INPUT_DIR
    read -p "Output directory: " OUTPUT_DIR
    read -p "Reference genome FASTA: " REFERENCE
    read -p "Sample ID: " SAMPLE_ID
    read -p "Number of threads [8]: " THREADS
    THREADS=${THREADS:-8}

    echo ""
    echo -e "${YELLOW}Starting batch processing...${NC}"
    echo ""

    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate cnvgen

    "${SCRIPT_DIR}/ont_timeseries_cnv_monitor.sh" \
        -i "$INPUT_DIR" \
        -o "$OUTPUT_DIR" \
        -r "$REFERENCE" \
        -s "$SAMPLE_ID" \
        --no-continuous \
        -t "$THREADS"

    echo ""
    echo -e "${GREEN}✓ Batch processing complete!${NC}"
}

view_results() {
    echo -e "${BLUE}[INFO] View analysis results${NC}"
    echo ""

    read -p "Output directory: " OUTPUT_DIR

    if [[ ! -d "$OUTPUT_DIR" ]]; then
        echo -e "${RED}✗${NC} Directory not found: $OUTPUT_DIR"
        return
    fi

    echo ""
    echo -e "${GREEN}Directory Structure:${NC}"
    tree -L 2 "$OUTPUT_DIR" 2>/dev/null || ls -R "$OUTPUT_DIR"

    echo ""
    echo -e "${GREEN}Processed BAMs:${NC}"
    if [[ -f "${OUTPUT_DIR}/logs/"*"_processed_bams.txt" ]]; then
        wc -l "${OUTPUT_DIR}/logs/"*"_processed_bams.txt"
    else
        echo "No processed BAMs found"
    fi

    echo ""
    echo -e "${GREEN}Latest Processing Log:${NC}"
    if [[ -f "${OUTPUT_DIR}/logs/"*"_processing_log.txt" ]]; then
        tail -20 "${OUTPUT_DIR}/logs/"*"_processing_log.txt"
    else
        echo "No processing log found"
    fi

    echo ""
    echo -e "${GREEN}Summary Reports:${NC}"
    find "$OUTPUT_DIR" -name "*summary*.txt" -o -name "*summary*.csv"

    echo ""
    echo -e "${GREEN}Visualizations:${NC}"
    find "$OUTPUT_DIR" -name "*.png" -o -name "*.pdf"
}

show_help() {
    if [[ -f "${SCRIPT_DIR}/ONT_TIMESERIES_CNV_README.md" ]]; then
        less "${SCRIPT_DIR}/ONT_TIMESERIES_CNV_README.md"
    else
        echo -e "${RED}✗${NC} README.md not found"
    fi
}

# Main menu loop
while true; do
    show_menu
    read -p "Select option [1-7]: " choice

    case $choice in
        1)
            check_dependencies
            echo ""
            read -p "Press Enter to continue..."
            ;;
        2)
            echo -e "${YELLOW}⚠ Test mode not yet implemented${NC}"
            echo "Please use option 4 to process your own BAM files"
            read -p "Press Enter to continue..."
            ;;
        3)
            run_continuous
            ;;
        4)
            run_batch
            echo ""
            read -p "Press Enter to continue..."
            ;;
        5)
            view_results
            echo ""
            read -p "Press Enter to continue..."
            ;;
        6)
            show_help
            ;;
        7)
            echo -e "${GREEN}Goodbye!${NC}"
            exit 0
            ;;
        *)
            echo -e "${RED}Invalid option${NC}"
            read -p "Press Enter to continue..."
            ;;
    esac

    clear
done
