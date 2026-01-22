#!/bin/bash

################################################################################
# ONT Individual BAM CNV Processing Script using CNVpytor
#
# Purpose: Process a single BAM file for CNV analysis using CNVpytor
#          Optimized for ONT long-read sequencing data
#
# Usage: ./ont_process_individual_bam_cnv.sh -b <bam> -o <output> -r <reference> -s <sample_id> [OPTIONS]
################################################################################

set -euo pipefail

# Default parameters
BAM_FILE=""
OUTPUT_DIR=""
REFERENCE=""
SAMPLE_ID=""
BIN_SIZES="1000 10000 100000"  # Multiple bin sizes for CNVpytor
CONDA_ENV="cnvgen"
THREADS=4

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

usage() {
    cat << EOF
Usage: $0 -b <bam> -o <output> -r <reference> -s <sample_id> [OPTIONS]

Required arguments:
    -b, --bam FILE               Input BAM file
    -o, --output-dir DIR         Output directory
    -r, --reference FILE         Reference genome FASTA file
    -s, --sample-id ID           Sample identifier

Optional arguments:
    --bin-sizes "SIZE1 SIZE2"    Bin sizes for CNVpytor (default: "1000 10000 100000")
    -e, --conda-env ENV          Conda environment (default: cnvgen)
    -t, --threads INT            Number of threads (default: 4)
    -h, --help                   Display this help message

Example:
    $0 -b sample.bam -o ./output -r genome.fa -s SAMPLE001
EOF
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bam)
            BAM_FILE="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -s|--sample-id)
            SAMPLE_ID="$2"
            shift 2
            ;;
        --bin-sizes)
            BIN_SIZES="$2"
            shift 2
            ;;
        -e|--conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$BAM_FILE" ]] || [[ -z "$OUTPUT_DIR" ]] || [[ -z "$REFERENCE" ]] || [[ -z "$SAMPLE_ID" ]]; then
    log_error "Missing required arguments"
    usage
fi

# Validate files exist
if [[ ! -f "$BAM_FILE" ]]; then
    log_error "BAM file does not exist: $BAM_FILE"
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    log_error "Reference file does not exist: $REFERENCE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Activate conda environment
# Temporarily disable unbound variable check for conda activation
set +u
source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"
set -u

log_info "Starting CNV analysis for: $BAM_FILE"
log_info "Sample ID: $SAMPLE_ID"
log_info "Output directory: $OUTPUT_DIR"
log_info "Bin sizes: $BIN_SIZES"

# CNVpytor analysis
PYTOR_FILE="${OUTPUT_DIR}/${SAMPLE_ID}.pytor"

# Check if BAM has reference sequences in header, if not fix it
log_info "Checking BAM header for reference sequences..."
if ! samtools view -H "$BAM_FILE" | grep -q "^@SQ"; then
    log_warn "BAM file missing @SQ lines, adding reference sequences to header..."
    FIXED_BAM="${OUTPUT_DIR}/$(basename ${BAM_FILE%.bam})_fixed.bam"

    # Create new header with @SQ lines
    TMP_HEADER="${OUTPUT_DIR}/tmp_header.sam"
    samtools view -H "$BAM_FILE" | head -1 > "$TMP_HEADER"
    samtools dict "$REFERENCE" | grep "^@SQ" >> "$TMP_HEADER"
    samtools view -H "$BAM_FILE" | tail -n +2 >> "$TMP_HEADER"

    # Reheader BAM
    samtools reheader "$TMP_HEADER" "$BAM_FILE" > "$FIXED_BAM"
    rm "$TMP_HEADER"

    # Index the fixed BAM
    log_info "Indexing fixed BAM file..."
    samtools index "$FIXED_BAM"

    # Use fixed BAM
    BAM_FILE="$FIXED_BAM"
    log_success "BAM header fixed and indexed successfully"
fi

# Step 1: Import read depth from BAM
log_info "Step 1: Importing read depth from BAM file..."
if cnvpytor -root "$PYTOR_FILE" -rd "$BAM_FILE" -T "$REFERENCE"; then
    log_success "Read depth imported successfully"
else
    log_error "Failed to import read depth"
    exit 1
fi

# Step 2: Generate histograms for each bin size
log_info "Step 2: Generating histograms for bin sizes: $BIN_SIZES"
if cnvpytor -root "$PYTOR_FILE" -his $BIN_SIZES; then
    log_success "Histograms generated successfully"
else
    log_error "Failed to generate histograms"
    exit 1
fi

# Step 3: Partition genome
log_info "Step 3: Partitioning genome for bin sizes: $BIN_SIZES"
if cnvpytor -root "$PYTOR_FILE" -partition $BIN_SIZES; then
    log_success "Genome partitioned successfully"
else
    log_error "Failed to partition genome"
    exit 1
fi

# Step 4: Call CNVs for each bin size
log_info "Step 4: Calling CNVs for bin sizes: $BIN_SIZES"
if cnvpytor -root "$PYTOR_FILE" -call $BIN_SIZES; then
    log_success "CNVs called successfully"
else
    log_error "Failed to call CNVs"
    exit 1
fi

# Step 5: Export results in multiple formats
log_info "Step 5: Exporting results..."

# Export VCF format
for binsize in $BIN_SIZES; do
    VCF_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_${binsize}bp.vcf"
    log_info "Exporting VCF for bin size ${binsize}bp: $VCF_FILE"

    if cnvpytor -root "$PYTOR_FILE" -view $binsize > "$VCF_FILE"; then
        log_success "VCF exported: $VCF_FILE"
    else
        log_warn "Failed to export VCF for bin size $binsize"
    fi
done

# Export TSV format with CNV calls
TSV_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_cnv_calls.tsv"
log_info "Exporting TSV format: $TSV_FILE"

# Extract CNV calls to TSV
for binsize in $BIN_SIZES; do
    echo "## Bin size: ${binsize}bp" >> "$TSV_FILE"
    cnvpytor -root "$PYTOR_FILE" -view $binsize | grep -v "^#" >> "$TSV_FILE" 2>/dev/null || true
done

log_success "TSV exported: $TSV_FILE"

# Step 6: Generate plots
log_info "Step 6: Generating CNV plots..."

# Generate Manhattan plot for each bin size
for binsize in $BIN_SIZES; do
    PLOT_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_manhattan_${binsize}bp.png"
    log_info "Generating Manhattan plot for bin size ${binsize}bp: $PLOT_FILE"

    # Use cnvpytor plot function
    cnvpytor -root "$PYTOR_FILE" -plot manhattan $binsize -o "$PLOT_FILE" 2>/dev/null || \
        log_warn "Failed to generate Manhattan plot for bin size $binsize"
done

# Generate circular plot (if available)
CIRCULAR_PLOT="${OUTPUT_DIR}/${SAMPLE_ID}_circular.png"
cnvpytor -root "$PYTOR_FILE" -plot circular 100000 -o "$CIRCULAR_PLOT" 2>/dev/null || \
    log_warn "Circular plot not generated"

# Step 7: Generate summary statistics
log_info "Step 7: Generating summary statistics..."

SUMMARY_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_summary.txt"

# Temporarily disable exit on error for summary generation
set +e

{
    echo "========================================"
    echo "CNV Analysis Summary"
    echo "========================================"
    echo "Sample ID: $SAMPLE_ID"
    echo "BAM File: $BAM_FILE"
    echo "Analysis Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Bin Sizes: $BIN_SIZES"
    echo "========================================"
    echo ""

    for binsize in $BIN_SIZES; do
        echo "CNV Calls for bin size ${binsize}bp:"
        cnvpytor -root "$PYTOR_FILE" -view $binsize 2>/dev/null | grep -v "^#" | wc -l | \
            awk '{print "  Total CNV events: " $1}'

        cnvpytor -root "$PYTOR_FILE" -view $binsize 2>/dev/null | grep -v "^#" | \
            awk '$2 == "duplication"' | wc -l | \
            awk '{print "  Duplications: " $1}'

        cnvpytor -root "$PYTOR_FILE" -view $binsize 2>/dev/null | grep -v "^#" | \
            awk '$2 == "deletion"' | wc -l | \
            awk '{print "  Deletions: " $1}'
        echo ""
    done
} > "$SUMMARY_FILE"

# Re-enable exit on error
set -e

log_success "Summary file created: $SUMMARY_FILE"

# Display summary
cat "$SUMMARY_FILE"

# Step 8: Create metadata file with timestamp
METADATA_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_metadata.json"

{
    echo "{"
    echo "  \"sample_id\": \"$SAMPLE_ID\","
    echo "  \"bam_file\": \"$BAM_FILE\","
    echo "  \"reference\": \"$REFERENCE\","
    echo "  \"timestamp\": \"$(date -Iseconds)\","
    echo "  \"bin_sizes\": [$(echo $BIN_SIZES | sed 's/ /, /g')],"
    echo "  \"pytor_file\": \"$PYTOR_FILE\","
    echo "  \"output_dir\": \"$OUTPUT_DIR\","
    echo "  \"vcf_files\": ["
    first=true
    for binsize in $BIN_SIZES; do
        if [ "$first" = true ]; then
            first=false
        else
            echo ","
        fi
        echo -n "    \"${OUTPUT_DIR}/${SAMPLE_ID}_${binsize}bp.vcf\""
    done
    echo ""
    echo "  ]"
    echo "}"
} > "$METADATA_FILE"

log_success "Metadata file created: $METADATA_FILE"

log_success "CNV analysis completed successfully!"
log_info "Results saved to: $OUTPUT_DIR"

exit 0
