#!/bin/bash

################################################################################
# ONT Time-Series CNV Monitoring Script
#
# Purpose: Monitor a directory for new BAM files generated from ONT basecalling
#          at intervals and process them for CNV analysis using CNVpytor
#
# Usage: ./ont_timeseries_cnv_monitor.sh -i <input_dir> -o <output_dir> -r <reference> -s <sample_id> [OPTIONS]
################################################################################

set -euo pipefail

# Default parameters
INPUT_DIR=""
OUTPUT_DIR=""
REFERENCE=""
SAMPLE_ID=""
CHECK_INTERVAL=300  # Check every 5 minutes
MIN_BAM_SIZE=10000000  # Minimum BAM size in bytes (10MB)
CONDA_ENV="cnvgen"
BIN_SIZES="1000 10000 100000"
THREADS=4
CONTINUOUS=true
MAX_ITERATIONS=0  # 0 = infinite

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
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

# Function to display usage
usage() {
    cat << EOF
Usage: $0 -i <input_dir> -o <output_dir> -r <reference> -s <sample_id> [OPTIONS]

Required arguments:
    -i, --input-dir DIR          Directory containing BAM files from ONT basecalling
    -o, --output-dir DIR         Output directory for CNV results
    -r, --reference FILE         Reference genome FASTA file
    -s, --sample-id ID           Sample identifier

Optional arguments:
    -c, --check-interval SEC     Interval to check for new BAM files (default: 300 seconds)
    -m, --min-size BYTES         Minimum BAM file size to process (default: 10000000)
    -e, --conda-env ENV          Conda environment for CNV analysis (default: cnvgen)
    --bin-sizes "SIZE1 SIZE2"    Bin sizes for CNVpytor (default: "1000 10000 100000")
    -t, --threads INT            Number of threads (default: 4)
    --max-iterations INT         Maximum number of iterations (0 = infinite, default: 0)
    --no-continuous              Process existing files once and exit
    -h, --help                   Display this help message

Example:
    # Continuous monitoring
    $0 -i /data/bam_dir -o /data/cnv_output -r /ref/hg38.fa -s SAMPLE001

    # Process existing files once
    $0 -i /data/bam_dir -o /data/cnv_output -r /ref/hg38.fa -s SAMPLE001 --no-continuous

    # Custom monitoring interval
    $0 -i /data/bam_dir -o /data/cnv_output -r /ref/hg38.fa -s SAMPLE001 -c 600
EOF
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
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
        -c|--check-interval)
            CHECK_INTERVAL="$2"
            shift 2
            ;;
        -m|--min-size)
            MIN_BAM_SIZE="$2"
            shift 2
            ;;
        -e|--conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --bin-sizes)
            BIN_SIZES="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --max-iterations)
            MAX_ITERATIONS="$2"
            shift 2
            ;;
        --no-continuous)
            CONTINUOUS=false
            shift
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
if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]] || [[ -z "$REFERENCE" ]] || [[ -z "$SAMPLE_ID" ]]; then
    log_error "Missing required arguments"
    usage
fi

# Validate input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
    log_error "Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Validate reference file exists
if [[ ! -f "$REFERENCE" ]]; then
    log_error "Reference file does not exist: $REFERENCE"
    exit 1
fi

# Create output directory structure
mkdir -p "$OUTPUT_DIR"/{individual_bams,timeseries_plots,combined_analysis,logs}
PROCESS_LOG="${OUTPUT_DIR}/logs/${SAMPLE_ID}_processing_log.txt"

# Initialize processed BAM tracking file
PROCESSED_BAMS="${OUTPUT_DIR}/logs/${SAMPLE_ID}_processed_bams.txt"
touch "$PROCESSED_BAMS"

# Log file for this session
exec > >(tee -a "$PROCESS_LOG")
exec 2>&1

log_info "========================================"
log_info "ONT Time-Series CNV Monitor Started"
log_info "========================================"
log_info "Input Directory: $INPUT_DIR"
log_info "Output Directory: $OUTPUT_DIR"
log_info "Reference: $REFERENCE"
log_info "Sample ID: $SAMPLE_ID"
log_info "Conda Environment: $CONDA_ENV"
log_info "Bin Sizes: $BIN_SIZES"
log_info "Check Interval: ${CHECK_INTERVAL}s"
log_info "Min BAM Size: $MIN_BAM_SIZE bytes"
log_info "Continuous Mode: $CONTINUOUS"
log_info "Processing log: $PROCESS_LOG"
log_info "========================================"

# Function to check if BAM is already processed
is_processed() {
    local bam_file="$1"
    grep -qxF "$bam_file" "$PROCESSED_BAMS" 2>/dev/null
}

# Function to mark BAM as processed
mark_processed() {
    local bam_file="$1"
    echo "$bam_file" >> "$PROCESSED_BAMS"
}

# Function to validate BAM file
validate_bam() {
    local bam_file="$1"

    # Check if file exists
    if [[ ! -f "$bam_file" ]]; then
        return 1
    fi

    # Check file size
    local file_size=$(stat -c%s "$bam_file" 2>/dev/null || echo 0)
    if [[ $file_size -lt $MIN_BAM_SIZE ]]; then
        log_warn "BAM file too small: $(basename "$bam_file") ($file_size bytes < $MIN_BAM_SIZE bytes)"
        return 1
    fi

    # Check if BAM index exists, create if not
    if [[ ! -f "${bam_file}.bai" ]]; then
        log_info "Creating BAM index: ${bam_file}.bai"
        if samtools index "$bam_file"; then
            log_success "BAM index created"
        else
            log_error "Failed to create BAM index"
            return 1
        fi
    fi

    return 0
}

# Function to process a single BAM file
process_bam() {
    local bam_file="$1"
    local bam_basename=$(basename "$bam_file" .bam)
    local timestamp=$(date '+%Y%m%d_%H%M%S')
    local output_prefix="${SAMPLE_ID}_${bam_basename}_${timestamp}"
    local bam_output_dir="${OUTPUT_DIR}/individual_bams/${output_prefix}"

    mkdir -p "$bam_output_dir"

    log_info "========================================" log_info "Processing BAM: $(basename "$bam_file")"
    log_info "Output prefix: $output_prefix"
    log_info "Output directory: $bam_output_dir"
    log_info "========================================"

    # Get script directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    # Call the individual CNV processing script
    if "${SCRIPT_DIR}/ont_process_individual_bam_cnv.sh" \
        -b "$bam_file" \
        -o "$bam_output_dir" \
        -r "$REFERENCE" \
        -s "$output_prefix" \
        --bin-sizes "$BIN_SIZES" \
        -e "$CONDA_ENV" \
        -t "$THREADS"; then

        log_success "Successfully processed: $(basename "$bam_file")"
        mark_processed "$bam_file"

        # Trigger time-series visualization update
        update_timeseries

        # Generate combined analysis
        generate_combined_analysis

        return 0
    else
        log_error "Failed to process: $(basename "$bam_file")"
        return 1
    fi
}

# Function to update time-series visualization
update_timeseries() {
    log_info "Updating time-series CNV visualization..."

    # Get script directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    if "${SCRIPT_DIR}/ont_timeseries_cnv_visualization.R" \
        --input-dir "${OUTPUT_DIR}/individual_bams" \
        --output-dir "${OUTPUT_DIR}/timeseries_plots" \
        --sample-id "$SAMPLE_ID"; then

        log_success "Time-series visualization updated"
    else
        log_warn "Failed to update time-series visualization (this is normal if only one timepoint exists)"
    fi
}

# Function to generate combined CNV analysis
generate_combined_analysis() {
    log_info "Generating combined CNV analysis..."

    # Get script directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    if "${SCRIPT_DIR}/ont_combined_cnv_analysis.R" \
        --input-dir "${OUTPUT_DIR}/individual_bams" \
        --output-dir "${OUTPUT_DIR}/combined_analysis" \
        --sample-id "$SAMPLE_ID"; then

        log_success "Combined CNV analysis completed"
    else
        log_warn "Failed to generate combined analysis (this is normal if only one timepoint exists)"
    fi
}

# Main monitoring loop
log_info "Starting monitoring loop..."

iteration=0
while true; do
    iteration=$((iteration + 1))
    log_info "========================================"
    log_info "Check iteration: $iteration"
    log_info "========================================"

    # Find all BAM files in the input directory
    shopt -s nullglob
    bam_files=("$INPUT_DIR"/*.bam)
    shopt -u nullglob

    if [[ ${#bam_files[@]} -eq 0 ]]; then
        log_info "No BAM files found in $INPUT_DIR"
    else
        log_info "Found ${#bam_files[@]} BAM file(s)"

        # Process each BAM file
        processed_this_iteration=0
        for bam_file in "${bam_files[@]}"; do
            # Skip if already processed
            if is_processed "$bam_file"; then
                log_info "Already processed: $(basename "$bam_file")"
                continue
            fi

            # Validate BAM file
            if validate_bam "$bam_file"; then
                process_bam "$bam_file"
                processed_this_iteration=$((processed_this_iteration + 1))
            else
                log_warn "Skipping invalid BAM: $(basename "$bam_file")"
            fi
        done

        log_info "Processed $processed_this_iteration new BAM file(s) in this iteration"
    fi

    # Check if we should continue
    if [[ "$CONTINUOUS" == "false" ]]; then
        log_info "Non-continuous mode: Exiting after processing existing files"
        break
    fi

    if [[ $MAX_ITERATIONS -gt 0 ]] && [[ $iteration -ge $MAX_ITERATIONS ]]; then
        log_info "Reached maximum iterations ($MAX_ITERATIONS): Exiting"
        break
    fi

    # Report processed count
    processed_count=$(wc -l < "$PROCESSED_BAMS" 2>/dev/null || echo 0)
    log_info "Total BAM files processed so far: $processed_count"

    log_info "Waiting ${CHECK_INTERVAL}s before next check..."
    sleep "$CHECK_INTERVAL"
done

log_success "========================================"
log_success "ONT Time-Series CNV Monitor Finished"
log_success "========================================"
log_success "Total iterations: $iteration"
log_success "Total BAMs processed: $(wc -l < "$PROCESSED_BAMS" 2>/dev/null || echo 0)"
log_success "Results directory: $OUTPUT_DIR"
log_success "========================================"

exit 0
