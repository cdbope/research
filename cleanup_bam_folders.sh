#!/bin/bash

# Script to check folder sizes and cleanup bam/fastq files
# Usage: ./cleanup_bam_folders.sh <sample_ids_file> <base_path>
# Example: ./cleanup_bam_folders.sh sample_ids.txt /home/chbope/data

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <sample_ids_file> <base_path>"
    echo "Example: $0 sample_ids.txt /home/chbope/data"
    exit 1
fi

SAMPLE_IDS_FILE="$1"
BASE_PATH="$2"
OUTPUT_FILE="folder_sizes_report.txt"
CLEANUP_LOG="cleanup_log.txt"

# Check if sample IDs file exists
if [ ! -f "$SAMPLE_IDS_FILE" ]; then
    echo "Error: Sample IDs file '$SAMPLE_IDS_FILE' not found!"
    exit 1
fi

# Check if base path exists
if [ ! -d "$BASE_PATH" ]; then
    echo "Error: Base path '$BASE_PATH' not found!"
    exit 1
fi

# Initialize output files
echo "=== Folder Size Report $(date) ===" > "$OUTPUT_FILE"
echo "=== Cleanup Log $(date) ===" > "$CLEANUP_LOG"

echo "Starting folder size check and cleanup process..."
echo "Sample IDs file: $SAMPLE_IDS_FILE"
echo "Base path: $BASE_PATH"
echo "Output report: $OUTPUT_FILE"
echo "Cleanup log: $CLEANUP_LOG"
echo ""

# Counter for processed samples
processed=0
total=$(wc -l < "$SAMPLE_IDS_FILE")

# Read sample IDs file line by line
while IFS= read -r sample_id || [ -n "$sample_id" ]; do
    # Skip empty lines and comments
    if [[ -z "$sample_id" || "$sample_id" =~ ^# ]]; then
        continue
    fi
    
    # Remove any whitespace
    sample_id=$(echo "$sample_id" | xargs)
    
    processed=$((processed + 1))
    echo "Processing sample $processed/$total: $sample_id"
    
    # Find all matching directories for this sample_id
    # Pattern: /sample_id/*/bam_pass/
    search_pattern="$BASE_PATH/$sample_id/*/bam_pass"
    
    found_dirs=$(find "$BASE_PATH" -path "$search_pattern" -type d 2>/dev/null)
    
    if [ -z "$found_dirs" ]; then
        echo "  Warning: No bam_pass directories found for sample: $sample_id" | tee -a "$CLEANUP_LOG"
        echo "$sample_id : No bam_pass directory found" >> "$OUTPUT_FILE"
        continue
    fi
    
    # Process each found directory
    total_size=0
    dir_count=0
    
    while IFS= read -r bam_pass_dir; do
        if [ -d "$bam_pass_dir" ]; then
            dir_count=$((dir_count + 1))
            echo "  Found directory: $bam_pass_dir"
            
            # 1. Calculate folder size (in bytes, then convert to human readable)
            size_bytes=$(du -sb "$bam_pass_dir" 2>/dev/null | cut -f1)
            if [ -n "$size_bytes" ]; then
                total_size=$((total_size + size_bytes))
                size_human=$(du -sh "$bam_pass_dir" 2>/dev/null | cut -f1)
                echo "    Directory size: $size_human" | tee -a "$CLEANUP_LOG"
            else
                echo "    Error calculating size for: $bam_pass_dir" | tee -a "$CLEANUP_LOG"
            fi
            
            # 3. Delete files with specified extensions
            echo "    Cleaning up files in: $bam_pass_dir" | tee -a "$CLEANUP_LOG"
            
            # Count files before deletion
            bam_files=$(find "$bam_pass_dir" -name "*.bam" 2>/dev/null | wc -l)
            bai_files=$(find "$bam_pass_dir" -name "*.bam.fai" 2>/dev/null | wc -l)
            fastq_files=$(find "$bam_pass_dir" -name "*.fastq" -o -name "*.fastq.gz" 2>/dev/null | wc -l)
            
            echo "    Files to delete: $bam_files *.bam, $bai_files *.bam.fai, $fastq_files *.fastq/*fastq.gz" | tee -a "$CLEANUP_LOG"
            
            # Perform deletion
            deleted_count=0
            
            # Delete .bam files
            if [ $bam_files -gt 0 ]; then
                find "$bam_pass_dir" -name "*.bam" -type f -delete 2>/dev/null
                deleted_count=$((deleted_count + bam_files))
            fi
            
            # Delete .bam.fai files
            if [ $bai_files -gt 0 ]; then
                find "$bam_pass_dir" -name "*.bam.fai" -type f -delete 2>/dev/null
                deleted_count=$((deleted_count + bai_files))
            fi
            
            # Delete .fastq and .fastq.gz files
            if [ $fastq_files -gt 0 ]; then
                find "$bam_pass_dir" -name "*.fastq" -type f -delete 2>/dev/null
                find "$bam_pass_dir" -name "*.fastq.gz" -type f -delete 2>/dev/null
                deleted_count=$((deleted_count + fastq_files))
            fi
            
            echo "    Deleted $deleted_count files total" | tee -a "$CLEANUP_LOG"
            
        fi
    done <<< "$found_dirs"
    
    # 2. Write total size to report file
    if [ $total_size -gt 0 ]; then
        # Convert bytes to human readable
        if [ $total_size -gt 1073741824 ]; then
            size_display=$(echo "scale=2; $total_size / 1073741824" | bc)"GB"
        elif [ $total_size -gt 1048576 ]; then
            size_display=$(echo "scale=2; $total_size / 1048576" | bc)"MB"
        elif [ $total_size -gt 1024 ]; then
            size_display=$(echo "scale=2; $total_size / 1024" | bc)"KB"
        else
            size_display="${total_size}B"
        fi
        
        echo "$sample_id : $size_display ($dir_count directories)" >> "$OUTPUT_FILE"
        echo "  Total size for $sample_id: $size_display" | tee -a "$CLEANUP_LOG"
    else
        echo "$sample_id : 0B (no valid directories)" >> "$OUTPUT_FILE"
    fi
    
    echo "" | tee -a "$CLEANUP_LOG"
    
done < "$SAMPLE_IDS_FILE"

echo "=== Process completed ===" | tee -a "$CLEANUP_LOG"
echo "Processed $processed samples"
echo "Size report saved to: $OUTPUT_FILE"
echo "Cleanup log saved to: $CLEANUP_LOG"
echo ""
echo "Summary report:"
cat "$OUTPUT_FILE"
