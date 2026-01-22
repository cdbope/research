#!/bin/bash

# Base directory containing subfolders
BASE_DIR="/data/GBM_WGS/occupdate17112025/OCC1"

# Check if base directory exists
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Directory $BASE_DIR does not exist"
    exit 1
fi

echo "Searching for 'output_clair3' folders in $BASE_DIR..."
echo "=========================================="

# Counter for deleted folders
count=0

# Find and delete output_clair3 folders
for subfolder in "$BASE_DIR"/*; do
    if [ -d "$subfolder" ]; then
        sample_name=$(basename "$subfolder")
        output_clair3_path="$subfolder/output_clair3"
        
        if [ -d "$output_clair3_path" ]; then
            echo "Found: $output_clair3_path"
            rm -rf "$output_clair3_path"
            if [ $? -eq 0 ]; then
                echo "  ✓ Deleted successfully"
                ((count++))
            else
                echo "  ✗ Failed to delete"
            fi
        fi
    fi
done

echo "=========================================="
echo "Total 'output_clair3' folders deleted: $count"
