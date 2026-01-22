#!/usr/bin/env bash
set -euo pipefail

# File containing sample IDs, one per line (e.g., sample_ids.txt)
SAMPLE_LIST="sample_ids.txt"

# Base directory where sample folders are located (edit as needed)
BASE_DIR="/path/to/samples"

while IFS= read -r SAMPLE_ID; do
    # skip empty lines
    [[ -z "$SAMPLE_ID" ]] && continue
    
    SAMPLE_DIR="${BASE_DIR}/${SAMPLE_ID}/clairsto_output/tmp"
    
    if [[ -d "$SAMPLE_DIR" ]]; then
        echo "Deleting tmp folder for sample ${SAMPLE_ID}..."
        rm -rf "$SAMPLE_DIR"
    else
        echo "No tmp folder found for sample ${SAMPLE_ID}, skipping."
    fi
done < "$SAMPLE_LIST"

echo "Done."

