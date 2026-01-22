#!/bin/bash

# Usage: ./extract_reads_by_regions.sh <input_folder> <regions_bed> <output_folder>

set -euo pipefail

INPUT_FOLDER="$1"
REGIONS_BED="$2"
OUTPUT_FOLDER="$3"

if [[ ! -d "$INPUT_FOLDER" ]]; then
    echo "Error: Input folder '$INPUT_FOLDER' does not exist."
    exit 1
fi

if [[ ! -f "$REGIONS_BED" ]]; then
    echo "Error: BED file '$REGIONS_BED' does not exist."
    exit 1
fi

mkdir -p "$OUTPUT_FOLDER"

for BAM_FILE in "$INPUT_FOLDER"/*.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    OUTPUT_BAM="$OUTPUT_FOLDER/${SAMPLE_NAME}.roi.bam"

    echo "Processing $SAMPLE_NAME..."

    samtools  view -@ 4 -b -L "$REGIONS_BED" "$BAM_FILE" |
        samtools sort -@ 4  -o "$OUTPUT_BAM"

    samtools index -@ 4 "$OUTPUT_BAM"
done

echo "All BAM Files were processed and saved in: $OUTPUT_FOLDER"
