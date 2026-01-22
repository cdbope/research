#!/usr/bin/env bash

# === CONFIG ===
# Path to your sample list
SAMPLE_LIST="sampleid.txt"

# Path to your R script
R_SCRIPT="/home/chbope/Documents/nanopore/clone/nextflow/bin/merge_annotations_prospective24june25.R"

# Path to your OCCgenes.rds reference file
REFERENCE="/media/chbope/Manga/images_reference/reference/OCCgenes.rds"

# Input folder where input CSVs are located
INPUT_FOLDER="/home/chbope/extension/trash/t19-008/"

# Output folder where final results will be saved
OUTPUT_FOLDER="/home/chbope/extension/trash/test2"

# Create output folder if it doesn't exist
mkdir -p "$OUTPUT_FOLDER"

# === Loop over sample IDs ===
while read -r line; do
  # Skip empty lines or lines starting with '#' (optional)
  [[ -z "$line" || "$line" =~ ^# ]] && continue

  # Get the first column as sample_id
  sample_id=$(echo "$line" | awk '{print $1}')

  echo "Running Rscript for $sample_id ..."

  Rscript "$R_SCRIPT" \
    "${INPUT_FOLDER}/${sample_id}_merge_annotateandfilter.csv" \
    "${INPUT_FOLDER}/${sample_id}_occ_pileup_annotateandfilter.csv" \
    "${INPUT_FOLDER}/${sample_id}_annotateandfilter_clairsto.csv" \
    "${sample_id}_merge_annotation_filter_snvs_allcall.csv" \
    "$REFERENCE" \
    "${OUTPUT_FOLDER}/${sample_id}_merge_annotation_filter_snvs_allcall_filter.csv"

done < "$SAMPLE_LIST"

echo "All done!"

