#!/bin/bash
###############################################################################
# Apply Stricter AF Filter to Exclude Germline Heterozygous Region
#
# Current filter: AF 10-90% (includes germline heterozygous at ~50%)
# Stricter filter: Exclude 40-60% region (removes most germline heterozygous)
#
# Trade-off: May lose some true somatic variants at AF ~50%
#
# Usage: ./apply_stricter_af_filter.sh
###############################################################################

set -euo pipefail

INPUT_VCF_DIR="results/prepared_vcfs/filter_vcf"
OUTPUT_VCF_DIR="results/prepared_vcfs/stricter_af_filter"

mkdir -p "$OUTPUT_VCF_DIR"

echo "============================================================================"
echo "STRICTER AF FILTER (EXCLUDE 40-60% GERMLINE REGION)"
echo "============================================================================"
echo "Filter: (AF >= 0.10 AND AF < 0.40) OR (AF > 0.60 AND AF <= 0.90)"
echo "This removes germline heterozygous variants (AF ~50%)"
echo "============================================================================"
echo

count=0
for vcf in "$INPUT_VCF_DIR"/*.vcf; do
    filename=$(basename "$vcf")
    output_vcf="$OUTPUT_VCF_DIR/$filename"

    echo "Processing: $filename"

    # Count before filter
    n_before=$(grep -v '^#' "$vcf" | wc -l)

    # Apply stricter AF filter: exclude 40-60% region
    bcftools view "$vcf" | \
      bcftools filter -i '(AF >= 0.10 && AF < 0.40) || (AF > 0.60 && AF <= 0.90)' \
      > "$output_vcf"

    # Count after filter
    n_after=$(grep -v '^#' "$output_vcf" | wc -l)
    n_removed=$((n_before - n_after))

    if [ "$n_before" -gt 0 ]; then
        pct_removed=$((100 * n_removed / n_before))
    else
        pct_removed=0
    fi

    echo "  Before: $n_before variants"
    echo "  After:  $n_after variants"
    echo "  Removed: $n_removed variants (${pct_removed}% - likely germline heterozygous)"
    echo "  ✓ Filtered: $output_vcf"
    echo

    count=$((count + 1))
done

echo "============================================================================"
echo "Stricter AF filtering complete!"
echo "Processed $count files in: $OUTPUT_VCF_DIR"
echo "============================================================================"
echo
echo "⚠️  IMPORTANT: This may remove some true somatic variants at AF ~50%"
echo "    Compare results to your current filter_vcf/ to assess trade-offs"
