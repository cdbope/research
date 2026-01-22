#!/bin/bash
###############################################################################
# Step 1: Prepare VCF files for SURVIVOR merge
#
# For structural variants from Sniffles, we don't need bcftools norm.
# We just need to:
# 1. Decompress .vcf.gz files to plain .vcf (SURVIVOR prefers uncompressed)
# 2. Put them in a common directory for SURVIVOR
#
# Usage: ./01_prepare_vcfs.sh
###############################################################################

set -euo pipefail  # Exit on error, undefined vars, and pipeline failures

# Configuration
INPUT_VCF_DIR="/media/chbope/Expansion/200gbmsv"
OUTPUT_VCF_DIR="/home/chbope/extension/script/svmeta/results/prepared_vcfs/filter_vcf"

# Create output directory
mkdir -p "$OUTPUT_VCF_DIR"

echo "============================================================================"
echo "STEP 1: PREPARE VCF FILES FOR SURVIVOR"
echo "============================================================================"
echo "Input directory:  $INPUT_VCF_DIR"
echo "Output directory: $OUTPUT_VCF_DIR"
echo "============================================================================"
echo
echo "Note: Sniffles VCFs are already valid SV VCFs."
echo "We will apply the following filters:"
echo "  1. PASS variants only (quality filter)"
echo "  2. Somatic-enrichment filters:"
echo "     - Allele frequency: 0.10 <= AF <= 0.90 (removes likely germline)"
echo "     - Read support: SUPPORT >= 5 (removes low-confidence calls)"
echo

shopt -s nullglob  # So the loop just skips if no matches instead of using literal '*.vcf.gz'

count=0
total_variants=0
total_pass_variants=0
total_somatic_enriched=0

for vcf_gz in "$INPUT_VCF_DIR"/*.vcf.gz; do
    filename=$(basename "$vcf_gz")         # e.g. KM00.wf_sv.vcf.gz
    base_no_gz="${filename%.gz}"           # e.g. KM00.wf_sv.vcf
    output_vcf="$OUTPUT_VCF_DIR/$base_no_gz"

    echo "Processing: $filename"
    echo "  -> Filtering PASS + somatic-enriched variants to: $(basename "$output_vcf")"

    # Count total variants before filtering
    n_total=$(gunzip -c "$vcf_gz" | grep -v '^#' | wc -l)

    # Count PASS variants
    n_pass=$(gunzip -c "$vcf_gz" | bcftools view -f PASS | grep -v '^#' | wc -l)

    # Apply PASS + somatic-enrichment filters
    # Note: bcftools isec requires bgzip-compressed input, so we use a different approach
    # Filter first, then remove variants that overlap with gnomAD
    gunzip -c "$vcf_gz" | \
        bcftools view -f PASS | \
        bcftools filter -i 'AF >= 0.10 && AF <= 0.90 && INFO/SUPPORT >= 5' | \
        bgzip -c > "$OUTPUT_VCF_DIR/temp_${filename}"

    # Index the temporary file
    bcftools index -t "$OUTPUT_VCF_DIR/temp_${filename}"

    # Remove gnomAD common variants using isec
    # -C: output variants present in first file but NOT in second
    # -w 1: write output from first file
    bcftools isec -C \
        -w 1 \
        "$OUTPUT_VCF_DIR/temp_${filename}" \
        /home/chbope/extension/script/svmeta/external_datasets/gnomad.v4.1.sv.sites.vcf.gz \
        > "$output_vcf"

    # Clean up temp files
    rm "$OUTPUT_VCF_DIR/temp_${filename}" "$OUTPUT_VCF_DIR/temp_${filename}.tbi"

    # Count somatic-enriched variants after all filters
    n_somatic=$(grep -v '^#' "$output_vcf" | wc -l)

    # Calculate percentages
    if [ "$n_total" -gt 0 ]; then
        pass_pct=$((n_pass * 100 / n_total))
        somatic_pct=$((n_somatic * 100 / n_total))
    else
        pass_pct=0
        somatic_pct=0
    fi

    echo "     Total variants:           $n_total"
    echo "     PASS variants:            $n_pass (${pass_pct}%)"
    echo "     Somatic-enriched (final): $n_somatic (${somatic_pct}%)"
    echo "     Filtered out:             $((n_total - n_somatic))"

    # Sanity check that file is non-empty
    if [ ! -s "$output_vcf" ]; then
        echo "  ⚠ Output VCF is empty: $output_vcf"
        exit 1
    fi

    echo "  ✓ Prepared: $output_vcf"
    count=$((count + 1))
    total_variants=$((total_variants + n_total))
    total_pass_variants=$((total_pass_variants + n_pass))
    total_somatic_enriched=$((total_somatic_enriched + n_somatic))
done

if [ "$count" -eq 0 ]; then
    echo "No *.vcf.gz files found in $INPUT_VCF_DIR"
    exit 1
fi

echo
echo "============================================================================"
echo "VCF preparation complete!"
echo "============================================================================"
echo "Prepared $count VCF files in: $OUTPUT_VCF_DIR"
echo
echo "FILTERING SUMMARY:"
echo "  Total variants across all samples:      $total_variants"
echo "  PASS variants:                           $total_pass_variants"
echo "  Somatic-enriched (final):                $total_somatic_enriched"
if [ "$total_variants" -gt 0 ]; then
    overall_pass_pct=$((total_pass_variants * 100 / total_variants))
    overall_somatic_pct=$((total_somatic_enriched * 100 / total_variants))
    pass_filtered=$((total_variants - total_pass_variants))
    somatic_filtered=$((total_pass_variants - total_somatic_enriched))
    total_filtered=$((total_variants - total_somatic_enriched))
    echo "  Overall PASS rate:                       ${overall_pass_pct}%"
    echo "  Overall somatic-enriched rate:           ${overall_somatic_pct}%"
    echo "  ------------------------------------------------"
    echo "  Variants filtered (quality):             $pass_filtered"
    echo "  Variants filtered (germline/low-conf):   $somatic_filtered"
    echo "  Total filtered out:                      $total_filtered"
fi
echo "============================================================================"
echo
echo "Next step: Run 02_merge_with_survivor.sh"
