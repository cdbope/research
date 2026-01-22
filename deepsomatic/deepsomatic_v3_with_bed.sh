#!/bin/bash

##############################################################################
# DeepSomatic Pipeline v3 with BED File Support - Batch Processing
# This version processes multiple samples from a sample list file
# and supports BED file region filtering
##############################################################################

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config_v3.sh"

if [ -f "$CONFIG_FILE" ]; then
    echo "Loading configuration from ${CONFIG_FILE}..."
    source "$CONFIG_FILE"
else
    echo "Error: Configuration file not found at ${CONFIG_FILE}"
    echo "Please create config_v3.sh with required settings"
    exit 1
fi

# Validate required paths
if [ "${ANNOVAR_DIR}" == "/path/to/annovar" ] || [ "${HUMANDB_DIR}" == "/path/to/annovar/humandb" ]; then
    echo "Error: Please update ANNOVAR paths in ${CONFIG_FILE}"
    exit 1
fi

if [ ! -d "${ANNOVAR_DIR}" ]; then
    echo "Error: ANNOVAR directory not found: ${ANNOVAR_DIR}"
    exit 1
fi

if [ ! -f "${ANNOVAR_DIR}/convert2annovar.pl" ]; then
    echo "Error: convert2annovar.pl not found in ${ANNOVAR_DIR}"
    exit 1
fi

# Check if sample list file exists
if [ ! -f "${SAMPLE_LIST_FILE}" ]; then
    echo "Error: Sample list file not found: ${SAMPLE_LIST_FILE}"
    echo "Please create a file with one sample ID per line"
    exit 1
fi

# Check BED file configuration
if [ -n "${BED_FILE}" ] && [ -f "${BED_FILE}" ]; then
    BED_DIR=$(dirname "${BED_FILE}")
    echo "✓ Using BED file for region filtering: ${BED_FILE}"
    echo "  BED regions: $(wc -l < ${BED_FILE}) regions"
    BED_MOUNT="-v ${BED_DIR}:${BED_DIR}:ro"
    REGIONS_PARAM="--regions=${BED_FILE}"
else
    echo "⚠ No BED file specified or file not found."
    echo "  Calling variants genome-wide (may find more variants than ClairS-TO)"
    echo "  To add BED file: Set BED_FILE=/path/to/file.bed in config_v3.sh"
    BED_MOUNT=""
    REGIONS_PARAM=""
fi

# Count samples
TOTAL_SAMPLES=$(grep -v '^#' "${SAMPLE_LIST_FILE}" | grep -v '^$' | wc -l)
echo ""
echo "=========================================================================="
echo "DeepSomatic Pipeline v3 - Batch Mode with BED Support"
echo "=========================================================================="
echo "Configuration:"
echo "  Samples to process: ${TOTAL_SAMPLES}"
echo "  Model Type: ${MODEL_TYPE}"
echo "  PON Filtering: ${USE_PON_FILTERING}"
echo "  BED File: ${BED_FILE:-None (genome-wide)}"
echo "  Input Directory: ${INPUT_DIR}"
echo "  Output Directory: ${OUTPUT_DIR}"
echo "=========================================================================="
echo ""

# Check if Docker image exists locally (only once for all samples)
echo "Checking DeepSomatic Docker image..."
if ! sudo docker image inspect google/deepsomatic:"${BIN_VERSION}" > /dev/null 2>&1; then
    echo "Image not found locally. Pulling DeepSomatic Docker image..."
    sudo docker pull google/deepsomatic:"${BIN_VERSION}"
else
    echo "DeepSomatic image ${BIN_VERSION} found locally. Skipping pull."
fi
echo ""

# Initialize counters
CURRENT_SAMPLE=0
SUCCESS_COUNT=0
FAILED_COUNT=0
FAILED_SAMPLES=()

# Process each sample
while IFS= read -r SAMPLE_ID || [ -n "$SAMPLE_ID" ]; do
    # Skip empty lines and comments
    [[ -z "$SAMPLE_ID" || "$SAMPLE_ID" =~ ^#.*$ ]] && continue

    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))

    echo ""
    echo "=========================================================================="
    echo "Processing sample ${CURRENT_SAMPLE}/${TOTAL_SAMPLES}: ${SAMPLE_ID}"
    echo "=========================================================================="

    # Construct BAM filename
    BAM_FILE="${SAMPLE_ID}${BAM_SUFFIX}"

    # Check if BAM file exists
    if [ ! -f "${INPUT_DIR}/${BAM_FILE}" ]; then
        echo "Error: BAM file not found: ${INPUT_DIR}/${BAM_FILE}"
        echo "Skipping sample ${SAMPLE_ID}..."
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_SAMPLES+=("${SAMPLE_ID} (BAM not found)")
        continue
    fi

    echo "BAM file: ${BAM_FILE}"
    echo "Start time: $(date)"

    # Create sample-specific output directory
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
    mkdir -p "${SAMPLE_OUTPUT_DIR}"

    # Clean up any previous intermediate files for this sample
    echo "Cleaning up previous intermediate files for ${SAMPLE_ID}..."
    rm -rf ${SAMPLE_OUTPUT_DIR}/intermediate_results_dir/*

    # Step 1: Run DeepSomatic variant calling
    echo ""
    echo "Step 1/5: Running DeepSomatic variant calling..."
    sudo docker run \
    -v "${INPUT_DIR}":"${INPUT_DIR}":ro \
    -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
    -v "${REF_DIR}":"${REF_DIR}":ro \
    ${BED_MOUNT} \
    google/deepsomatic:"${BIN_VERSION}" \
    run_deepsomatic \
    --model_type=${MODEL_TYPE} \
    --ref="${REF_DIR}/${REF_GENOME}" \
    --reads_tumor="${INPUT_DIR}/${BAM_FILE}" \
    ${REGIONS_PARAM} \
    --output_vcf="${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz" \
    --sample_name_tumor="${SAMPLE_ID}" \
    --num_shards=${NUM_SHARDS} \
    --logging_dir="${SAMPLE_OUTPUT_DIR}/logs" \
    --intermediate_results_dir="${SAMPLE_OUTPUT_DIR}/intermediate_results_dir" \
    --use_default_pon_filtering=${USE_PON_FILTERING}

    # Check if DeepSomatic completed successfully
    if [ $? -ne 0 ]; then
        echo "Error: DeepSomatic variant calling failed for ${SAMPLE_ID}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_SAMPLES+=("${SAMPLE_ID} (DeepSomatic failed)")
        continue
    fi

    echo "✓ DeepSomatic variant calling completed for ${SAMPLE_ID}"

    # Navigate to sample output directory
    cd ${SAMPLE_OUTPUT_DIR}

    # Count total variants before filtering
    TOTAL_VARIANTS=$(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz 2>/dev/null | wc -l)
    echo "  Total variants called: ${TOTAL_VARIANTS}"

    # Step 2: Filter VCF for PASS variants only
    echo ""
    echo "Step 2/5: Filtering VCF for PASS variants only..."
    bcftools view -f PASS -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz

    if [ ! -f ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz ]; then
        echo "Error: Failed to create filtered VCF for ${SAMPLE_ID}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_SAMPLES+=("${SAMPLE_ID} (VCF filtering failed)")
        continue
    fi

    PASS_VARIANTS=$(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz 2>/dev/null | wc -l)
    echo "  ✓ PASS variants: ${PASS_VARIANTS}"

    # Step 3: Convert to ANNOVAR format
    echo ""
    echo "Step 3/5: Converting VCF to ANNOVAR input format..."
    ${ANNOVAR_DIR}/convert2annovar.pl ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz \
      --format vcf4 \
      --filter pass \
      --includeinfo \
      --outfile ${SAMPLE_ID}_deepsomatic_to_snv_avinput

    if [ ! -f ${SAMPLE_ID}_deepsomatic_to_snv_avinput ]; then
        echo "Error: Failed to convert VCF to ANNOVAR format for ${SAMPLE_ID}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_SAMPLES+=("${SAMPLE_ID} (ANNOVAR conversion failed)")
        continue
    fi

    echo "  ✓ Converted to ANNOVAR format"

    # Step 4: Run ANNOVAR annotation
    echo ""
    echo "Step 4/5: Running ANNOVAR annotation..."
    ${ANNOVAR_DIR}/table_annovar.pl ${SAMPLE_ID}_deepsomatic_to_snv_avinput \
      -outfile ${SAMPLE_ID}_deepsomatic_annotated \
      -buildver hg38 \
      -protocol refGene,clinvar_20240611,cosmic100coding2024 \
      -operation g,f,f \
      ${HUMANDB_DIR} \
      -otherinfo

    if [ ! -f ${SAMPLE_ID}_deepsomatic_annotated.hg38_multianno.txt ]; then
        echo "Error: ANNOVAR annotation failed for ${SAMPLE_ID}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_SAMPLES+=("${SAMPLE_ID} (ANNOVAR annotation failed)")
        continue
    fi

    echo "  ✓ Annotation completed"

    # Step 5: Filter and extract annotated variants
    echo ""
    echo "Step 5/5: Filtering and extracting annotated variants..."
    awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
      ${SAMPLE_ID}_deepsomatic_annotated.hg38_multianno.txt \
      | awk '/exonic/ || /TERT/ || /Func.refGene/' \
      | awk '!/dist=166/' \
      | cut -f1-16,25,26 > ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv

    if [ -f ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv ]; then
        VARIANT_COUNT=$(($(wc -l < ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv) - 1))
        echo "  ✓ Success! Found ${VARIANT_COUNT} filtered variants for ${SAMPLE_ID}"
        echo ""
        echo "  Variant count summary for ${SAMPLE_ID}:"
        echo "    Total called:     ${TOTAL_VARIANTS}"
        echo "    PASS filter:      ${PASS_VARIANTS}"
        echo "    Final annotated:  ${VARIANT_COUNT}"
        echo ""
        echo "  Output: ${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}_annotateandfilter_deep_somatic.csv"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "Warning: Final filtered output file was not created for ${SAMPLE_ID}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        FAILED_SAMPLES+=("${SAMPLE_ID} (Final filtering failed)")
        continue
    fi

    echo "Completed processing ${SAMPLE_ID} at $(date)"
    echo "=========================================================================="

done < "${SAMPLE_LIST_FILE}"

# Print final summary
echo ""
echo "=========================================================================="
echo "BATCH PROCESSING SUMMARY"
echo "=========================================================================="
echo "Total samples:          ${TOTAL_SAMPLES}"
echo "Successfully processed: ${SUCCESS_COUNT}"
echo "Failed:                 ${FAILED_COUNT}"
echo ""

if [ ${FAILED_COUNT} -gt 0 ]; then
    echo "Failed samples:"
    for failed_sample in "${FAILED_SAMPLES[@]}"; do
        echo "  - ${failed_sample}"
    done
    echo ""
fi

echo "Configuration used:"
echo "  BED File: ${BED_FILE:-None (genome-wide)}"
echo "  PON Filtering: ${USE_PON_FILTERING}"
echo "  Model: ${MODEL_TYPE}"
echo ""
echo "All results saved in: ${OUTPUT_DIR}"
echo "=========================================================================="

if [ -n "${BED_FILE}" ]; then
    echo ""
    echo "Note: Variants were called only in BED file regions"
    echo "      This should produce results comparable to ClairS-TO"
else
    echo ""
    echo "⚠ Warning: No BED file used - genome-wide calling may find more variants"
    echo "           To match ClairS-TO, add BED_FILE to config_v3.sh"
fi

echo ""
echo "Batch processing completed at: $(date)"

# Exit with error if any samples failed
if [ ${FAILED_COUNT} -gt 0 ]; then
    exit 1
else
    echo "All samples processed successfully!"
    exit 0
fi
