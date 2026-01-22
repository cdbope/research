#!/bin/bash

##############################################################################
# DeepSomatic v2 with BED File Support
#
# This version adds support for BED file region filtering to match ClairS-TO
##############################################################################

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config_with_bed_example.sh"

if [ -f "$CONFIG_FILE" ]; then
    echo "Loading configuration from ${CONFIG_FILE}..."
    source "$CONFIG_FILE"
else
    echo "Error: Configuration file not found at ${CONFIG_FILE}"
    echo "Please create config.sh with required settings"
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
    echo "  To add BED file: Set BED_FILE=/path/to/file.bed in config.sh"
    BED_MOUNT=""
    REGIONS_PARAM=""
fi

# Clean up any previous intermediate files to avoid corruption
echo "Cleaning up previous intermediate files..."
sudo rm -rf ${OUTPUT_DIR}/intermediate_results_dir/*

# Check if Docker image exists locally
echo "Checking DeepSomatic Docker image..."
if ! sudo docker image inspect google/deepsomatic:"${BIN_VERSION}" > /dev/null 2>&1; then
    echo "Image not found locally. Pulling DeepSomatic Docker image..."
    sudo docker pull google/deepsomatic:"${BIN_VERSION}"
else
    echo "DeepSomatic image ${BIN_VERSION} found locally. Skipping pull."
fi

echo ""
echo "=========================================================================="
echo "DeepSomatic Variant Calling Configuration"
echo "=========================================================================="
echo "Sample ID:     ${SAMPLE_ID}"
echo "Model Type:    ${MODEL_TYPE}"
echo "PON Filtering: ${USE_PON_FILTERING}"
echo "BED Regions:   ${BED_FILE:-None (genome-wide)}"
echo "Input BAM:     ${INPUT_DIR}/${BAM_FILE}"
echo "Output Dir:    ${OUTPUT_DIR}"
echo "=========================================================================="
echo ""

echo "Running DeepSomatic variant calling..."
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
--output_vcf="${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz" \
--sample_name_tumor="${SAMPLE_ID}" \
--num_shards=${NUM_SHARDS} \
--logging_dir="${OUTPUT_DIR}/logs" \
--intermediate_results_dir="${OUTPUT_DIR}/intermediate_results_dir" \
--use_default_pon_filtering=${USE_PON_FILTERING}

# Check if DeepSomatic completed successfully
if [ $? -ne 0 ]; then
    echo "Error: DeepSomatic variant calling failed"
    exit 1
fi

echo "✓ DeepSomatic variant calling completed successfully"

# Navigate to output directory for annotation
cd ${OUTPUT_DIR}

# Count total variants before filtering
TOTAL_VARIANTS=$(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz | wc -l)
echo "Total variants called: ${TOTAL_VARIANTS}"

echo "Filtering VCF for PASS variants only..."
bcftools view -f PASS -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz

# Check if filtered VCF was created
if [ ! -f ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz ]; then
    echo "Error: Failed to create filtered VCF"
    exit 1
fi

PASS_VARIANTS=$(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz | wc -l)
echo "✓ PASS variants: ${PASS_VARIANTS} (filtered out $((TOTAL_VARIANTS - PASS_VARIANTS)))"

echo "Converting VCF to ANNOVAR input format..."
${ANNOVAR_DIR}/convert2annovar.pl ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz \
  --format vcf4 \
  --filter pass \
  --includeinfo \
  --outfile deepsomatic_To_snv_avinput

# Check if conversion was successful
if [ ! -f deepsomatic_To_snv_avinput ]; then
    echo "Error: Failed to convert VCF to ANNOVAR format"
    exit 1
fi

echo "Running ANNOVAR annotation..."
${ANNOVAR_DIR}/table_annovar.pl deepsomatic_To_snv_avinput \
  -outfile deepsomatic_To_snv_avinput_snv \
  -buildver hg38 \
  -protocol refGene,clinvar_20240611,cosmic100coding2024 \
  -operation g,f,f \
  ${HUMANDB_DIR} \
  -otherinfo

# Check if ANNOVAR annotation was successful
if [ ! -f deepsomatic_To_snv_avinput_snv.hg38_multianno.txt ]; then
    echo "Error: ANNOVAR annotation failed"
    exit 1
fi

echo "Filtering and extracting annotated variants..."
awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
  deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | awk '/exonic/ || /TERT/ || /Func.refGene/' \
  | awk '!/dist=166/' \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv

# Check if final output was created
if [ -f ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv ]; then
    FINAL_COUNT=$(($(wc -l < ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv) - 1))
    echo "✓ Annotated and filtered variants saved to: ${OUTPUT_DIR}/${SAMPLE_ID}_annotateandfilter_deep_somatic.csv"
    echo "  Final filtered variants: ${FINAL_COUNT}"
else
    echo "Warning: Final filtered output file was not created"
    exit 1
fi

echo ""
echo "=========================================================================="
echo "Pipeline Completed Successfully!"
echo "=========================================================================="
echo "Variant Count Summary:"
echo "  Total variants called:     ${TOTAL_VARIANTS}"
echo "  PASS filter:               ${PASS_VARIANTS}"
echo "  Final annotated variants:  ${FINAL_COUNT}"
echo ""
echo "Output files:"
echo "  - Raw VCF:           ${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz"
echo "  - PASS-filtered VCF: ${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz"
echo "  - ANNOVAR annotated: ${OUTPUT_DIR}/deepsomatic_To_snv_avinput_snv.hg38_multianno.txt"
echo "  - Final CSV:         ${OUTPUT_DIR}/${SAMPLE_ID}_annotateandfilter_deep_somatic.csv"
echo "=========================================================================="

if [ -n "${BED_FILE}" ]; then
    echo ""
    echo "Note: Variants were called only in BED file regions"
    echo "      This should produce results comparable to ClairS-TO"
else
    echo ""
    echo "⚠ Warning: No BED file used - genome-wide calling may find more variants"
    echo "           To match ClairS-TO, add BED_FILE to config.sh"
fi
