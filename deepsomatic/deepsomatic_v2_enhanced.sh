#!/bin/bash

##############################################################################
# DeepSomatic v2 Enhanced - Single Sample with Artifact Filtering
#
# Enhancements:
# - Quality metric filtering (GQ, DP, VAF)
# - PON filtering enabled
# - Tiered variant output (High/Medium/Low confidence)
# - Summary statistics
##############################################################################

# Exit on error
set -e

# Configuration
BIN_VERSION="1.8.0"
MODEL_TYPE="ONT_TUMOR_ONLY"
NUM_SHARDS=$(nproc)
USE_PON_FILTERING=true  # Enable PON filtering to remove artifacts

# Directories
REF_DIR="/home/chbope/extension/nWGS_manuscript_data/data/reference"
REF_GENOME="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
ANNOVAR_DIR="/home/chbope/extension/annovar"
HUMANDB_DIR="/home/chbope/extension/annovar/humandb"

# Quality thresholds for artifact filtering
MIN_GQ=30        # Minimum genotype quality
MIN_DP=20        # Minimum read depth
MIN_VAF=0.05     # Minimum variant allele frequency (5%)

# Input/Output (provide as arguments or set here)
if [ $# -lt 3 ]; then
    echo "Usage: $0 <SAMPLE_ID> <INPUT_DIR> <BAM_FILE> [OUTPUT_DIR]"
    echo ""
    echo "Example: $0 T25-152 /path/to/bams sample.bam /path/to/output"
    echo ""
    echo "Arguments:"
    echo "  SAMPLE_ID  : Sample identifier (e.g., T25-152)"
    echo "  INPUT_DIR  : Directory containing BAM file"
    echo "  BAM_FILE   : BAM filename (e.g., sample.bam)"
    echo "  OUTPUT_DIR : Output directory (optional, default: ./deepsomatic_output_SAMPLE_ID)"
    exit 1
fi

SAMPLE_ID=$1
INPUT_DIR=$2
BAM_FILE=$3
OUTPUT_DIR=${4:-"./deepsomatic_output_${SAMPLE_ID}"}

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=========================================================================="
echo "DeepSomatic Enhanced Pipeline with Artifact Filtering"
echo "=========================================================================="
echo "Sample ID: ${SAMPLE_ID}"
echo "Input BAM: ${INPUT_DIR}/${BAM_FILE}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Reference: ${REF_DIR}/${REF_GENOME}"
echo ""
echo "Quality Thresholds:"
echo "  - Minimum GQ: ${MIN_GQ}"
echo "  - Minimum DP: ${MIN_DP}"
echo "  - Minimum VAF: ${MIN_VAF}"
echo "  - PON Filtering: ${USE_PON_FILTERING}"
echo "=========================================================================="

# Verify input files exist
if [ ! -f "${INPUT_DIR}/${BAM_FILE}" ]; then
    echo "Error: BAM file not found: ${INPUT_DIR}/${BAM_FILE}"
    exit 1
fi

if [ ! -f "${REF_DIR}/${REF_GENOME}" ]; then
    echo "Error: Reference genome not found: ${REF_DIR}/${REF_GENOME}"
    exit 1
fi

# Check if Docker image exists, if not pull it
if ! sudo docker image inspect google/deepsomatic:"${BIN_VERSION}" > /dev/null 2>&1; then
    echo "Pulling DeepSomatic Docker image..."
    sudo docker pull google/deepsomatic:"${BIN_VERSION}"
fi

# Clean up previous intermediate results to avoid corrupted files
if [ -d "${OUTPUT_DIR}/intermediate_results_dir" ]; then
    echo "Cleaning up previous intermediate results..."
    sudo rm -rf ${OUTPUT_DIR}/intermediate_results_dir/*
fi

# Step 1: Run DeepSomatic variant calling
echo ""
echo "Step 1/7: Running DeepSomatic variant calling..."
echo "----------------------------------------------------------------------"

sudo docker run \
  -v "${INPUT_DIR}:${INPUT_DIR}:ro" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  -v "${REF_DIR}:${REF_DIR}:ro" \
  google/deepsomatic:"${BIN_VERSION}" \
  run_deepsomatic \
  --model_type=${MODEL_TYPE} \
  --ref="${REF_DIR}/${REF_GENOME}" \
  --reads_tumor="${INPUT_DIR}/${BAM_FILE}" \
  --output_vcf="${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz" \
  --sample_name_tumor="${SAMPLE_ID}" \
  --num_shards=${NUM_SHARDS} \
  --logging_dir="${OUTPUT_DIR}/logs" \
  --intermediate_results_dir="${OUTPUT_DIR}/intermediate_results_dir" \
  --use_default_pon_filtering=${USE_PON_FILTERING}

if [ $? -ne 0 ]; then
    echo "Error: DeepSomatic variant calling failed"
    exit 1
fi

echo "✓ DeepSomatic variant calling completed"

cd ${OUTPUT_DIR}

# Step 2: Filter for PASS variants
echo ""
echo "Step 2/7: Filtering for PASS variants..."
echo "----------------------------------------------------------------------"

bcftools view -f PASS -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz \
  ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz

if [ ! -f ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz ]; then
    echo "Error: Failed to create PASS-filtered VCF"
    exit 1
fi

PASS_COUNT=$(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz | wc -l)
echo "✓ PASS variants: ${PASS_COUNT}"

# Step 3: Apply quality metric filtering
echo ""
echo "Step 3/7: Applying quality metric filters (GQ≥${MIN_GQ}, DP≥${MIN_DP})..."
echo "----------------------------------------------------------------------"

bcftools view ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz | \
bcftools filter -i "GQ>=${MIN_GQ} && INFO/DP>=${MIN_DP}" \
  -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz

bcftools index -t ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz

QC_COUNT=$(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz | wc -l)
echo "✓ Quality-filtered variants: ${QC_COUNT} (removed $((PASS_COUNT - QC_COUNT)) low-quality)"

# Step 4: Filter for adequate VAF
echo ""
echo "Step 4/7: Filtering for VAF ≥ ${MIN_VAF} (removing ultra-low frequency artifacts)..."
echo "----------------------------------------------------------------------"

# Create a filter expression for VAF
# Note: VAF might be in different FORMAT fields depending on DeepSomatic version
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GQ\t%DP\t%AD\t%VAF]\n' \
  ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz > ${SAMPLE_ID}_variants_with_metrics.txt

# Extract variants meeting VAF threshold
awk -F'\t' -v min_vaf=${MIN_VAF} 'NR==1 || $9 >= min_vaf' \
  ${SAMPLE_ID}_variants_with_metrics.txt > ${SAMPLE_ID}_vaf_pass.txt

# Create position file for filtering
awk 'NR>1 {print $1":"$2}' ${SAMPLE_ID}_vaf_pass.txt > ${SAMPLE_ID}_vaf_positions.txt

# Filter VCF to keep only variants with adequate VAF
if [ -s ${SAMPLE_ID}_vaf_positions.txt ]; then
    bcftools view -T ${SAMPLE_ID}_vaf_positions.txt \
      ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz \
      -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vaf.vcf.gz

    VAF_COUNT=$(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vaf.vcf.gz | wc -l)
    echo "✓ VAF-filtered variants: ${VAF_COUNT} (removed $((QC_COUNT - VAF_COUNT)) ultra-low VAF)"
else
    echo "⚠ No variants passed VAF filter. Using quality-filtered VCF instead."
    cp ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz \
       ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vaf.vcf.gz
    VAF_COUNT=${QC_COUNT}
fi

# Step 5: Convert to ANNOVAR format
echo ""
echo "Step 5/7: Converting to ANNOVAR format..."
echo "----------------------------------------------------------------------"

${ANNOVAR_DIR}/convert2annovar.pl ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vaf.vcf.gz \
  --format vcf4 \
  --filter pass \
  --includeinfo \
  --outfile deepsomatic_To_snv_avinput

if [ ! -f deepsomatic_To_snv_avinput ]; then
    echo "Error: Failed to convert VCF to ANNOVAR format"
    exit 1
fi

echo "✓ Converted to ANNOVAR format"

# Step 6: ANNOVAR annotation
echo ""
echo "Step 6/7: Running ANNOVAR annotation..."
echo "----------------------------------------------------------------------"

${ANNOVAR_DIR}/table_annovar.pl deepsomatic_To_snv_avinput \
  -outfile deepsomatic_To_snv_avinput_snv \
  -buildver hg38 \
  -protocol refGene,clinvar_20240611,cosmic100coding2024 \
  -operation g,f,f \
  ${HUMANDB_DIR} \
  -otherinfo

if [ ! -f deepsomatic_To_snv_avinput_snv.hg38_multianno.txt ]; then
    echo "Error: ANNOVAR annotation failed"
    exit 1
fi

echo "✓ ANNOVAR annotation completed"

# Step 7: Tiered filtering based on evidence level
echo ""
echo "Step 7/7: Creating tiered variant sets..."
echo "----------------------------------------------------------------------"

# Tier 1: HIGH CONFIDENCE - COSMIC + Pathogenic/Likely pathogenic
awk -F'\t' '
  /^Chr\t/ ||
  (/COSM/ && (/Pathogenic/ || /Likely_pathogenic/))
' deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_tier1_high_confidence.csv

TIER1_COUNT=$(($(wc -l < ${SAMPLE_ID}_tier1_high_confidence.csv) - 1))

# Tier 2: MEDIUM-HIGH CONFIDENCE - Database evidence (COSMIC OR ClinVar)
awk -F'\t' '
  /^Chr\t/ ||
  /COSM/ ||
  /Pathogenic/ ||
  /Likely_pathogenic/ ||
  /Uncertain_significance/
' deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_tier2_database_evidence.csv

TIER2_COUNT=$(($(wc -l < ${SAMPLE_ID}_tier2_database_evidence.csv) - 1))

# Tier 3: MEDIUM CONFIDENCE - Functional variants (exonic, TERT promoter)
awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
  deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | awk '/exonic/ || /TERT/ || /Func.refGene/' \
  | awk '!/dist=166/' \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_tier3_functional.csv

TIER3_COUNT=$(($(wc -l < ${SAMPLE_ID}_tier3_functional.csv) - 1))

# All annotated variants (original filtering for comparison)
awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
  deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | awk '/exonic/ || /TERT/ || /Func.refGene/' \
  | awk '!/dist=166/' \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv

ALL_COUNT=$(($(wc -l < ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv) - 1))

# Create summary report
echo ""
echo "=========================================================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "=========================================================================="
echo ""
echo "Filtering Summary:"
echo "----------------------------------------------------------------------"
echo "  Total raw variants:           $(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz | wc -l)"
echo "  PASS filter:                  ${PASS_COUNT}"
echo "  Quality filter (GQ≥${MIN_GQ}, DP≥${MIN_DP}): ${QC_COUNT}"
echo "  VAF filter (VAF≥${MIN_VAF}):         ${VAF_COUNT}"
echo ""
echo "Annotated Variant Tiers:"
echo "----------------------------------------------------------------------"
echo "  Tier 1 (HIGH - COSMIC+Pathogenic):     ${TIER1_COUNT} variants"
echo "  Tier 2 (MEDIUM-HIGH - Database):       ${TIER2_COUNT} variants"
echo "  Tier 3 (MEDIUM - Functional):          ${TIER3_COUNT} variants"
echo "  All filtered variants:                 ${ALL_COUNT} variants"
echo ""
echo "Output Files:"
echo "----------------------------------------------------------------------"
echo "  Raw VCF:"
echo "    ${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz"
echo ""
echo "  Filtered VCFs:"
echo "    ${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz"
echo "    ${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz"
echo "    ${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vaf.vcf.gz"
echo ""
echo "  Annotated Variants:"
echo "    ${OUTPUT_DIR}/${SAMPLE_ID}_tier1_high_confidence.csv          (RECOMMENDED)"
echo "    ${OUTPUT_DIR}/${SAMPLE_ID}_tier2_database_evidence.csv"
echo "    ${OUTPUT_DIR}/${SAMPLE_ID}_tier3_functional.csv"
echo "    ${OUTPUT_DIR}/${SAMPLE_ID}_annotateandfilter_deep_somatic.csv (all)"
echo ""
echo "  Full ANNOVAR output:"
echo "    ${OUTPUT_DIR}/deepsomatic_To_snv_avinput_snv.hg38_multianno.txt"
echo ""
echo "Recommendations:"
echo "----------------------------------------------------------------------"
echo "  • For clinical reporting: Use Tier 1 variants (highest confidence)"
echo "  • For research: Include Tier 2 variants (database-supported)"
echo "  • For validation: Review Tier 3 variants"
echo "  • Compare with ClairS-TO results for additional confidence"
echo "=========================================================================="

# Create filtering statistics file
cat > ${SAMPLE_ID}_filtering_stats.txt << EOF
DeepSomatic Enhanced Pipeline - Filtering Statistics
Sample: ${SAMPLE_ID}
Date: $(date)

Quality Thresholds Applied:
  - Minimum GQ: ${MIN_GQ}
  - Minimum DP: ${MIN_DP}
  - Minimum VAF: ${MIN_VAF}
  - PON Filtering: ${USE_PON_FILTERING}

Variant Counts by Filter Stage:
  Total raw variants:                 $(bcftools view -H ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz | wc -l)
  After PASS filter:                  ${PASS_COUNT}
  After quality filter (GQ≥${MIN_GQ}, DP≥${MIN_DP}):  ${QC_COUNT}
  After VAF filter (VAF≥${MIN_VAF}):             ${VAF_COUNT}

Annotated Variants by Evidence Tier:
  Tier 1 (HIGH - COSMIC+Pathogenic):     ${TIER1_COUNT}
  Tier 2 (MEDIUM-HIGH - Database):       ${TIER2_COUNT}
  Tier 3 (MEDIUM - Functional):          ${TIER3_COUNT}
  All filtered variants:                 ${ALL_COUNT}

Artifact Removal Estimate:
  Removed by quality filters: $((PASS_COUNT - VAF_COUNT)) variants
  Remaining artifacts estimate: Review Tier 3 variants without database support
EOF

echo "✓ Filtering statistics saved to: ${SAMPLE_ID}_filtering_stats.txt"
echo ""
echo "Pipeline completed at: $(date)"
