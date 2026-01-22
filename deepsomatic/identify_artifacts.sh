#!/bin/bash

##############################################################################
# Artifact Identification Tool
# Compares DeepSomatic and ClairS-TO to identify potential artifacts
##############################################################################

if [ $# -lt 2 ]; then
    echo "Usage: $0 <deepsomatic.vcf.gz> <clairsto.vcf.gz> [output_prefix]"
    echo ""
    echo "This script identifies:"
    echo "  1. Variants unique to DeepSomatic (potential artifacts)"
    echo "  2. Variants common to both callers (high confidence)"
    echo "  3. Quality metrics for unique variants"
    exit 1
fi

DS_VCF=$1
CLAIR_VCF=$2
PREFIX=${3:-"comparison"}

echo "=========================================================================="
echo "Artifact Identification Analysis"
echo "=========================================================================="
echo "DeepSomatic VCF: ${DS_VCF}"
echo "ClairS-TO VCF: ${CLAIR_VCF}"
echo "Output prefix: ${PREFIX}"
echo ""

# Create temporary directory
TMP_DIR=$(mktemp -d)
trap "rm -rf ${TMP_DIR}" EXIT

# Extract variant positions
echo "Step 1: Extracting variant positions..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%CHROM:%POS:%REF:%ALT\n' ${DS_VCF} \
  | sort -k5,5 > ${TMP_DIR}/ds_variants.txt

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%CHROM:%POS:%REF:%ALT\n' ${CLAIR_VCF} \
  | sort -k5,5 > ${TMP_DIR}/clair_variants.txt

DS_COUNT=$(wc -l < ${TMP_DIR}/ds_variants.txt)
CLAIR_COUNT=$(wc -l < ${TMP_DIR}/clair_variants.txt)

echo "  DeepSomatic variants: ${DS_COUNT}"
echo "  ClairS-TO variants: ${CLAIR_COUNT}"

# Find common and unique variants
echo ""
echo "Step 2: Identifying common and unique variants..."

# Common variants (in both)
join -t $'\t' -1 5 -2 5 ${TMP_DIR}/ds_variants.txt ${TMP_DIR}/clair_variants.txt \
  | awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5}' \
  > ${PREFIX}_common_variants.txt

# DeepSomatic unique
comm -23 <(cut -f5 ${TMP_DIR}/ds_variants.txt) <(cut -f5 ${TMP_DIR}/clair_variants.txt) \
  > ${TMP_DIR}/ds_unique_keys.txt

grep -F -f ${TMP_DIR}/ds_unique_keys.txt ${TMP_DIR}/ds_variants.txt \
  | cut -f1-4 > ${PREFIX}_deepsomatic_unique.txt

# ClairS-TO unique
comm -13 <(cut -f5 ${TMP_DIR}/ds_variants.txt) <(cut -f5 ${TMP_DIR}/clair_variants.txt) \
  > ${TMP_DIR}/clair_unique_keys.txt

grep -F -f ${TMP_DIR}/clair_unique_keys.txt ${TMP_DIR}/clair_variants.txt \
  | cut -f1-4 > ${PREFIX}_clairsto_unique.txt

COMMON_COUNT=$(wc -l < ${PREFIX}_common_variants.txt)
DS_UNIQUE_COUNT=$(wc -l < ${PREFIX}_deepsomatic_unique.txt)
CLAIR_UNIQUE_COUNT=$(wc -l < ${PREFIX}_clairsto_unique.txt)

echo "  Common variants (both callers): ${COMMON_COUNT}"
echo "  DeepSomatic unique: ${DS_UNIQUE_COUNT}"
echo "  ClairS-TO unique: ${CLAIR_UNIQUE_COUNT}"

# Extract quality metrics for DeepSomatic unique variants
echo ""
echo "Step 3: Extracting quality metrics for DeepSomatic-unique variants..."

if [ ${DS_UNIQUE_COUNT} -gt 0 ]; then
    # Create BED file for unique variants
    awk '{print $1"\t"$2-1"\t"$2"\t"$1":"$2":"$3":"$4}' ${PREFIX}_deepsomatic_unique.txt \
      > ${TMP_DIR}/ds_unique.bed

    # Extract quality metrics
    bcftools view -R ${TMP_DIR}/ds_unique.bed ${DS_VCF} | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GQ\t%DP\t%AD\t%VAF]\n' \
      > ${PREFIX}_deepsomatic_unique_metrics.txt

    # Add header
    echo -e "Chr\tPos\tRef\tAlt\tGQ\tDP\tAD\tVAF" | cat - ${PREFIX}_deepsomatic_unique_metrics.txt \
      > ${PREFIX}_deepsomatic_unique_with_metrics.txt

    echo "  ✓ Metrics extracted for ${DS_UNIQUE_COUNT} DeepSomatic-unique variants"
fi

# Analyze quality metrics
echo ""
echo "Step 4: Analyzing artifact likelihood..."

if [ ${DS_UNIQUE_COUNT} -gt 0 ] && [ -f ${PREFIX}_deepsomatic_unique_with_metrics.txt ]; then
    # Count potential artifacts based on thresholds
    LOW_GQ=$(awk 'NR>1 && $5 < 30' ${PREFIX}_deepsomatic_unique_with_metrics.txt | wc -l)
    LOW_DP=$(awk 'NR>1 && $6 < 20' ${PREFIX}_deepsomatic_unique_with_metrics.txt | wc -l)
    LOW_VAF=$(awk 'NR>1 && $8 < 0.05' ${PREFIX}_deepsomatic_unique_with_metrics.txt | wc -l)

    echo "  Potential artifact flags in DeepSomatic-unique variants:"
    echo "    - Low GQ (<30): ${LOW_GQ} variants"
    echo "    - Low Depth (<20): ${LOW_DP} variants"
    echo "    - Low VAF (<5%): ${LOW_VAF} variants"
fi

# Create summary report
echo ""
echo "Step 5: Creating summary report..."

cat > ${PREFIX}_artifact_analysis_report.txt << EOF
Artifact Identification Analysis Report
Generated: $(date)

Input Files:
  DeepSomatic VCF: ${DS_VCF}
  ClairS-TO VCF:   ${CLAIR_VCF}

Variant Count Summary:
  DeepSomatic total:        ${DS_COUNT}
  ClairS-TO total:          ${CLAIR_COUNT}
  Common (both callers):    ${COMMON_COUNT}
  DeepSomatic unique:       ${DS_UNIQUE_COUNT}
  ClairS-TO unique:         ${CLAIR_UNIQUE_COUNT}

Concordance:
  Agreement rate: $((COMMON_COUNT * 100 / DS_COUNT))% of DeepSomatic variants
  Excess variants: ${DS_UNIQUE_COUNT} ($(((DS_UNIQUE_COUNT * 100) / DS_COUNT))% of DeepSomatic)

Artifact Likelihood (DeepSomatic-unique variants):
EOF

if [ ${DS_UNIQUE_COUNT} -gt 0 ] && [ -f ${PREFIX}_deepsomatic_unique_with_metrics.txt ]; then
    cat >> ${PREFIX}_artifact_analysis_report.txt << EOF
  Low GQ (<30):             ${LOW_GQ} variants ($((LOW_GQ * 100 / DS_UNIQUE_COUNT))%)
  Low Depth (<20):          ${LOW_DP} variants ($((LOW_DP * 100 / DS_UNIQUE_COUNT))%)
  Low VAF (<5%):            ${LOW_VAF} variants ($((LOW_VAF * 100 / DS_UNIQUE_COUNT))%)

Recommendation:
  - HIGH CONFIDENCE: ${COMMON_COUNT} variants (both callers agree)
  - REVIEW REQUIRED: ${DS_UNIQUE_COUNT} DeepSomatic-unique variants
  - LIKELY ARTIFACTS: ~$((LOW_GQ + LOW_DP + LOW_VAF)) variants with quality flags

Output Files:
  ${PREFIX}_common_variants.txt                    - High confidence variants
  ${PREFIX}_deepsomatic_unique.txt                 - Review for artifacts
  ${PREFIX}_deepsomatic_unique_with_metrics.txt    - With quality metrics
  ${PREFIX}_clairsto_unique.txt                    - Variants missed by DeepSomatic
  ${PREFIX}_artifact_analysis_report.txt           - This report
EOF
else
    echo "  No unique variants to analyze" >> ${PREFIX}_artifact_analysis_report.txt
fi

echo ""
echo "=========================================================================="
echo "ANALYSIS COMPLETE"
echo "=========================================================================="
echo ""
echo "Summary:"
echo "  • HIGH CONFIDENCE: ${COMMON_COUNT} variants (both callers)"
echo "  • REVIEW NEEDED: ${DS_UNIQUE_COUNT} DeepSomatic-unique variants"
if [ ${DS_UNIQUE_COUNT} -gt 0 ]; then
    echo "  • LIKELY ARTIFACTS: ~$((LOW_GQ + LOW_DP + LOW_VAF)) with quality flags"
fi
echo ""
echo "Output files:"
echo "  ${PREFIX}_common_variants.txt"
echo "  ${PREFIX}_deepsomatic_unique_with_metrics.txt"
echo "  ${PREFIX}_artifact_analysis_report.txt"
echo ""
echo "Next steps:"
echo "  1. Review ${PREFIX}_artifact_analysis_report.txt for summary"
echo "  2. Use ${PREFIX}_common_variants.txt for high-confidence set"
echo "  3. Manually review DeepSomatic-unique variants in IGV"
echo "  4. Apply enhanced filtering (deepsomatic_v2_enhanced.sh)"
echo "=========================================================================="

cat ${PREFIX}_artifact_analysis_report.txt
