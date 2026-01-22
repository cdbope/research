# Artifact Filtering Strategy for DeepSomatic Variants

## Problem Statement

Benchmark results show DeepSomatic detects **100% more variants** than ClairS-TO (66 vs 33).
This raises the critical question: **Are the extra 35 variants real somatic mutations or artifacts?**

## Current Filtering (Both Pipelines)

### 1. VCF-level Filtering
```bash
bcftools view -f PASS
```
- Only keeps variants that passed all quality filters
- ✅ Good baseline filtering

### 2. Functional Filtering
```bash
awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/'
awk '/exonic/ || /TERT/ || /Func.refGene/'
awk '!/dist=166/'
```
- Keeps: exonic nonsynonymous, upstream, TERT promoter
- Removes: benign variants, specific distance variants
- ⚠️ **Problem:** Too permissive, no quality thresholds

## Recommended Additional Filtering Strategies

### Strategy 1: Quality Metrics Filtering (RECOMMENDED)

Add filtering based on VCF quality fields **BEFORE** ANNOVAR annotation:

```bash
# Extract quality metrics from VCF and filter
bcftools view -f PASS -O z ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz | \
bcftools filter -i 'GQ>=30 && DP>=20 && (FORMAT/VAF)>=0.05' \
  -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.filtered.vcf.gz
```

**Recommended Thresholds:**
- `GQ >= 30` - Genotype Quality (high confidence)
- `DP >= 20` - Minimum read depth (adequate coverage)
- `VAF >= 0.05` (5%) - Minimum allele frequency (avoid ultra-low VAF artifacts)

**Rationale:**
- From benchmark: DeepSomatic detects down to 5.36% VAF
- Many artifacts occur at very low VAF (<5%)
- Setting VAF ≥5% removes noise while keeping real low-frequency variants

### Strategy 2: Strand Bias Filtering

ONT sequencing can have strand-specific artifacts. Add strand bias filter:

```bash
# Filter out variants with severe strand bias
bcftools filter -i 'SB < 0.001 || SB == "."' \
  -O z -o ${SAMPLE_ID}_filtered.vcf.gz
```

### Strategy 3: Multi-Caller Consensus (HIGH CONFIDENCE SET)

Use variants detected by **BOTH** callers as high-confidence set:

```bash
# Create intersection of DeepSomatic and ClairS-TO
bcftools isec -p intersection_dir \
  deepsomatic.vcf.gz \
  clairsto.vcf.gz

# Variants in both callers (high confidence)
cp intersection_dir/0002.vcf.gz ${SAMPLE_ID}_high_confidence.vcf.gz
```

**From your benchmark:**
- 31 variants agreed by both = **High confidence**
- 35 DeepSomatic-only = **Need additional validation**

### Strategy 4: Database-Based Filtering

Keep variants with evidence in cancer databases:

```bash
# After ANNOVAR, prioritize variants with database support
awk '
  /Func.refGene/ ||
  /COSM/ ||
  /Pathogenic/ ||
  /Likely_pathogenic/ ||
  (/exonic/ && /nonsynonymous/ && $GQ >= 30)
' annotated.txt > high_confidence.csv
```

**Priority levels:**
1. **Tier 1 (Highest confidence):** COSMIC + Pathogenic/Likely pathogenic
2. **Tier 2 (High confidence):** COSMIC OR ClinVar
3. **Tier 3 (Medium confidence):** Exonic nonsynonymous with GQ≥30, DP≥20
4. **Tier 4 (Low confidence):** Other variants (requires validation)

### Strategy 5: Panel of Normals (PON) Filtering

You're already using this, but verify it's enabled:

```bash
--use_default_pon_filtering=true
```

**What it does:**
- Removes common germline variants
- Filters recurrent sequencing artifacts
- Based on normal samples database

⚠️ **Check:** Your script has `${USE_PON_FILTERING}` - make sure this is set to `true`!

## Recommended Enhanced Pipeline

Here's an improved version of your DeepSomatic pipeline with artifact filtering:

```bash
#!/bin/bash

# [Previous setup code...]

# Step 1: Run DeepSomatic with PON filtering
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
  --use_default_pon_filtering=true  # ← ENSURE THIS IS TRUE

cd ${OUTPUT_DIR}

# Step 2: Filter for PASS + Quality metrics
echo "Filtering VCF for PASS variants with quality thresholds..."
bcftools view -f PASS ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz | \
bcftools filter -i 'GQ>=30 && INFO/DP>=20' \
  -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz

# Index the filtered VCF
bcftools index -t ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz

# Step 3: Extract variants with adequate VAF (≥5%)
echo "Filtering for VAF >= 5%..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GQ\t%DP\t%AD\t%VAF]\n' \
  ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz | \
awk -F'\t' '$8 >= 0.05' > ${SAMPLE_ID}_vaf_filtered.txt

# Create filtered VCF with VAF >= 5%
bcftools view -i 'FORMAT/VAF >= 0.05' \
  ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vcf.gz \
  -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vaf.vcf.gz

# Step 4: ANNOVAR annotation
echo "Converting to ANNOVAR format..."
${ANNOVAR_DIR}/convert2annovar.pl ${SAMPLE_ID}_ont_deepsomatic_output.pass.qc.vaf.vcf.gz \
  --format vcf4 \
  --filter pass \
  --includeinfo \
  --outfile deepsomatic_To_snv_avinput

echo "Running ANNOVAR annotation..."
${ANNOVAR_DIR}/table_annovar.pl deepsomatic_To_snv_avinput \
  -outfile deepsomatic_To_snv_avinput_snv \
  -buildver hg38 \
  -protocol refGene,clinvar_20240611,cosmic100coding2024 \
  -operation g,f,f \
  ${HUMANDB_DIR} \
  -otherinfo

# Step 5: Multi-tier filtering based on evidence
echo "Creating tiered variant sets..."

# Tier 1: High confidence (COSMIC + Pathogenic)
awk -F'\t' '
  /Func.refGene/ ||
  (/COSM/ && (/Pathogenic/ || /Likely_pathogenic/))
' deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_tier1_high_confidence.csv

# Tier 2: Database-supported (COSMIC OR ClinVar)
awk -F'\t' '
  /Func.refGene/ ||
  /COSM/ ||
  /Pathogenic/ ||
  /Likely_pathogenic/
' deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_tier2_database_supported.csv

# Tier 3: Quality filtered exonic variants
awk -F'\t' '
  /Func.refGene/ ||
  (/exonic/ && /nonsynonymous/ && !/Benign/) ||
  (/upstream/ && /TERT/)
' deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | awk '!/dist=166/' \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_tier3_quality_filtered.csv

# Original filtering (for comparison)
awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
  deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
  | awk '/exonic/ || /TERT/ || /Func.refGene/' \
  | awk '!/dist=166/' \
  | cut -f1-16,25,26 > ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv

# Summary report
echo "Variant filtering summary:"
echo "  Tier 1 (High confidence - COSMIC+Pathogenic): $(wc -l < ${SAMPLE_ID}_tier1_high_confidence.csv)"
echo "  Tier 2 (Database supported): $(wc -l < ${SAMPLE_ID}_tier2_database_supported.csv)"
echo "  Tier 3 (Quality filtered): $(wc -l < ${SAMPLE_ID}_tier3_quality_filtered.csv)"
echo "  Original filtering: $(wc -l < ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv)"
```

## Validation Approaches

### 1. Compare with ClairS-TO

Create a comparison script:

```bash
#!/bin/bash
# Compare DeepSomatic vs ClairS-TO variants

# Extract variant positions
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' deepsomatic.vcf.gz | sort > ds_variants.txt
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' clairsto.vcf.gz | sort > clair_variants.txt

# Find variants unique to DeepSomatic
comm -23 ds_variants.txt clair_variants.txt > deepsomatic_unique.txt

echo "DeepSomatic unique variants: $(wc -l < deepsomatic_unique.txt)"
echo "Review these for potential artifacts"
```

### 2. Visual Validation with IGV

For DeepSomatic-unique variants:

```bash
# Create IGV batch script for manual review
awk '{print $1":"$2"-"$2}' deepsomatic_unique_variants.bed > igv_review.txt
```

**Look for in IGV:**
- ❌ Variants at read ends (alignment artifacts)
- ❌ Variants in homopolymer regions
- ❌ Variants with strand bias (all on one strand)
- ❌ Variants in low-complexity regions
- ✅ Variants with balanced strand support
- ✅ Variants with consistent VAF across reads

### 3. Statistical Metrics from Benchmark

From your quality metrics comparison, flag suspicious variants:

**Red flags (likely artifacts):**
- GQ < 20 (low confidence)
- DP < 10 (low coverage)
- VAF < 5% (ultra-low frequency)
- Large depth discrepancy vs expected coverage

**Good signs (likely real):**
- GQ ≥ 30
- DP ≥ 20
- VAF 5-50% (typical somatic range)
- Presence in COSMIC/ClinVar

## Specific Recommendations for Your Data

Based on benchmark results:

### 1. Focus on the 31 Common Variants
These are **high confidence** - both callers agree:
```bash
# Extract only variants found by both callers
# Use quality_metrics_comparison.csv from your benchmark
```

### 2. Review the 35 DeepSomatic-Only Variants

From your benchmark, these genes were ONLY found by DeepSomatic:
```
AIF1L, ART1, CAMTA1, CREBBP, ERBB2, IGFLR1, KCNMB3, KMT2D,
MARCHF9, MECOM, MRPS15, MYBL1, MYH11, NOTCH4, NTRK2, NUP98,
NUTM2B, PGAP3, SGSM3, SH2D2A, SKIDA1, SMO, TENM1
```

**Prioritize for validation:**
- ✅ **ERBB2, CREBBP, KMT2D, NOTCH4, NTRK2** - Known cancer genes
- ⚠️ **NUP98** - Found in 5 samples (suspicious if recurrent)
- ⚠️ Other genes - Check COSMIC frequency

### 3. Apply Tiered Filtering

**For clinical reporting:**
1. **Report Tier 1 only** (COSMIC + Pathogenic) - Highest confidence
2. **Mention Tier 2** (Database evidence) - High confidence
3. **Flag Tier 3** for validation - Medium confidence
4. **Exclude** variants without quality metrics or database support

## Quick Implementation: Enhanced v2 Script

Would you like me to create an enhanced version of `deepsomatic_v2.sh` with these filtering strategies implemented?

It would include:
1. ✅ Quality metric filtering (GQ≥30, DP≥20, VAF≥5%)
2. ✅ PON filtering verification
3. ✅ Tiered output (Tier 1, 2, 3)
4. ✅ Comparison with ClairS-TO results
5. ✅ Summary statistics

## Summary: Answering Your Question

**How to ensure DeepSomatic variants are not artifacts?**

1. **Use multi-caller consensus:** The 31 variants both callers found = high confidence
2. **Apply quality thresholds:** GQ≥30, DP≥20, VAF≥5%
3. **Require database evidence:** Prioritize COSMIC/ClinVar variants
4. **Enable PON filtering:** Remove common germline/artifacts
5. **Manual review:** Use IGV for DeepSomatic-unique variants in known cancer genes
6. **Tier variants:** Report different confidence levels separately

**Expected outcome:**
- Current: 66 variants → Filtered: ~35-45 high-confidence variants
- This should bring DeepSomatic closer to ClairS-TO's 33 variants
- Remaining differences = DeepSomatic's better sensitivity vs artifacts

Let me know if you want me to create the enhanced filtering scripts!
