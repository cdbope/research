# Filtering Comparison: DeepSomatic vs ClairS-TO Scripts

## Side-by-Side Script Analysis

### ClairS-TO Pipeline (Your Current Script)

```bash
# 1. Variant Calling
/opt/bin/run_clairs_to \
    --tumor_bam_fn=${occ_bam} \
    --ref_fn=${reference_genome} \
    --threads=${task.cpus} \
    --platform="ont_r10_dorado_4khz" \
    --output_dir=clairsto_output \
    --bed_fn=${occ_protein_coding_bed} \
    --conda_prefix /opt/micromamba/envs/clairs-to

# 2. Merge SNV and INDEL VCFs
bcftools merge --force-samples clairsto_output/snv.vcf.gz clairsto_output/indel.vcf.gz \
    -o ${sample_id}_merge_snv_indel_claisto.vcf.gz

# 3. Convert to ANNOVAR (with PASS filter)
convert2annovar.pl ${sample_id}_merge_snv_indel_claisto.vcf.gz \
    --format vcf4 \
    --filter pass \
    --includeinfo \
    --outfile ${sample_id}_clairS_To_snv_avinput

# 4. ANNOVAR annotation
table_annovar.pl ${sample_id}_clairS_To_snv_avinput \
    -outfile ClairS_TO_snv \
    -buildver hg38 \
    -protocol refGene,clinvar_20240611,cosmic100coding2024 \
    -operation g,f,f \
    ${params.humandb_dir} \
    -otherinfo

# 5. Functional filtering (same as DeepSomatic)
awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
    ClairS_TO_snv.hg38_multianno.txt \
    | awk '/exonic/ || /TERT/ || /Func.refGene/' \
    | awk '!/dist=166/' \
    | cut -f1-16,25,26 > ${sample_id}_annotateandfilter_clairsto.csv
```

### DeepSomatic Pipeline (Your Current Script)

```bash
# 1. Variant Calling
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
    --use_default_pon_filtering=${USE_PON_FILTERING}  # ← PON filtering

# 2. Filter PASS variants
bcftools view -f PASS -O z \
    -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz \
    ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz

# 3. Convert to ANNOVAR (with PASS filter)
${ANNOVAR_DIR}/convert2annovar.pl ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz \
    --format vcf4 \
    --filter pass \
    --includeinfo \
    --outfile deepsomatic_To_snv_avinput

# 4. ANNOVAR annotation
${ANNOVAR_DIR}/table_annovar.pl deepsomatic_To_snv_avinput \
    -outfile deepsomatic_To_snv_avinput_snv \
    -buildver hg38 \
    -protocol refGene,clinvar_20240611,cosmic100coding2024 \
    -operation g,f,f \
    ${HUMANDB_DIR} \
    -otherinfo

# 5. Functional filtering (identical to ClairS-TO)
awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
    deepsomatic_To_snv_avinput_snv.hg38_multianno.txt \
    | awk '/exonic/ || /TERT/ || /Func.refGene/' \
    | awk '!/dist=166/' \
    | cut -f1-16,25,26 > ${SAMPLE_ID}_annotateandfilter_deep_somatic.csv
```

## Key Differences in Filtering

### ✅ IDENTICAL FILTERING (Both Scripts)

| Filter Stage | Filter Applied | Effect |
|-------------|----------------|--------|
| **VCF Level** | `--filter pass` in convert2annovar | Only PASS variants |
| **Annotation** | Same ANNOVAR databases | Same annotation |
| **Functional** | `/exonic/ && /nonsynonymous/ && !/Benign/` | Exonic nonsynonymous, not benign |
| **Functional** | `/upstream/ \|\| /TERT/` | TERT promoter variants |
| **Functional** | `!/dist=166/` | Exclude specific distance |

### ⚠️ DIFFERENT FILTERING

| Feature | DeepSomatic | ClairS-TO | Impact |
|---------|-------------|-----------|--------|
| **PON Filtering** | `--use_default_pon_filtering=true` | ❌ NOT VISIBLE | May filter germline/artifacts |
| **BED Region** | ❌ None | `--bed_fn=${occ_protein_coding_bed}` | ClairS-TO only calls in BED regions |
| **SNV+INDEL Merge** | Single VCF output | Separate, then merged | ClairS-TO processes separately |
| **Quality Thresholds** | ❌ None (just PASS) | ❌ None (just PASS) | **BOTH lack GQ/DP/VAF filters** |

## Critical Finding

### 🚨 BOTH PIPELINES LACK QUALITY METRIC FILTERING

Neither your DeepSomatic nor ClairS-TO script applies:
- ❌ Genotype Quality (GQ) thresholds
- ❌ Read Depth (DP) thresholds
- ❌ Allele Frequency (VAF) thresholds

**Both rely ONLY on:**
1. VCF PASS filter (from variant caller)
2. Functional annotation filtering (exonic, nonsynonymous, etc.)

## Why DeepSomatic Has More Variants

### Factor 1: BED Region Filtering ⭐ MAJOR DIFFERENCE

**ClairS-TO:**
```bash
--bed_fn=${occ_protein_coding_bed}  # Only calls variants in specified regions
```

**DeepSomatic:**
```bash
# NO BED file specified → Calls variants genome-wide (or at least more broadly)
```

**Impact:** If ClairS-TO BED file only contains specific genes/regions, it will miss variants outside those regions.

**Example:**
- BED contains 500 cancer genes
- DeepSomatic calls genome-wide
- DeepSomatic will find variants in genes NOT in the BED file

### Factor 2: PON Filtering (Unclear)

**DeepSomatic:**
```bash
--use_default_pon_filtering=${USE_PON_FILTERING}
```

**Question:** What is `${USE_PON_FILTERING}` set to in your script?
- If `true`: Removes germline/artifacts → FEWER variants
- If `false`: Keeps everything → MORE variants

**ClairS-TO:** No visible PON filtering in your script

### Factor 3: Variant Caller Sensitivity

**DeepSomatic:**
- Deep learning model
- May be more sensitive to low VAF variants
- Benchmark shows min VAF: 5.36%

**ClairS-TO:**
- Also deep learning, but different architecture
- Benchmark shows min VAF: 13.04%
- Less sensitive to very low VAF variants

### Factor 4: PASS Filter Stringency

**DeepSomatic PASS criteria:** Unknown (internal to tool)
**ClairS-TO PASS criteria:** Unknown (internal to tool)

These may differ significantly!

## Recommendations for Equal Filtering

### Option 1: Add BED File to DeepSomatic (Match ClairS-TO)

Make DeepSomatic call only in the same regions as ClairS-TO:

```bash
# Add to DeepSomatic
google/deepsomatic:"${BIN_VERSION}" \
run_deepsomatic \
    --regions=${BED_FILE} \  # ← ADD THIS
    ...other parameters...
```

This will likely reduce DeepSomatic variants to be closer to ClairS-TO.

### Option 2: Add Quality Filters to BOTH Pipelines (Recommended)

Apply the same quality thresholds to both:

**For DeepSomatic (Enhanced version):**
```bash
bcftools view -f PASS ${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz | \
bcftools filter -i 'GQ>=30 && INFO/DP>=20 && FORMAT/VAF>=0.05' \
    -O z -o ${SAMPLE_ID}_filtered.vcf.gz
```

**For ClairS-TO (NEW - add to your pipeline):**
```bash
# After bcftools merge, before convert2annovar
bcftools filter -i 'GQ>=30 && INFO/DP>=20 && FORMAT/AF>=0.05' \
    ${sample_id}_merge_snv_indel_claisto.vcf.gz \
    -O z -o ${sample_id}_merge_snv_indel_claisto.filtered.vcf.gz

# Then use the filtered VCF for convert2annovar
convert2annovar.pl ${sample_id}_merge_snv_indel_claisto.filtered.vcf.gz \
    --format vcf4 \
    --filter pass \
    --includeinfo \
    --outfile ${sample_id}_clairS_To_snv_avinput
```

### Option 3: Check ClairS-TO Internal Filtering

ClairS-TO may have built-in quality filtering. Check:

```bash
# Look at ClairS-TO output VCF header for filters
bcftools view -h clairsto_output/snv.vcf.gz | grep "##FILTER"

# Check what variants are marked as PASS vs filtered
bcftools view -H clairsto_output/snv.vcf.gz | cut -f7 | sort | uniq -c
```

## Enhanced ClairS-TO Script with Quality Filtering

Here's your ClairS-TO script with added quality filtering to match DeepSomatic enhanced:

```bash
#!/bin/bash

# 1. Run ClairS-TO
/opt/bin/run_clairs_to \
    --tumor_bam_fn=${occ_bam} \
    --ref_fn=${reference_genome} \
    --threads=${task.cpus} \
    --platform="ont_r10_dorado_4khz" \
    --output_dir=clairsto_output \
    --bed_fn=${occ_protein_coding_bed} \
    --conda_prefix /opt/micromamba/envs/clairs-to

# 2. Merge SNV and INDEL
bcftools merge --force-samples \
    clairsto_output/snv.vcf.gz \
    clairsto_output/indel.vcf.gz \
    -O z -o ${sample_id}_merge_snv_indel_claisto.vcf.gz

# 3. NEW: Apply quality filtering
bcftools filter -i 'FILTER="PASS" && GQ>=30 && DP>=20 && AF>=0.05' \
    ${sample_id}_merge_snv_indel_claisto.vcf.gz \
    -O z -o ${sample_id}_merge_snv_indel_claisto.qc.vcf.gz

bcftools index -t ${sample_id}_merge_snv_indel_claisto.qc.vcf.gz

# 4. Convert to ANNOVAR (using QC-filtered VCF)
convert2annovar.pl ${sample_id}_merge_snv_indel_claisto.qc.vcf.gz \
    --format vcf4 \
    --filter pass \
    --includeinfo \
    --outfile ${sample_id}_clairS_To_snv_avinput

# 5. ANNOVAR annotation
table_annovar.pl ${sample_id}_clairS_To_snv_avinput \
    -outfile ClairS_TO_snv \
    -buildver hg38 \
    -protocol refGene,clinvar_20240611,cosmic100coding2024 \
    -operation g,f,f \
    ${params.humandb_dir} \
    -otherinfo

# 6. Tiered filtering (like DeepSomatic enhanced)

# Tier 1: High confidence
awk -F'\t' '
    /^Chr\t/ ||
    (/COSM/ && (/Pathogenic/ || /Likely_pathogenic/))
' ClairS_TO_snv.hg38_multianno.txt \
    | cut -f1-16,25,26 > ${sample_id}_tier1_high_confidence_clairsto.csv

# Tier 2: Database evidence
awk -F'\t' '
    /^Chr\t/ ||
    /COSM/ ||
    /Pathogenic/ ||
    /Likely_pathogenic/
' ClairS_TO_snv.hg38_multianno.txt \
    | cut -f1-16,25,26 > ${sample_id}_tier2_database_evidence_clairsto.csv

# Tier 3: Functional variants (original filtering)
awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
    ClairS_TO_snv.hg38_multianno.txt \
    | awk '/exonic/ || /TERT/ || /Func.refGene/' \
    | awk '!/dist=166/' \
    | cut -f1-16,25,26 > ${sample_id}_annotateandfilter_clairsto.csv
```

## Summary: Answer to Your Question

**Q: Does ClairS-TO do the same filtering as DeepSomatic?**

**A: NO, but the differences are subtle:**

| Filtering Stage | DeepSomatic | ClairS-TO | Same? |
|----------------|-------------|-----------|-------|
| **PON Filtering** | Yes (if enabled) | Not visible | ❌ |
| **BED Region Limiting** | No | Yes | ❌ MAJOR |
| **PASS Filter** | Yes | Yes | ✅ |
| **Quality Metrics (GQ/DP/VAF)** | ❌ No (in current script) | ❌ No | ✅ Both lacking |
| **Functional Annotation** | Yes | Yes | ✅ Identical |
| **Database Annotation** | Same databases | Same databases | ✅ |

**The BED file is likely the MAIN reason ClairS-TO has fewer variants.**

If `${occ_protein_coding_bed}` only contains specific genes/regions, ClairS-TO will NEVER call variants outside those regions, while DeepSomatic calls genome-wide.

## Action Items

1. **Check BED File Coverage:**
   ```bash
   wc -l ${occ_protein_coding_bed}
   # How many regions/genes?
   ```

2. **Add Same BED to DeepSomatic:**
   ```bash
   --regions=${BED_FILE}
   ```

3. **Add Quality Filtering to Both:**
   - Use enhanced scripts provided
   - Apply same GQ/DP/VAF thresholds to both callers

4. **Compare Apples-to-Apples:**
   - Same regions (BED file)
   - Same quality thresholds
   - Then see true sensitivity differences
