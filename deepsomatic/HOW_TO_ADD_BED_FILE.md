# How to Add BED File to DeepSomatic

## Quick Answer

Add this line to your DeepSomatic `run_deepsomatic` command:

```bash
--regions="${BED_FILE}" \
```

## Step-by-Step Instructions

### Method 1: Modify deepsomatic_v2.sh (Recommended)

#### 1. Edit the config.sh file

Add BED file path to your config:

```bash
# Add this line to config.sh
BED_FILE="/path/to/your/bedfile.bed"
```

#### 2. Add volume mount for BED file

In `deepsomatic_v2.sh` line 46-49, add BED directory mount:

**Before:**
```bash
sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}":ro \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${REF_DIR}":"${REF_DIR}":ro \
```

**After:**
```bash
# Get BED directory
BED_DIR=$(dirname "${BED_FILE}")

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}":ro \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${REF_DIR}":"${REF_DIR}":ro \
-v "${BED_DIR}":"${BED_DIR}":ro \
```

#### 3. Add --regions parameter

In `deepsomatic_v2.sh` around line 51-60, add the regions parameter:

**Before:**
```bash
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
```

**After:**
```bash
google/deepsomatic:"${BIN_VERSION}" \
run_deepsomatic \
--model_type=${MODEL_TYPE} \
--ref="${REF_DIR}/${REF_GENOME}" \
--reads_tumor="${INPUT_DIR}/${BAM_FILE}" \
--regions="${BED_FILE}" \
--output_vcf="${OUTPUT_DIR}/${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz" \
--sample_name_tumor="${SAMPLE_ID}" \
--num_shards=${NUM_SHARDS} \
--logging_dir="${OUTPUT_DIR}/logs" \
--intermediate_results_dir="${OUTPUT_DIR}/intermediate_results_dir" \
--use_default_pon_filtering=${USE_PON_FILTERING}
```

### Method 2: Create New Script with BED Support

I can create a new script `deepsomatic_v2_with_bed.sh` that includes BED file support.

### Method 3: Use the Same BED as ClairS-TO

If you want to match your ClairS-TO analysis exactly:

```bash
# In config.sh, add:
BED_FILE="${occ_protein_coding_bed}"  # Use same BED as ClairS-TO
```

## Complete Modified Script Section

Here's the exact change needed in deepsomatic_v2.sh:

```bash
# After line 43, add:
# Get BED directory if BED_FILE is set
if [ -n "${BED_FILE}" ] && [ -f "${BED_FILE}" ]; then
    BED_DIR=$(dirname "${BED_FILE}")
    echo "Using BED file for region filtering: ${BED_FILE}"
    BED_MOUNT="-v ${BED_DIR}:${BED_DIR}:ro"
    REGIONS_PARAM="--regions=${BED_FILE}"
else
    echo "No BED file specified. Calling variants genome-wide."
    BED_MOUNT=""
    REGIONS_PARAM=""
fi

# Then modify docker run command (line 46):
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
```

## BED File Format

Your BED file should be in standard format:

```
chr1    1000000    1001000    GENE1
chr1    2000000    2001000    GENE2
chr2    3000000    3001000    GENE3
```

Or minimal format (3 columns):
```
chr1    1000000    1001000
chr1    2000000    2001000
chr2    3000000    3001000
```

## Verify BED File

Before running, check your BED file:

```bash
# Check format
head ${BED_FILE}

# Count regions
wc -l ${BED_FILE}

# Check if it's sorted (required)
sort -k1,1 -k2,2n ${BED_FILE} | diff - ${BED_FILE}
# If output is empty, file is sorted
# If output shows differences, sort it:
sort -k1,1 -k2,2n ${BED_FILE} > ${BED_FILE}.sorted
```

## Testing

After modification, test with a small sample:

```bash
# Set BED file in config.sh
echo 'BED_FILE="/path/to/your/regions.bed"' >> config.sh

# Run the script
./deepsomatic_v2.sh

# Check output - should have fewer variants if BED is restrictive
```

## Expected Results

**Without BED file:**
- DeepSomatic calls variants genome-wide
- ~66 variants (from your benchmark)

**With BED file (same as ClairS-TO):**
- DeepSomatic calls only in BED regions
- Expected: ~35-45 variants (closer to ClairS-TO's 33)

## Troubleshooting

### Error: "Could not open BED file"

```bash
# Check file exists
ls -lh ${BED_FILE}

# Check Docker can access it (volume mount issue)
# Make sure BED_DIR is mounted with -v flag
```

### Error: "BED file not sorted"

```bash
# Sort the BED file
sort -k1,1 -k2,2n ${BED_FILE} > ${BED_FILE}.sorted
mv ${BED_FILE}.sorted ${BED_FILE}
```

### Error: "Chromosome naming mismatch"

```bash
# Check chromosome naming
head ${BED_FILE}
# Should match reference: chr1, chr2... or 1, 2...

# Check reference chromosome names
samtools faidx ${REF_GENOME}
head ${REF_GENOME}.fai
```

## DeepSomatic BED Documentation

According to DeepSomatic documentation, the `--regions` parameter accepts:

1. **BED file:** `--regions=/path/to/file.bed`
2. **Region string:** `--regions="chr1:1000-2000"`
3. **Multiple regions:** `--regions="chr1:1000-2000 chr2:3000-4000"`

## Alternative: Post-Calling BED Filtering

If you don't want to modify the calling step, you can filter the VCF after calling:

```bash
# After DeepSomatic completes, filter VCF to BED regions
bcftools view -R ${BED_FILE} \
    ${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz \
    -O z -o ${SAMPLE_ID}_ont_deepsomatic_output.pass.bed_filtered.vcf.gz

# Then proceed with ANNOVAR annotation
```

**Pros:** Simple, no Docker changes
**Cons:** Wastes compute time calling variants outside regions

## Recommendation

**For fair comparison with ClairS-TO:**
1. Use the SAME BED file as ClairS-TO
2. Modify deepsomatic_v2.sh as shown above
3. Re-run benchmark comparison
4. Compare apples-to-apples

This will show you the TRUE difference in caller sensitivity, not just region coverage differences.
