# DeepSomatic Pipeline v3 - Batch Processing

This version processes multiple samples automatically from a sample list file with ANNOVAR annotation.

## Differences Between Versions

| Feature | v2 | v3 |
|---------|----|----|
| Sample Mode | Single sample | Multiple samples (batch) |
| Config File | config.sh | config_v3.sh |
| Sample Input | Hardcoded in config | Text file with list |
| Output Structure | Single directory | Per-sample subdirectories |
| Progress Tracking | No | Yes (with summary) |
| Error Handling | Stops on error | Continues with remaining samples |

## Quick Start

### 1. Create Sample List File

Create a file `sample_id.txt` in your INPUT_DIR with one sample ID per line:

```bash
cd /home/chbope/extension/nWGS_manuscript_data/data/testdata/single_bam_folder/merge_bam
nano sample_id.txt
```

Example content:
```
T25-152
T001
T002
```

### 2. Verify BAM Files Exist

The pipeline will look for BAM files with the naming pattern: `{SAMPLE_ID}{BAM_SUFFIX}`

For example, with `BAM_SUFFIX=".occ.bam"`:
- T25-152 → T25-152.occ.bam
- T001 → T001.occ.bam

Verify your BAM files:
```bash
cd /home/chbope/extension/nWGS_manuscript_data/data/testdata/single_bam_folder/merge_bam
ls *.bam
```

### 3. Configure the Pipeline

Edit [config_v3.sh](config_v3.sh) if needed (already pre-configured with your paths):

```bash
cd /home/chbope/extension/script/deepsomatic
nano config_v3.sh
```

Key settings:
- `SAMPLE_LIST_FILE` - Path to your sample_id.txt file
- `BAM_SUFFIX` - Suffix to append to sample IDs (default: `.occ.bam`)
- `MODEL_TYPE` - DeepSomatic model (default: `ONT_TUMOR_ONLY`)

### 4. Run the Pipeline

```bash
cd /home/chbope/extension/script/deepsomatic
bash deepsomatic_v3.sh
```

## Output Structure

For each sample, outputs are organized in separate subdirectories:

```
OUTPUT_DIR/
├── T25-152/
│   ├── T25-152_ont_deepsomatic_output.vcf.gz          # Raw VCF
│   ├── T25-152_ont_deepsomatic_output.pass.vcf.gz     # PASS-filtered VCF
│   ├── T25-152_deepsomatic_to_snv_avinput             # ANNOVAR input
│   ├── T25-152_deepsomatic_annotated.hg38_multianno.txt  # Full annotation
│   ├── T25-152_annotateandfilter_deep_somatic.csv     # Final filtered variants
│   ├── logs/                                           # DeepSomatic logs
│   └── intermediate_results_dir/                       # Temporary files
├── T001/
│   ├── T001_ont_deepsomatic_output.vcf.gz
│   ├── ...
└── T002/
    └── ...
```

## Pipeline Steps for Each Sample

1. **Variant Calling**: DeepSomatic somatic variant detection
2. **Filtering**: Extract PASS-quality variants with bcftools
3. **Conversion**: Convert VCF to ANNOVAR format
4. **Annotation**: Annotate with refGene, ClinVar, COSMIC
5. **Post-filtering**: Filter for clinically relevant variants

## Features

### Progress Tracking

The pipeline displays:
- Current sample number and total (e.g., "Processing sample 2/10")
- Step-by-step progress for each sample
- Timestamp for each sample start/completion
- Real-time error notifications

### Error Handling

- Validates BAM file existence before processing
- Checks success at each step
- Continues with remaining samples if one fails
- Provides detailed error messages
- Final summary shows all failed samples

### Summary Report

At completion, displays:
```
========================================
BATCH PROCESSING SUMMARY
========================================
Total samples: 3
Successfully processed: 2
Failed: 1

Failed samples:
  - T002 (BAM not found)

All results saved in: /path/to/output
========================================
```

## Sample List File Format

The sample list file supports:

- **One sample ID per line**: Simple text format
- **Comments**: Lines starting with `#` are ignored
- **Empty lines**: Blank lines are skipped
- **No header**: Start directly with sample IDs

Example:
```
# Batch 1 - Oxford Nanopore samples
T25-152
T001

# Batch 2 - FFPE samples
T002
T003
```

## Customizing BAM File Names

The pipeline constructs BAM filenames as: `{SAMPLE_ID}{BAM_SUFFIX}`

**Examples:**

| BAM_SUFFIX | Sample ID | Resulting BAM Filename |
|------------|-----------|------------------------|
| `.occ.bam` | T25-152 | T25-152.occ.bam |
| `.bam` | T001 | T001.bam |
| `_sorted.bam` | Sample1 | Sample1_sorted.bam |

Edit `BAM_SUFFIX` in [config_v3.sh](config_v3.sh) to match your naming convention.

## Advanced Usage

### Processing Specific Samples

Create a temporary sample list:
```bash
echo "T25-152" > temp_samples.txt
echo "T001" >> temp_samples.txt
```

Update `SAMPLE_LIST_FILE` in config_v3.sh:
```bash
SAMPLE_LIST_FILE="temp_samples.txt"
```

### Parallel Processing

Currently, samples are processed sequentially. Each sample uses `NUM_SHARDS` (default: 12) parallel threads internally.

### Resume After Failure

The pipeline continues processing remaining samples after a failure. To reprocess only failed samples:

1. Check the summary for failed samples
2. Create a new sample list with only failed samples
3. Update `SAMPLE_LIST_FILE` and rerun

## Troubleshooting

### "BAM file not found" Error

**Cause**: BAM filename doesn't match expected pattern

**Solution**:
1. Check actual BAM filenames: `ls ${INPUT_DIR}/*.bam`
2. Adjust `BAM_SUFFIX` in config_v3.sh to match

Example:
```bash
# If your files are named T25-152.bam instead of T25-152.occ.bam
BAM_SUFFIX=".bam"
```

### "Sample list file not found" Error

**Cause**: `SAMPLE_LIST_FILE` path is incorrect

**Solution**:
```bash
# Verify the file exists
ls -l ${INPUT_DIR}/sample_id.txt

# Or use absolute path in config_v3.sh
SAMPLE_LIST_FILE="/absolute/path/to/sample_id.txt"
```

### No Variants in Final Output

**Cause**: All variants filtered out by stringent criteria

**Solution**: Check intermediate files:
```bash
# Check raw variant count
zcat ${OUTPUT_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz | grep -v "^#" | wc -l

# Check PASS variants
zcat ${OUTPUT_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz | grep -v "^#" | wc -l

# Check annotated variants
wc -l ${OUTPUT_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_deepsomatic_annotated.hg38_multianno.txt
```

### Pipeline Stops on First Sample

**Cause**: v2 script being used instead of v3

**Solution**: Ensure you're running the correct script:
```bash
bash deepsomatic_v3.sh  # NOT deepsomatic_v2.sh
```

## Performance Tips

1. **Use SSD storage** for OUTPUT_DIR to speed up I/O
2. **Adjust NUM_SHARDS** based on CPU cores and memory
3. **Process samples in batches** if you have many samples
4. **Monitor disk space** - intermediate files can be large

## Comparison with v2

**Use v2 when:**
- Processing a single sample
- Need immediate feedback for one sample
- Testing pipeline settings

**Use v3 when:**
- Processing multiple samples
- Need batch processing automation
- Want per-sample error isolation
- Require processing summary

## Files

- `deepsomatic_v3.sh` - Main batch processing script
- `config_v3.sh` - Configuration file for v3
- `sample_id.txt.example` - Example sample list file
- `README_v3.md` - This documentation

## Example Workflow

```bash
# 1. Create sample list
cd /home/chbope/extension/nWGS_manuscript_data/data/testdata/single_bam_folder/merge_bam
cat > sample_id.txt << EOF
T25-152
T001
EOF

# 2. Verify configuration
cd /home/chbope/extension/script/deepsomatic
cat config_v3.sh

# 3. Run pipeline
bash deepsomatic_v3.sh

# 4. Check results
ls -lh /home/chbope/extension/data/200GMBs/pcgr/starsigndna/refit/tmp/*/
```

## Support

For issues specific to:
- **DeepSomatic**: https://github.com/google/deepvariant
- **ANNOVAR**: http://annovar.openbioinformatics.org/
- **This Pipeline**: Check logs in `${OUTPUT_DIR}/${SAMPLE_ID}/logs/`
