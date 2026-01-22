# DeepSomatic Pipeline Version Comparison

## Quick Reference

| Feature | v1 (Original) | v2 (Single Sample + ANNOVAR) | v3 (Batch + ANNOVAR) |
|---------|---------------|------------------------------|----------------------|
| **Purpose** | Basic variant calling | Single sample with annotation | Multiple samples with annotation |
| **Script** | deepsomatic.sh | deepsomatic_v2.sh | deepsomatic_v3.sh |
| **Config** | Hardcoded in script | config.sh | config_v3.sh |
| **Sample Input** | Hardcoded | Hardcoded in config | Text file (sample_id.txt) |
| **ANNOVAR** | ❌ No | ✅ Yes | ✅ Yes |
| **Filtering** | ❌ No | ✅ Yes (bcftools + custom) | ✅ Yes (bcftools + custom) |
| **Output Structure** | Single directory | Single directory | Per-sample subdirectories |
| **Progress Tracking** | ❌ No | ❌ No | ✅ Yes |
| **Error Handling** | Stops on error | Stops on error | Continues processing |
| **Summary Report** | ❌ No | ❌ No | ✅ Yes |
| **Reference Mount Fix** | ✅ Yes | ✅ Yes | ✅ Yes |
| **Cleanup** | ✅ Yes | ✅ Yes | ✅ Yes (per-sample) |

## Which Version Should I Use?

### Use v1 (deepsomatic.sh) if:
- You only need basic variant calling
- You don't need annotation
- You want the simplest setup
- You'll handle annotation separately

### Use v2 (deepsomatic_v2.sh) if:
- Processing a **single sample**
- Need ANNOVAR annotation
- Want filtered, clinically relevant variants
- Testing pipeline on one sample first

### Use v3 (deepsomatic_v3.sh) if:
- Processing **multiple samples** (2+)
- Need batch processing automation
- Want separate output directories per sample
- Need progress tracking and summary reports
- One sample failure shouldn't stop others

## Example Use Cases

### Scenario 1: Quick Test
**Use v2**
```bash
# Edit config.sh with your sample
SAMPLE_ID="test_sample"
BAM_FILE="test_sample.bam"

# Run
bash deepsomatic_v2.sh
```

### Scenario 2: Single Patient Analysis
**Use v2**
```bash
# For one patient with detailed annotation
SAMPLE_ID="Patient001"
bash deepsomatic_v2.sh
```

### Scenario 3: Cohort Study (10+ samples)
**Use v3**
```bash
# Create sample list
cat > sample_id.txt << EOF
Patient001
Patient002
Patient003
...
Patient100
EOF

# Run batch
bash deepsomatic_v3.sh
```

### Scenario 4: Production Pipeline
**Use v3**
- Automated processing
- Multiple samples
- Need error resilience
- Want processing logs per sample

## Output Comparison

### v2 Output Structure
```
OUTPUT_DIR/
├── T25-152_ont_deepsomatic_output.vcf.gz
├── T25-152_ont_deepsomatic_output.pass.vcf.gz
├── deepsomatic_to_snv_avinput
├── deepsomatic_annotated.hg38_multianno.txt
├── T25-152_annotateandfilter_deep_somatic.csv
├── logs/
└── intermediate_results_dir/
```

### v3 Output Structure
```
OUTPUT_DIR/
├── T25-152/
│   ├── T25-152_ont_deepsomatic_output.vcf.gz
│   ├── T25-152_ont_deepsomatic_output.pass.vcf.gz
│   ├── T25-152_deepsomatic_to_snv_avinput
│   ├── T25-152_deepsomatic_annotated.hg38_multianno.txt
│   ├── T25-152_annotateandfilter_deep_somatic.csv
│   ├── logs/
│   └── intermediate_results_dir/
├── T001/
│   └── [same structure]
└── T002/
    └── [same structure]
```

## Migration Guide

### From v1 to v2
1. Install bcftools and ANNOVAR
2. Create config.sh with your paths
3. Run deepsomatic_v2.sh instead of deepsomatic.sh

### From v2 to v3
1. Create sample_id.txt file with your sample IDs
2. Copy config.sh to config_v3.sh
3. Update SAMPLE_LIST_FILE in config_v3.sh
4. Run deepsomatic_v3.sh

### From v1 to v3
1. Install bcftools and ANNOVAR
2. Create config_v3.sh with your paths
3. Create sample_id.txt with sample IDs
4. Run deepsomatic_v3.sh

## Common Configurations

### ONT (Oxford Nanopore) - Tumor Only
```bash
# All versions
MODEL_TYPE="ONT_TUMOR_ONLY"
BAM_SUFFIX=".occ.bam"  # or your naming convention
```

### FFPE WGS - Tumor Only
```bash
MODEL_TYPE="FFPE_WGS_TUMOR_ONLY"
BAM_SUFFIX=".bam"
```

### FFPE WES - Tumor Only
```bash
MODEL_TYPE="FFPE_WES_TUMOR_ONLY"
BAM_SUFFIX=".bam"
```

## Performance Comparison

| Aspect | v2 | v3 |
|--------|----|----|
| **Setup Time** | 5 min | 10 min (first time) |
| **Per-Sample Time** | ~30-60 min | ~30-60 min |
| **Total Time (10 samples)** | 5-10 hours (manual) | 5-10 hours (automatic) |
| **Manual Intervention** | High (per sample) | Low (only setup) |
| **Error Recovery** | Manual restart | Automatic continue |

## File Locations

```
/home/chbope/extension/script/deepsomatic/
├── deepsomatic.sh           # v1 - Basic variant calling
├── deepsomatic_v2.sh        # v2 - Single sample + ANNOVAR
├── deepsomatic_v3.sh        # v3 - Batch + ANNOVAR
├── config.sh                # v2 configuration
├── config_v3.sh             # v3 configuration
├── sample_id.txt.example    # Example sample list for v3
├── README_v2.md             # v2 documentation
├── README_v3.md             # v3 documentation
└── VERSION_COMPARISON.md    # This file
```

## Getting Help

1. **v2 Issues**: See [README_v2.md](README_v2.md)
2. **v3 Issues**: See [README_v3.md](README_v3.md)
3. **General Issues**: Check logs in respective output directories
