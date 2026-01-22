# CNV Pipeline Troubleshooting Guide

## Common Issues and Solutions

### Issue 1: No CNV Data Extracted - Empty VCF Files

**Error Message:**
```
Warning message:
No VCF found in /path/to/directory
Error: No CNV data could be extracted from any time point
```

**Root Cause:**
The BAM files have insufficient read coverage for CNV detection. CNVpytor requires adequate genomic coverage (typically 10-30x) to reliably detect copy number variations.

**Diagnosis:**
Check your BAM file coverage:
```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate cnvgen
samtools flagstat your_file.bam
samtools idxstats your_file.bam | head -20
```

If you see only a few thousand reads (e.g., < 100,000 reads), this is insufficient for CNV calling.

**Solutions:**

#### Option 1: Use Higher Coverage Data
- Process BAM files with adequate coverage (millions of reads)
- Ensure your sequencing depth is at least 5-10x genome-wide
- For cancer samples, 30x or higher is recommended

#### Option 2: Generate Synthetic Test Data
For testing the visualization pipeline with synthetic data:

```bash
cd /home/chbope/extension/script/cnvmonitor

# Generate synthetic test data
source ~/miniconda3/etc/profile.d/conda.sh
conda activate cnvgen

Rscript generate_test_cnv_data.R \
  --output-dir /path/to/output/synthetic_test \
  --sample-id YOUR_SAMPLE_ID \
  --bin-size 100000 \
  --num-timepoints 5

# Run visualization
Rscript ont_timeseries_cnv_visualization.R \
  --input-dir /path/to/output/synthetic_test \
  --output-dir /path/to/output/synthetic_test/visualization \
  --sample-id YOUR_SAMPLE_ID \
  --bin-size 100000
```

#### Option 3: Adjust CNVpytor Parameters
For low-coverage data, try using larger bin sizes:

```bash
./ont_process_individual_bam_cnv.sh \
  --bam your_file.bam \
  --output-dir ./output \
  --reference reference.fa \
  --sample-id SAMPLE_ID \
  --bin-sizes "100000 1000000" \
  --threads 4
```

### Issue 2: VCF File Naming Pattern Not Recognized

**Error Message:**
```
Warning message:
No VCF found in /path/to/directory
```

**Root Cause:**
The visualization script expects VCF files named with the pattern: `*_<binsize>bp.vcf`

For example:
- `T001_100000bp.vcf` ✓ Correct
- `T001.vcf` ✗ Wrong (missing bin size)
- `T001.vcf.gz` ✗ Wrong (compressed, and missing bin size)

**Solution:**
Ensure your VCF files follow the CNVpytor naming convention. The [ont_process_individual_bam_cnv.sh](ont_process_individual_bam_cnv.sh) script automatically generates properly named files.

### Issue 3: Wrong VCF Format (Sniffles2 vs CNVpytor)

**Error Message:**
May fail silently or report "No CNV data extracted"

**Root Cause:**
The visualization script expects CNVpytor VCF format, not Sniffles2 or other SV caller formats.

**Solution:**
Use the CNVpytor processing pipeline:
```bash
./ont_process_individual_bam_cnv.sh \
  --bam your_file.bam \
  --output-dir ./output \
  --reference reference.fa \
  --sample-id SAMPLE_ID
```

### Issue 4: Missing Conda Environment

**Error Message:**
```
conda: command not found
```

**Solution:**
Activate conda and create the environment:
```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda env create -f cnvgen_environment.yml
conda activate cnvgen
```

### Issue 5: Missing R Packages

**Error Message:**
```
Error in library(packagename) : there is no package called 'packagename'
```

**Solution:**
Install missing R packages in the cnvgen environment:
```bash
conda activate cnvgen
Rscript -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'data.table', 'cowplot', 'viridis', 'scales', 'gridExtra', 'ggpubr', 'optparse'), repos='https://cloud.r-project.org/')"
```

## Expected Directory Structure

For the visualization script to work properly, your directory structure should look like:

```
input_directory/
├── SAMPLE_timepoint_1/
│   ├── SAMPLE_timepoint_1_100000bp.vcf
│   └── SAMPLE_timepoint_1_metadata.json
├── SAMPLE_timepoint_2/
│   ├── SAMPLE_timepoint_2_100000bp.vcf
│   └── SAMPLE_timepoint_2_metadata.json
└── SAMPLE_timepoint_3/
    ├── SAMPLE_timepoint_3_100000bp.vcf
    └── SAMPLE_timepoint_3_metadata.json
```

Each subdirectory must:
1. Contain the sample ID in its name
2. Have a VCF file with the pattern `*_<binsize>bp.vcf`
3. Optionally have a metadata JSON file with timestamp

## Testing Your Installation

### Quick Test with Synthetic Data

```bash
cd /home/chbope/extension/script/cnvmonitor

# Generate test data
source ~/miniconda3/etc/profile.d/conda.sh
conda activate cnvgen

Rscript generate_test_cnv_data.R \
  --output-dir ./test_output \
  --sample-id TEST \
  --bin-size 100000 \
  --num-timepoints 3

# Run visualization
Rscript ont_timeseries_cnv_visualization.R \
  --input-dir ./test_output \
  --output-dir ./test_output/visualization \
  --sample-id TEST \
  --bin-size 100000

# Check results
ls -lh ./test_output/visualization/
```

Expected output files:
- `TEST_cnv_count_timeseries.png` - CNV counts over time
- `TEST_genome_wide_cnv_heatmap.png` - Genome-wide heatmap
- `TEST_per_chromosome_cnv.png` - Per-chromosome profiles
- `TEST_cnv_size_distribution.png` - Size distributions
- `TEST_timeseries_summary.txt` - Summary statistics

## Getting Help

If you continue to experience issues:

1. Check the log files in the output directory
2. Verify your BAM file quality with `samtools flagstat`
3. Ensure your reference genome matches your BAM file
4. Test with synthetic data to isolate the issue
5. Review the script documentation in the README

## Coverage Requirements Summary

| Analysis Type | Minimum Coverage | Recommended Coverage |
|--------------|------------------|---------------------|
| Germline CNV | 10x | 30x |
| Somatic CNV | 30x | 60x+ |
| Low-confidence regions | 5x | 10x |

**Note:** Coverage is genome-wide average. Local coverage may vary.
