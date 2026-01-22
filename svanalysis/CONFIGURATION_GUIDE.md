# Configuration Guide for SV Analysis Script

## Overview

The `sv_analysis_v2.py` script now includes a clear **CONFIGURATION SECTION** at the top of the file (lines 28-59) that allows you to easily customize paths for different systems and directories.

## How to Configure Paths

### Step 1: Open the Script

Open [sv_analysis_v2.py](sv_analysis_v2.py) in a text editor.

### Step 2: Locate the Configuration Section

Find the configuration section at the top of the script (after the imports):

```python
# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE PATHS FOR YOUR SYSTEM
# ============================================================================

# Path to sample.txt file (contains sample_id and purity values)
# Format: sample_id\tpurity (tab-separated)
SAMPLE_FILE = "/home/chbope/extension/script/svanalysis/sample.txt"

# Directory containing VCF files (*.vcf or *.vcf.gz)
# VCF files should be named: {sample_id}.vcf or {sample_id}.vcf.gz
VCF_DIR = "/home/chbope/extension/script/svanalysis"

# Output directory for all results (tables, plots, reports)
# Will be created if it doesn't exist
OUTPUT_DIR = "/home/chbope/extension/script/svanalysis/results_v2"

# ============================================================================
# QUALITY FILTERING PARAMETERS (OPTIONAL - Advanced users only)
# ============================================================================

# Minimum read support for a variant to be considered high-quality
MIN_SUPPORT = 10

# Minimum VAF for a variant to be included in analysis
MIN_VAF = 0.1

# Clonality threshold: variants with VAF >= this are considered "clonal"
CLONAL_THRESHOLD = 0.3
```

### Step 3: Modify the Paths

Change the three main paths to match your system:

#### SAMPLE_FILE
- **Purpose**: Location of your sample.txt file
- **Format**: Tab-separated file with sample_id and purity
- **Example**:
  ```python
  SAMPLE_FILE = "/path/to/your/sample.txt"
  ```

#### VCF_DIR
- **Purpose**: Directory containing your VCF files
- **Naming**: VCF files should be named `{sample_id}.wf_sv.vcf.gz`
- **Example**:
  ```python
  VCF_DIR = "/data/gbm_project/vcf_files"
  ```

#### OUTPUT_DIR
- **Purpose**: Where all results will be saved
- **Auto-created**: Directory will be created if it doesn't exist
- **Example**:
  ```python
  OUTPUT_DIR = "/data/gbm_project/sv_analysis_results"
  ```

### Step 4: (Optional) Adjust Quality Filtering Parameters

Advanced users can modify the quality filtering thresholds:

- **MIN_SUPPORT**: Minimum number of reads supporting a variant (default: 10)
- **MIN_VAF**: Minimum variant allele frequency to include (default: 0.1 = 10%)
- **CLONAL_THRESHOLD**: VAF threshold to classify as "clonal" (default: 0.3 = 30%)

## Example Configurations

### Example 1: Different Server
```python
SAMPLE_FILE = "/mnt/data/cancer_study/samples.txt"
VCF_DIR = "/mnt/data/cancer_study/sniffles_output"
OUTPUT_DIR = "/mnt/data/cancer_study/results"
```

### Example 2: Local Analysis
```python
SAMPLE_FILE = "/home/user/project/sample.txt"
VCF_DIR = "/home/user/project/vcf"
OUTPUT_DIR = "/home/user/project/results"
```

### Example 3: Shared Cluster
```python
SAMPLE_FILE = "/scratch/lab/gbm_cohort/sample_info.txt"
VCF_DIR = "/scratch/lab/gbm_cohort/variants"
OUTPUT_DIR = "/scratch/lab/gbm_cohort/analysis/sv_results"
```

## File Requirements

### sample.txt Format
```
sample_id	purity
Sample1	0.85
Sample2	0.64
Sample3	0.50
```

- **Tab-separated** (not spaces)
- **Header line required**: `sample_id` and `purity`
- **Sample IDs**: Must match VCF filenames (without `.wf_sv.vcf.gz`)

### VCF File Naming
VCF files must be named: `{sample_id}.wf_sv.vcf.gz`

**Examples:**
- Sample ID: `Sample1` → VCF file: `Sample1.wf_sv.vcf.gz`
- Sample ID: `KM21_566` → VCF file: `KM21_566.wf_sv.vcf.gz`

## Running the Analysis

After configuring the paths, run the script:

```bash
cd /path/to/script/directory
python sv_analysis_v2.py
```

The script will:
1. Read paths from the configuration section
2. Print the paths being used
3. Process all VCF files
4. Save results to OUTPUT_DIR

## Output Files

All results will be saved in `OUTPUT_DIR`:

### Tables
- `sv_summary_display.csv` - Display table (Sample, Purity, SV counts)
- `sv_summary_table.csv` - Full table with all VAF metrics
- `clustering_results.csv` - Cluster assignments

### Plots
- `clustering_dendrogram_high_quality.png` - Clustering (recommended)
- `clustering_dendrogram_clonal.png` - Clonal clustering
- `clustering_dendrogram_weighted.png` - Weighted clustering
- `pca_plot.png` - PCA visualization
- `heatmap.png` - Feature heatmap
- `sv_distribution.png` - SV type distributions
- `sv_stacked_histogram.png` - Aggregate SV type histogram
- `vaf_analysis.png` - VAF vs purity plots
- `comprehensive_vaf_analysis.png` - 6-panel VAF comparison

### Reports
- `analysis_report.txt` - Statistical summary and cluster characteristics

## Troubleshooting

### Error: "VCF file not found"
**Problem**: Script can't find VCF files
**Solution**: Check that:
1. `VCF_DIR` points to correct directory
2. VCF filenames match pattern: `{sample_id}.wf_sv.vcf.gz`
3. Sample IDs in sample.txt match VCF filenames

### Error: "No such file or directory" (sample.txt)
**Problem**: Script can't find sample.txt
**Solution**: Verify `SAMPLE_FILE` path is correct and file exists

### Error: "Permission denied" (output directory)
**Problem**: No write permissions for output directory
**Solution**: Choose `OUTPUT_DIR` where you have write permissions

## Advanced: Adjusting Quality Thresholds

### Stricter Quality Filtering
For higher confidence (fewer variants):
```python
MIN_SUPPORT = 15      # More reads required
MIN_VAF = 0.15        # Higher VAF threshold
CLONAL_THRESHOLD = 0.4  # More stringent clonal definition
```

### More Permissive Filtering
For more sensitive detection (more variants):
```python
MIN_SUPPORT = 5       # Fewer reads required
MIN_VAF = 0.05        # Lower VAF threshold
CLONAL_THRESHOLD = 0.2  # More permissive clonal definition
```

### GBM-Optimized (Default)
Current settings are optimized for GBM analysis:
```python
MIN_SUPPORT = 10      # Balances sensitivity/specificity
MIN_VAF = 0.1         # Standard for somatic variants
CLONAL_THRESHOLD = 0.3  # Standard clonality threshold
```

## Summary

✅ All paths are now configurable in one place
✅ No need to modify code in multiple locations
✅ Easy to move script to different systems
✅ Quality thresholds can be adjusted for different cancer types
✅ Clear documentation in the configuration section

For questions about methodology, see [VAF_ANALYSIS_GUIDE.md](VAF_ANALYSIS_GUIDE.md).
