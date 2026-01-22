# Portability Update - Summary

## What Was Changed

The `sv_analysis_v2.py` script has been updated to make it fully portable and easy to configure for different systems and directories.

## Key Changes

### 1. Configuration Section Added (Lines 28-59)

A clear **CONFIGURATION SECTION** has been added at the top of the script with all configurable parameters:

```python
# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE PATHS FOR YOUR SYSTEM
# ============================================================================

# Path to sample.txt file (contains sample_id and purity values)
SAMPLE_FILE = "/home/chbope/extension/script/svanalysis/sample.txt"

# Directory containing VCF files (*.vcf or *.vcf.gz)
VCF_DIR = "/home/chbope/extension/script/svanalysis"

# Output directory for all results (tables, plots, reports)
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

### 2. Script Updated to Use Configuration

- **main() function**: Now uses `SAMPLE_FILE`, `VCF_DIR`, `OUTPUT_DIR` instead of hardcoded paths
- **analyze_all_samples()**: Uses `MIN_SUPPORT` and `MIN_VAF` parameters
- **parse_vcf()**: Uses `CLONAL_THRESHOLD` for clonality classification

### 3. Documentation Created

New file: [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md)
- Complete instructions for configuring paths
- Example configurations for different systems
- Troubleshooting guide
- Advanced parameter tuning guide

## How to Use

### Quick Start

1. **Open the script**:
   ```bash
   nano sv_analysis_v2.py
   # or use your preferred editor
   ```

2. **Find the CONFIGURATION SECTION** (lines 28-59)

3. **Update the three main paths**:
   ```python
   SAMPLE_FILE = "/your/path/to/sample.txt"
   VCF_DIR = "/your/path/to/vcf/files"
   OUTPUT_DIR = "/your/path/to/output/directory"
   ```

4. **Run the analysis**:
   ```bash
   python sv_analysis_v2.py
   ```

### Example: Moving to a Different Server

**Before** (hardcoded):
```python
base_dir = "/home/chbope/extension/script/svanalysis"  # Line 742 in old version
```

**After** (configurable):
```python
SAMPLE_FILE = "/mnt/data/project/samples.txt"
VCF_DIR = "/mnt/data/project/vcf"
OUTPUT_DIR = "/mnt/data/project/results"
```

## Benefits

✅ **One place to configure**: All paths in configuration section at top
✅ **No code modification needed**: Just change configuration variables
✅ **Clear documentation**: Well-commented configuration section
✅ **Flexible quality thresholds**: Can adjust for different analyses
✅ **Easy troubleshooting**: Clear path information printed at runtime

## Files Modified

1. **sv_analysis_v2.py**:
   - Added CONFIGURATION SECTION (lines 28-59)
   - Updated main() to use configuration variables
   - Updated analyze_all_samples() to use MIN_SUPPORT, MIN_VAF
   - Updated parse_vcf() to use CLONAL_THRESHOLD

2. **CHANGES.md**:
   - Documented portability improvements

3. **CONFIGURATION_GUIDE.md** (NEW):
   - Comprehensive configuration guide

4. **PORTABILITY_UPDATE.md** (NEW):
   - This summary document

## Advanced: Adjusting Quality Parameters

For different cancer types or sequencing depths, you can adjust:

### Stricter Filtering (Higher Confidence)
```python
MIN_SUPPORT = 15       # More reads required
MIN_VAF = 0.15         # Higher VAF threshold
CLONAL_THRESHOLD = 0.4  # More stringent clonal
```

### More Permissive (Higher Sensitivity)
```python
MIN_SUPPORT = 5        # Fewer reads required
MIN_VAF = 0.05         # Lower VAF threshold
CLONAL_THRESHOLD = 0.2  # More permissive clonal
```

### GBM-Optimized (Current Default)
```python
MIN_SUPPORT = 10       # Balanced
MIN_VAF = 0.1          # Standard somatic
CLONAL_THRESHOLD = 0.3  # Standard clonal
```

## Verification

When you run the script, it will print the configuration being used:

```
================================================================================
STRUCTURAL VARIANT ANALYSIS FOR GBM SAMPLES (VERSION 2)
FOCUS: SIGNIFICANT VAF METRICS ONLY
================================================================================
Sample file: /your/path/to/sample.txt
VCF directory: /your/path/to/vcf/files
Output directory: /your/path/to/output/directory
================================================================================
```

Verify these paths are correct before the analysis proceeds.

## Backward Compatibility

The default configuration uses the same paths as before, so existing users don't need to make any changes unless they want to move the analysis to a different location.

## Next Steps

1. Review [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md) for detailed instructions
2. Configure paths for your system
3. Run the analysis: `python sv_analysis_v2.py`
4. Check results in your configured OUTPUT_DIR

## Summary

The script is now **fully portable** and can be easily moved between systems, users, and projects by simply modifying the configuration variables at the top of the file. No code changes are needed!
