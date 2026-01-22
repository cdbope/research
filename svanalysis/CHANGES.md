# Recent Changes

## Latest Update - Portability Improvements

### Changes Made:

1. **Added Configurable Paths for Portability**
   - Added CONFIGURATION SECTION at top of [sv_analysis_v2.py](sv_analysis_v2.py) (lines 28-59)
   - Three main configurable paths:
     - `SAMPLE_FILE`: Path to sample.txt file
     - `VCF_DIR`: Directory containing VCF files
     - `OUTPUT_DIR`: Directory for results output
   - Optional quality filtering parameters:
     - `MIN_SUPPORT`: Minimum read support (default: 10)
     - `MIN_VAF`: Minimum VAF threshold (default: 0.1)
     - `CLONAL_THRESHOLD`: Clonality threshold (default: 0.3)

2. **Updated Script to Use Configuration Variables**
   - Modified `main()` function to use configuration paths
   - Updated `analyze_all_samples()` to use MIN_SUPPORT and MIN_VAF
   - Updated `parse_vcf()` to use CLONAL_THRESHOLD
   - All hardcoded paths removed

3. **Created Configuration Guide**
   - New file: [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md)
   - Comprehensive instructions for path configuration
   - Example configurations for different systems
   - Troubleshooting guide
   - Advanced quality threshold adjustment guide

### Benefits:

✅ **Easy portability**: Change paths in one place
✅ **Clear documentation**: Configuration section is well-commented
✅ **No code changes needed**: Just modify configuration variables
✅ **Flexible**: Can adjust quality thresholds for different analyses

### How to Use:

1. Open [sv_analysis_v2.py](sv_analysis_v2.py)
2. Locate CONFIGURATION SECTION (lines 28-59)
3. Modify paths to match your system
4. Run: `python sv_analysis_v2.py`

See [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md) for detailed instructions.

---

## Previous Update - SV Stacked Histogram

### Changes Made:

1. **Replaced VAF Distribution Histograms with SV Type Stacked Histogram**
   - Removed: `vaf_distribution_histograms.png` (4-panel VAF distribution)
   - Added: `sv_stacked_histogram.png` (2-panel SV type histogram)

2. **New SV Stacked Histogram Features:**
   - **Left panel**: All structural variants stacked by type (DEL, DUP, INV, BND, INS)
   - **Right panel**: High-quality structural variants stacked by type
   - Shows TOTAL aggregate counts across ALL samples combined
   - Easy visualization of SV type composition 

3. **Anonymized All Sample IDs**
   - Replaced all real sample IDs with generic names (Sample1-Sample8)
   - Updated in all markdown documentation files:
     - README.md
     - README_V2.md
     - QUICK_START.md
     - VAF_ANALYSIS_GUIDE.md
     - FINAL_SUMMARY.md
   - Sample IDs in actual analysis output files (CSV) remain unchanged

### Files Updated:

- [sv_analysis_v2.py](file:///home/chbope/extension/script/svanalysis/sv_analysis_v2.py) - Replaced `plot_vaf_distribution_histogram()` with `plot_sv_stacked_histogram()`
- All *.md files - Sample ID anonymization

### New Output:

**sv_stacked_histogram.png** contains:

**Panel 1 (Left) - All SVs:**
- Stacked bar chart showing total SV counts
- Each color represents a different SV type (DEL, DUP, INV, BND, INS)
- X-axis: Samples
- Y-axis: Total SV count
- Shows overall SV burden and composition

**Panel 2 (Right) - High-Quality SVs:**
- Same format as Panel 1
- Only includes high-quality variants (PASS + PRECISE + Support ≥10)
- Shows filtered SV burden and composition
- Compares to left panel to see filtering effect

### Usage:

Run the analysis as before:

```bash
cd /home/chbope/extension/script/svanalysis
python sv_analysis_v2.py
```

The SV stacked histogram will be generated in:
`results_v2/sv_stacked_histogram.png`

### Why This Change:

1. **More informative for SV analysis**: Direct visualization of SV type distribution
2. **Cleaner presentation**: 2 panels instead of 4
3. **Easier interpretation**: Stacked bars show both total burden and type composition
4. **Focused on main analysis goal**: SV characterization, not VAF distribution

### Note on Confidentiality:

All markdown documentation files now use generic sample names (Sample1, Sample2, etc.) to protect confidential information. The actual analysis results files (CSV, PNG) still contain original sample IDs for your use, but documentation examples are anonymized.
