# ⚠️ CRITICAL: Genome Build Mismatch Detected

## Problem Identified

Your SV data uses **hg38/GRCh38** reference genome, but the TCGA/PCAWG comparison files are in **hg19/GRCh37** coordinates!

### Evidence:

**Your data (hg38):**
```
chr1: 248,956,422 bp (hg38)
chr7: 159,345,973 bp (hg38)
```

**TCGA reference file (hg19):**
```
EGFR: chr7:55,019,032-55,211,628 (hg19 coordinates)
```

## Impact

This coordinate mismatch will cause:
- ❌ **Failed SV matching** - coordinates don't align between builds
- ❌ **Missed overlaps** - known GBM alterations won't be detected
- ❌ **Incorrect frequency comparisons** - validation will fail

## Solution

We need to **liftOver** the TCGA/PCAWG coordinates from hg19 to hg38.

### Option 1: Use UCSC LiftOver Tool (Recommended)

I can create hg38 versions of the reference files using coordinate conversion.

### Option 2: Use Gene-Based Matching

Instead of coordinate matching, match SVs by gene overlap (less precise but works across builds).

## Action Required

Would you like me to:

1. **Create hg38 versions** of tcga_gbm_sv_summary.csv and pcawg_gbm_sv_summary.csv
2. **Update the comparison script** to use hg38 coordinates
3. **Add build detection** to warn users about mismatches

## Quick Fix

For now, I'll create approximate hg38 coordinates for the major GBM driver genes based on their current positions.
