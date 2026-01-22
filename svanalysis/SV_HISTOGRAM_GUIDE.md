# SV Aggregate Stacked Histogram Guide

## Overview

The **sv_stacked_histogram.png** shows the TOTAL aggregate counts of each structural variant type across ALL samples combined.

## Visualization Description

The plot contains **TWO side-by-side aggregate bars**:

### Left Panel: All Structural Variants
- Single stacked bar showing total SV counts aggregated across all 8 samples
- Each color segment represents a different SV type:
  - **DEL** (Deletions) - typically the largest segment
  - **DUP** (Duplications) - usually small
  - **INV** (Inversions) - usually small
  - **BND** (Breakends/Translocations) - usually small
  - **INS** (Insertions) - typically large, second to DEL

- **Count labels** are shown within each segment
- **Y-axis**: Total SV count (all samples combined)

### Right Panel: High-Quality Structural Variants
- Same format as left panel
- Only includes variants passing quality filters:
  - FILTER = PASS
  - PRECISE (not IMPRECISE)
  - Read support ≥ 10

- Shows effect of quality filtering
- Compare to left panel to see how many variants are filtered out

## Example Results (Your Data)

### All SVs (Total across 8 samples):
```
Total: 279,776 SVs
- DEL: 124,932 (44.7%)
- INS: 153,379 (54.8%)
- DUP: 323 (0.1%)
- INV: 550 (0.2%)
- BND: 592 (0.2%)
```

### High-Quality SVs (Total across 8 samples):
```
Total: 214,670 SVs
- DEL_HQ: 98,525 (45.9%)
- INS_HQ: ~115,000 (53.6%)
- DUP_HQ: 160 (0.1%)
- INV_HQ: 347 (0.2%)
- BND_HQ: ~600 (0.3%)
```

## Interpretation

### SV Type Distribution
- **Deletions and Insertions dominate** (~99% of all SVs)
  - This is typical for GBM and long-read sequencing
- **Duplications, Inversions, Breakends are rare** (~1%)
  - These complex rearrangements are less frequent but potentially important

### Quality Filtering Effect
Compare left vs right panels:
- **~76.7% of variants pass quality filters** (214,670 / 279,776)
- Filtering removes ~65,000 low-quality variants
- DEL and INS maintain similar proportions after filtering
- Proportion changes indicate which SV types have more low-quality calls

### GBM-Specific Insights
1. **High DEL count**: Characteristic of GBM (e.g., CDKN2A deletions)
2. **High INS count**: Long-read sequencing advantage (better insertion detection)
3. **Low complex SV count**: DUP, INV, BND are rarer but may be drivers
4. **Good quality**: 76.7% passing strict filters indicates good data quality

## Differences from Per-Sample Plots

This aggregate histogram is **different** from the per-sample plots elsewhere:

| Feature | Aggregate Histogram | Per-Sample Plots |
|---------|-------------------|------------------|
| **Granularity** | Total across ALL samples | Individual sample breakdown |
| **Purpose** | Overall cohort composition | Sample-to-sample variation |
| **X-axis** | Single bar (aggregate) | Multiple bars (one per sample) |
| **When to use** | Cohort summary, publication | Sample QC, outlier detection |

## When to Use This Plot

### Use Aggregate Histogram For:
✅ **Cohort-level summary**: "Our GBM cohort had 280k SVs total"
✅ **SV type distribution**: "Deletions and insertions comprised 99% of SVs"
✅ **Quality filtering report**: "77% of variants passed quality filters"
✅ **Publication figure**: Clean, simple visualization of overall composition

### Use Per-Sample Plots For:
✅ **Sample QC**: Identify outliers with unusual SV counts
✅ **Batch effects**: Compare SV distributions across batches
✅ **Sample-specific analysis**: Focus on individual sample characteristics
✅ **Detailed exploration**: See which samples drive the aggregate patterns

## Comparison to Other Plots

The analysis generates multiple SV distribution plots:

1. **sv_stacked_histogram.png** (THIS PLOT)
   - AGGREGATE total counts
   - 2 bars (all vs high-quality)
   - Cohort-level summary

2. **sv_distribution.png**
   - PER-SAMPLE stacked and grouped bar plots
   - Shows individual sample contributions
   - Sample-level detail

3. **sv_heatmap.png**
   - Per-sample normalized features
   - Includes purity, VAF, and clonal fraction
   - Multi-dimensional comparison

## Example Use in Publication

**Methods:**
> "Across the 8 GBM samples, we identified 279,776 structural variants, of which 214,670 (77%) passed quality filters (FILTER=PASS, PRECISE, read support ≥10). Deletions (124,932, 45%) and insertions (153,379, 55%) comprised the majority of variants, with duplications, inversions, and breakends each representing <1% of total SVs (Figure X)."

**Figure Legend:**
> "Figure X: Aggregate distribution of structural variant types across GBM cohort. Left panel shows all detected SVs (n=279,776); right panel shows high-quality filtered SVs (n=214,670). DEL: deletions, DUP: duplications, INV: inversions, BND: breakends, INS: insertions."

## Customization

To modify the appearance, edit `sv_analysis_v2.py` lines 439-509:

### Change Colors
```python
colors = plt.cm.Spectral(np.linspace(0, 1, len(sv_types)))  # Different colormap
```

### Adjust Bar Width
```python
axes[0].bar(0, count, bottom=bottom, width=0.8)  # Wider bar (default 0.6)
```

### Modify Labels
```python
if count > 500:  # Label smaller segments
    axes[0].text(0, bottom + count/2, f'{count:,}', ...)
```

## Technical Details

**Data Source**: Summed from `sv_summary_table.csv`
- All SVs: Sum of DEL, DUP, INV, BND, INS columns
- HQ SVs: Sum of DEL_HQ, DUP_HQ, INV_HQ, BND_HQ, INS_HQ columns

**File Location**: `results_v2/sv_stacked_histogram.png`

**Resolution**: 300 DPI (publication quality)

**Format**: PNG with tight bounding box
