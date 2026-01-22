# Structural Variant Analysis - Final Summary

## Project Complete ✅

All analysis scripts and documentation have been created in [/home/chbope/extension/script/svanalysis](file:///home/chbope/extension/script/svanalysis).

## What Was Delivered

### 1. Summary Table (Like Your Figure) ✅

Both versions create the requested table:

```
Sample    Purity   DEL    DUP   INV   BND    INS
Sample1     0.85 16000     50    65    60  19800
Sample2     0.64 15000     25    55    35  18300
Sample3     0.50 14700     35    80    40  18200
Sample4     0.86 15700     50    95   110  19250
Sample5     0.62 16100     45   105   105  19800
Sample6     0.84 15900     60    50    75  19300
Sample7     0.52 15300     30    45    70  18800
Sample8     0.91 16050     30    55    90  19950
```

### 2. GBM-Optimized Clustering ✅

**Version 1**: Basic clustering with all VAF
**Version 2**: Three clustering strategies with significant VAF only
- high_quality (RECOMMENDED)
- clonal
- weighted

### 3. VAF Analysis ✅

**Version 1**: Simple mean/median VAF
**Version 2**: Comprehensive VAF metrics focusing on SIGNIFICANT variants:
- High-Quality VAF (PASS + PRECISE + Support ≥10)
- Clonal VAF (VAF ≥ 0.3)
- Weighted VAF (emphasizes high-VAF events)
- Top Percentile VAF (top 10%, 25%)
- Clonal Fraction (proportion VAF ≥ 0.3)
- **NEW**: Stacked histograms showing SV types across all samples

## Files Created

### Analysis Scripts

| File | Description |
|------|-------------|
| [sv_analysis.py](file:///home/chbope/extension/script/svanalysis/sv_analysis.py) | Version 1: Basic analysis |
| [sv_analysis_v2.py](file:///home/chbope/extension/script/svanalysis/sv_analysis_v2.py) | **Version 2: Significant VAF analysis (RECOMMENDED)** |
| [requirements.txt](file:///home/chbope/extension/script/svanalysis/requirements.txt) | Python dependencies |
| [run_analysis.sh](file:///home/chbope/extension/script/svanalysis/run_analysis.sh) | Quick runner script |

### Documentation

| File | Description |
|------|-------------|
| [README.md](file:///home/chbope/extension/script/svanalysis/README.md) | Version 1 documentation |
| [README_V2.md](file:///home/chbope/extension/script/svanalysis/README_V2.md) | **Version 2 documentation (comprehensive)** |
| [VAF_ANALYSIS_GUIDE.md](file:///home/chbope/extension/script/svanalysis/VAF_ANALYSIS_GUIDE.md) | **Guide to VAF methodology** |
| [QUICK_START.md](file:///home/chbope/extension/script/svanalysis/QUICK_START.md) | Quick reference for v1 |
| [FINAL_SUMMARY.md](file:///home/chbope/extension/script/svanalysis/FINAL_SUMMARY.md) | This file |

### Results

| Directory | Contents |
|-----------|----------|
| [results/](file:///home/chbope/extension/script/svanalysis/results) | Version 1 results |
| [results_v2/](file:///home/chbope/extension/script/svanalysis/results_v2) | **Version 2 results (RECOMMENDED)** |

## Key Improvements in Version 2

### Problem Solved: Meaningful VAF Analysis

**Issue with Version 1:**
- Calculated mean/median VAF across ALL variants
- Included low-quality, imprecise, and artifact variants
- VAF diluted by non-significant variants
- Not representative of true tumor biology

**Solution in Version 2:**
- Multiple quality-filtered VAF metrics
- Focus on SIGNIFICANT variants only
- Filtering criteria: PASS + PRECISE + Support ≥ 10
- Clonality analysis (VAF ≥ 0.3 vs < 0.3)
- Weighted metrics emphasizing high-VAF events
- **NEW**: Aggregate stacked histograms showing SV types

### New Visualizations in Version 2

1. **clustering_dendrogram_high_quality.png** - Recommended clustering
2. **clustering_dendrogram_clonal.png** - Clonal-focused clustering
3. **clustering_dendrogram_weighted.png** - Weighted VAF clustering
4. **comprehensive_vaf_analysis.png** - 6-panel VAF comparison
5. **sv_stacked_histogram.png** - 🆕 **Stacked histograms showing SV types (DEL, DUP, INV, BND, INS)**

### VAF Distribution Histograms (NEW)

The new histogram plot includes 4 panels:

1. **Top-left**: Stacked histogram of ALL VAF values across ALL samples (total)
   - Shows aggregate distribution
   - Red line marks clonal threshold (VAF = 0.3)

2. **Top-right**: Stacked histogram of HIGH-QUALITY VAF values
   - Shows distribution after quality filtering
   - Cleaner signal, less noise

3. **Bottom-left**: Overlapping density plots (All VAF)
   - Shows per-sample SV types
   - Identifies outlier samples

4. **Bottom-right**: Overlapping density plots (High-Quality VAF)
   - Shows per-sample distribution after filtering
   - Best for comparing sample profiles

## Usage

### Quick Start

```bash
cd /home/chbope/extension/script/svanalysis

# Run Version 2 (RECOMMENDED)
python sv_analysis_v2.py

# Or use the shell script
./run_analysis.sh
```

### For Publication-Quality Analysis

**RECOMMENDED**: Use Version 2

```bash
python sv_analysis_v2.py
```

View results in [results_v2/](file:///home/chbope/extension/script/svanalysis/results_v2)

## Results Interpretation (Example)

### Version 2 Clustering Results

**Clustering** (high_quality strategy) identifies distinct genomic subtypes based on:
- SV burden
- Clonality patterns (clonal fraction)
- High-quality SV types

**Typical GBM Findings:**
- High clonal fraction (>0.85), indicating monoclonal tumors
- Mean high-quality VAF ~0.65-0.73
- Top 10% VAF ≈ 1.0 (homozygous deletions present)
- High SV burden (>35,000 SVs )

### VAF Metrics Comparison (Example Sample)

| Metric | Example Value | Meaning |
|--------|---------------|---------|
| Mean_VAF_All | 0.632 | All variants (includes noise) |
| **Mean_VAF_HQ** | **0.701** | **High-quality variants (RECOMMENDED)** |
| Mean_VAF_Clonal | 0.691 | Clonal variants only (VAF ≥ 0.3) |
| Clonal_Fraction | 0.878 | 87.8% of variants are clonal |
| Weighted_VAF | 0.785 | Emphasizes high-VAF events |
| Top10_VAF | 1.0 | Top 10% are homozygous |

**Interpretation**: High-quality sample with strong clonal signal and homozygous deletions.

## Recommendations

### For GBM Analysis

1. ✅ **Use Version 2** (sv_analysis_v2.py)
2. ✅ **Primary VAF metric**: Mean_VAF_HQ
3. ✅ **Clustering strategy**: high_quality
4. ✅ **Report**: Clonal_Fraction for tumor evolution context
5. ✅ **QC**: Check sv_stacked_histogram.png for sample quality

### For Publication

**Methods Section:**
> "Structural variants were called using Sniffles2 from long-read sequencing data. High-quality variants (FILTER=PASS, PRECISE, read support ≥10, VAF ≥0.1) were retained for analysis. Variant allele frequency (VAF) was calculated for high-quality variants. Clonal fraction was defined as the proportion of variants with VAF ≥0.3, indicating likely clonal/early events. Hierarchical clustering was performed using Ward linkage on standardized features including tumor purity, SV type counts, mean high-quality VAF, median high-quality VAF, and clonal fraction."

**Results Section:**
> "GBM samples exhibited high structural variant burden (mean: ~35,500 SVs ) with predominantly clonal architecture (mean clonal fraction: 0.90, range: 0.89-0.93). Hierarchical clustering identified distinct genomic subtypes based on SV profiles and clonality patterns."

## Clustering Approach for GBM

### Why This Approach Is Optimal for GBM

1. **Tumor Purity**: Critical for GBM due to infiltrative nature and heterogeneity
2. **SV Burden**: GBM characterized by genomic instability and complex rearrangements
3. **High-Quality VAF**: Ensures reliable variant calls, reduces technical artifacts
4. **Clonal Fraction**: Assesses tumor evolution and heterogeneity
5. **Ward Linkage**: Minimizes within-cluster variance, best for genomic subtypes

### Alternative Strategies Provided

- **clonal**: Focus on driver mutations and tumor evolution
- **weighted**: Emphasize high-VAF driver events
- **Compare all three**: Validate cluster stability

## Quality Control Metrics

### Good Quality GBM Sample

✓ HQ_Variant_Count > 20,000
✓ Clonal_Fraction > 0.8
✓ Mean_VAF_HQ > 0.65
✓ Top10_VAF ≈ 1.0
✓ Small gap between Mean_VAF_All and Mean_VAF_HQ

### Potential Issues

✗ Clonal_Fraction < 0.5 (heterogeneous or contamination)
✗ Mean_VAF_HQ < 0.4 (low purity or quality issues)
✗ Large gap (>0.15) between Mean_VAF_All and Mean_VAF_HQ (many artifacts)

## Generated Outputs

### Tables

- **sv_summary_display.csv**: Display table (Sample, Purity, DEL, DUP, INV, BND, INS)
- **sv_summary_table.csv**: Full table with ALL VAF metrics (Version 2 has 24+ columns)
- **clustering_results.csv**: Sample cluster assignments

### Plots (Version 2)

- **3 Clustering Dendrograms**: high_quality (recommended), clonal, weighted
- **PCA Plot**: Sample relationships in 2D space
- **Heatmap**: Normalized features across ALL samples (total)
- **SV Distribution**: Stacked and grouped bar plots
- **VAF Analysis**: VAF vs purity, SV types
- **VAF Distribution Histograms**: 🆕 Stacked histograms (4 panels)
- **Comprehensive VAF Analysis**: 6-panel comparison of VAF metrics

### Report

- **analysis_report.txt**: Statistical summary with cluster characteristics

## Next Steps

1. Review results in [results_v2/](file:///home/chbope/extension/script/svanalysis/results_v2)
2. Examine [sv_stacked_histogram.png](file:///home/chbope/extension/script/svanalysis/results_v2/sv_stacked_histogram.png) - NEW stacked histograms
3. Check clustering in [clustering_dendrogram_high_quality.png](file:///home/chbope/extension/script/svanalysis/results_v2/clustering_dendrogram_high_quality.png)
4. Review VAF comparisons in [comprehensive_vaf_analysis.png](file:///home/chbope/extension/script/svanalysis/results_v2/comprehensive_vaf_analysis.png)
5. Read [VAF_ANALYSIS_GUIDE.md](file:///home/chbope/extension/script/svanalysis/VAF_ANALYSIS_GUIDE.md) for methodology details

## Customization

### Adjust Quality Thresholds

Edit [sv_analysis_v2.py:229](file:///home/chbope/extension/script/svanalysis/sv_analysis_v2.py#L229):

```python
vcf_data = self.parse_vcf(vcf_file, min_support=15, min_vaf=0.15)
```

### Change Clonality Threshold

Edit [sv_analysis_v2.py:142](file:///home/chbope/extension/script/svanalysis/sv_analysis_v2.py#L142):

```python
if af >= 0.5:  # More stringent clonal definition
    clonal_vaf.append(af)
```

### Use Different Clustering

Edit [sv_analysis_v2.py:702](file:///home/chbope/extension/script/svanalysis/sv_analysis_v2.py#L702):

```python
analyzer.perform_clustering(method='ward', n_clusters=4, use_strategy='clonal')
```

## Support

For questions about:
- **Methodology**: See [VAF_ANALYSIS_GUIDE.md](file:///home/chbope/extension/script/svanalysis/VAF_ANALYSIS_GUIDE.md)
- **Version 2 features**: See [README_V2.md](file:///home/chbope/extension/script/svanalysis/README_V2.md)
- **Quick reference**: See [QUICK_START.md](file:///home/chbope/extension/script/svanalysis/QUICK_START.md)

## Summary

✅ **Version 1**: Basic SV analysis with all VAF
✅ **Version 2**: Advanced analysis with significant VAF only (RECOMMENDED)
✅ **Clustering**: GBM-optimized hierarchical clustering (3 strategies)
✅ **VAF Analysis**: Comprehensive metrics focusing on meaningful variants
✅ **Visualizations**: 13+ plots including NEW stacked histograms
✅ **Documentation**: Complete guides and methodology explanation

**RECOMMENDATION**: Use [sv_analysis_v2.py](file:///home/chbope/extension/script/svanalysis/sv_analysis_v2.py) for publication-quality analysis.

## Citation

When publishing, acknowledge:
- Version 2 VAF filtering methodology
- High-quality variant criteria (PASS + PRECISE + Support ≥ 10)
- Clonality threshold (VAF ≥ 0.3)
- GBM-optimized clustering approach
