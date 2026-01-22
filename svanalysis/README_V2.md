# Structural Variant Analysis for GBM Samples - VERSION 2

## NEW in Version 2: Focus on SIGNIFICANT VAF Only

Version 2 addresses the critical issue of **meaningful VAF analysis** by filtering out low-quality and non-significant variants that would otherwise skew the average VAF.

### Key Problem Solved

**Version 1 Issue**: Taking mean/median VAF across ALL variants includes:
- Low-quality variants (imprecise, low support)
- Passenger mutations with low VAF
- Chromosome-specific artifacts
- Non-significant variants that dilute the signal

**Version 2 Solution**: Multiple sophisticated VAF metrics focusing on **significant variants only**:

1. **High-Quality VAF** (RECOMMENDED)
   - Filters: PASS + PRECISE + Support ≥ 10
   - Uses only reliable, well-supported variants
   - Best for clustering and analysis

2. **Clonal VAF** (VAF ≥ 0.3)
   - Focuses on likely driver events
   - Represents early/clonal mutations
   - Indicates tumor evolution

3. **Weighted VAF**
   - Gives more weight to high-VAF variants
   - Emphasizes significant events
   - De-emphasizes noise

4. **Top Percentile VAF** (Top 10%, 25%)
   - Analyzes only the most significant variants
   - Focuses on high-confidence calls
   - Reduces artifact influence

## What's New in Version 2

### Enhanced VAF Metrics

#### Quality Filtering
```python
# High-Quality Variant Criteria:
- FILTER = "PASS"
- PRECISE (not IMPRECISE)
- SUPPORT ≥ 10 reads
- VAF ≥ 0.1 (minimum threshold)
```

#### Clonality Analysis
```python
# Clonal Variants: VAF ≥ 0.3 (likely driver events)
# Subclonal Variants: 0.1 ≤ VAF < 0.3 (heterogeneous/late events)
# Clonal Fraction = proportion of variants with VAF ≥ 0.3
```

### Multiple Clustering Strategies

Version 2 performs clustering with **3 different VAF strategies** to compare results:

1. **high_quality** (RECOMMENDED)
   - Uses: Mean_VAF_HQ, Median_VAF_HQ, Clonal_Fraction
   - Best for: Reliable genomic subtype identification

2. **clonal**
   - Uses: Mean_VAF_Clonal, Median_VAF_Clonal, Clonal_Fraction
   - Best for: Tumor evolution analysis

3. **weighted**
   - Uses: Weighted_VAF, Clonal_Fraction
   - Best for: Emphasizing high-VAF driver events

### Comprehensive VAF Analysis Plot

New 6-panel visualization comparing:
1. All VAF vs High-Quality VAF
2. Clonal Fraction vs Purity
3. Comparison of VAF Metrics
4. Clonal vs Subclonal Counts
5. High-Quality SV Distribution
6. Top Percentile VAF Comparison

## Installation

Same as Version 1:

```bash
pip install -r requirements.txt
```

## Usage

### Run Version 2 Analysis

```bash
python sv_analysis_v2.py
```

### Output Location

Results are saved to `results_v2/` directory (separate from v1).

## Output Files

### Tables

1. **sv_summary_table.csv** - Complete table with ALL VAF metrics:
   - Standard SV counts (DEL, DUP, INV, BND, INS)
   - High-quality SV counts (DEL_HQ, DUP_HQ, etc.)
   - Mean_VAF_All, Median_VAF_All (all variants, for reference)
   - **Mean_VAF_HQ, Median_VAF_HQ** (RECOMMENDED metrics)
   - Mean_VAF_Clonal, Median_VAF_Clonal
   - Clonal_Fraction (proportion VAF ≥ 0.3)
   - Clonal_Count, Subclonal_Count
   - Weighted_VAF
   - Top10_VAF, Top25_VAF
   - HQ_Variant_Count

2. **sv_summary_display.csv** - Display table (same as v1)

3. **clustering_results.csv** - Cluster assignments (using high_quality strategy)

### Clustering Dendrograms (3 strategies)

- **clustering_dendrogram_high_quality.png** ✅ RECOMMENDED
- clustering_dendrogram_clonal.png
- clustering_dendrogram_weighted.png

### Visualizations

- **comprehensive_vaf_analysis.png** 🆕 NEW - 6-panel VAF comparison
- pca_clustering.png (using high-quality VAF)
- sv_heatmap.png (using high-quality VAF)
- sv_distribution.png
- vaf_analysis.png (using high-quality VAF)

### Report

- analysis_report.txt (with enhanced VAF statistics)

## Understanding the Results

### VAF Metric Comparison

| Metric | Description | Best For |
|--------|-------------|----------|
| **Mean_VAF_HQ** | Mean of high-quality variants | Clustering, reliable analysis |
| Mean_VAF_All | Mean of all variants | Reference only |
| Mean_VAF_Clonal | Mean of VAF ≥ 0.3 | Identifying driver events |
| Weighted_VAF | VAF-weighted mean | Emphasizing significant variants |
| Top10_VAF | Mean of top 10% VAF | Most confident calls |
| Clonal_Fraction | Proportion VAF ≥ 0.3 | Tumor clonality |

### Example Results

```
Sample: Sample1 (Purity = 0.84)
- Mean_VAF_All: 0.632 (all 36,070 variants)
- Mean_VAF_HQ: 0.701 (25,908 high-quality variants)
- Mean_VAF_Clonal: 0.691 (31,569 clonal variants)
- Clonal_Fraction: 0.878 (87.8% clonal)
- Weighted_VAF: 0.785
- Top10_VAF: 1.0 (top variants are homozygous)
```

**Interpretation**:
- High purity sample with mostly clonal variants
- Strong signal (high clonal fraction)
- High-quality VAF (0.701) more accurate than all VAF (0.632)
- Top variants show homozygous deletions/events

### Clustering Results Comparison

Version 2 provides different clustering with each strategy:

**High-Quality Strategy** (RECOMMENDED):
- Cluster 1: Sample1, Sample4, Sample5, Sample6, Sample8
- Cluster 2: Sample2, Sample3
- Cluster 3: Sample7

**Clonal Strategy**:
- Cluster 1: Sample2, Sample3, Sample7
- Cluster 2: Sample1, Sample6, Sample8
- Cluster 3: Sample4, Sample5

**Weighted Strategy**:
- Same as High-Quality

### Which Strategy to Use?

**For GBM Analysis, we RECOMMEND: high_quality strategy**

Reasons:
1. ✅ Filters out artifacts and low-quality calls
2. ✅ Balances sensitivity and specificity
3. ✅ Includes clonal fraction for evolution analysis
4. ✅ Most robust to technical variation
5. ✅ Best for identifying genomic subtypes

**Use clonal strategy when**:
- Focusing specifically on driver mutations
- Studying tumor evolution
- Comparing clonal architecture

**Use weighted strategy when**:
- High-VAF events are most important
- Want to de-emphasize subclonal events
- Homogeneous tumors

## Version 1 vs Version 2 Comparison

| Feature | Version 1 | Version 2 |
|---------|-----------|-----------|
| VAF Calculation | Mean/Median of ALL variants | Multiple quality-filtered metrics |
| Quality Filter | None | PASS + PRECISE + Support ≥ 10 |
| Clonality Analysis | No | Yes (VAF ≥ 0.3) |
| Clustering Strategies | 1 (basic) | 3 (high_quality, clonal, weighted) |
| VAF Plots | 2 panels | 8 panels (2 basic + 6 comprehensive) |
| High-Quality SV Counts | No | Yes (DEL_HQ, DUP_HQ, etc.) |
| Weighted VAF | No | Yes |
| Top Percentile VAF | No | Yes (Top 10%, 25%) |
| **Recommendation** | Good for initial exploration | **BETTER for publication-quality analysis** |

## Recommendations for GBM Analysis

### For Clustering
1. **Primary analysis**: Use `high_quality` strategy
2. **Validation**: Compare with `clonal` and `weighted` strategies
3. **Key metrics**: Mean_VAF_HQ, Clonal_Fraction

### For Interpretation
1. Focus on High-Quality VAF metrics (ignore "All VAF")
2. Use Clonal_Fraction to assess tumor homogeneity
3. Compare Top10_VAF across samples to identify homozygous events
4. Check comprehensive_vaf_analysis.png for quality assessment

### Quality Control
- High Clonal_Fraction (>0.8) = homogeneous tumor, good quality
- Low Clonal_Fraction (<0.5) = heterogeneous tumor or quality issues
- Large difference between Mean_VAF_All and Mean_VAF_HQ = many low-quality variants
- Top10_VAF near 1.0 = presence of homozygous deletions/amplifications

## Advanced Customization

### Adjust Quality Thresholds

Edit `sv_analysis_v2.py`:

```python
# Line 229: Change min_support and min_vaf
vcf_data = self.parse_vcf(vcf_file, min_support=15, min_vaf=0.15)
```

### Adjust Clonality Threshold

Edit `sv_analysis_v2.py`:

```python
# Line 142: Change clonal VAF threshold from 0.3 to custom value
if af >= 0.5:  # More stringent clonal definition
    clonal_vaf.append(af)
```

### Use Different Clustering Strategy

Edit `sv_analysis_v2.py`:

```python
# Line 702: Change use_strategy parameter
analyzer.perform_clustering(method='ward', n_clusters=3, use_strategy='clonal')
```

## Citation

If using this analysis for publication, please acknowledge:
- Version 2 VAF filtering methodology
- High-quality variant criteria (PASS + PRECISE + Support ≥ 10)
- Clonality analysis (VAF ≥ 0.3 threshold)

## Troubleshooting

**Q: Why are some VAF metrics zero?**
A: Sample may not have enough high-quality or clonal variants passing filters. Check HQ_Variant_Count and Clonal_Count.

**Q: Which VAF metric should I report?**
A: Use **Mean_VAF_HQ** for primary analysis. Report **Clonal_Fraction** for tumor evolution context.

**Q: Different clustering results between strategies?**
A: Expected. High_quality is most robust. Use it for primary analysis.

**Q: How to interpret Clonal_Fraction?**
A: Proportion of variants with VAF ≥ 0.3.
- High (>0.8): Clonal/homogeneous tumor
- Medium (0.5-0.8): Mixed clonality
- Low (<0.5): Highly heterogeneous or subclonal

## Contact

For questions about Version 2 methodology or results interpretation, consult the analysis team.
