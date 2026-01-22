# Quick Start Guide

## Run Analysis

```bash
cd /home/chbope/extension/script/svanalysis
python sv_analysis.py
```

## View Results

All outputs are in the `results/` directory:

### Key Tables
- **sv_summary_display.csv** - Main table (like your figure)
- **sv_summary_table.csv** - Full table with VAF metrics
- **clustering_results.csv** - Cluster assignments

### Visualizations
- **clustering_dendrogram.png** - Hierarchical clustering tree
- **pca_clustering.png** - PCA plot with clusters
- **sv_heatmap.png** - Heatmap of all features
- **sv_distribution.png** - SV counts by type
- **vaf_analysis.png** - VAF vs purity analysis

### Report
- **analysis_report.txt** - Complete statistical summary

## Current Results Summary

### Sample Table (like your figure)
```
Sample      Purity   DEL    DUP   INV   BND    INS
Sample1      0.84 16083     53    66    63  19805
Sample2      0.64 15064     24    54    36  18305
Sample3      0.50 14669     33    80    40  18160
Sample4       0.86 15727     50    96   110  19250
Sample5       0.62 16129     43   103   103  19801
Sample6       0.84 15878     60    50    76  19302
Sample7       0.52 15330     28    47    72  18817
Sample8       0.91 16052     32    54    92  19939
```

### Clustering Results (3 clusters)

**Cluster 1** (Lower purity, highest VAF)
- Samples: Sample2, Sample3, Sample7
- Mean Purity: 0.55
- Mean VAF: 0.661

**Cluster 2** (High purity, mid VAF)
- Samples: Sample1, Sample6, Sample8
- Mean Purity: 0.86
- Mean VAF: 0.641

**Cluster 3** (Mid purity, mid VAF)
- Samples: Sample4, Sample5
- Mean Purity: 0.74
- Mean VAF: 0.645

### Key Findings

1. **SV Burden**: All samples show high structural variant burden typical of GBM
   - Deletions: ~15,000-16,000 per sample
   - Insertions: ~18,000-20,000 per sample
   - Low duplications, inversions, and breakends

2. **VAF Analysis**: Mean VAF ranges from 0.630-0.674
   - Relatively high VAF suggests clonal events
   - Consistent with aggressive GBM phenotype

3. **Clustering Insights**:
   - Cluster 1: Lower purity samples (potential higher normal contamination)
   - Cluster 2: Highest purity samples (best tumor representation)
   - Cluster 3: Intermediate characteristics

## Clustering Approach for GBM

The analysis uses an optimal approach for GBM samples:

1. **Hierarchical clustering with Ward linkage**
   - Best for identifying distinct genomic subtypes
   - Minimizes within-cluster variance

2. **Features considered**:
   - Tumor purity (critical for GBM heterogeneity)
   - All SV types (DEL, DUP, INV, BND, INS)
   - Mean and median VAF (clonality indicators)

3. **Why this works for GBM**:
   - Accounts for tumor purity variations
   - Captures genomic instability patterns
   - Incorporates clonality information via VAF
   - Robust to GBM's high mutational burden

## Next Steps

1. Review the visualizations in `results/`
2. Check clustering assignments in `clustering_results.csv`
3. Examine VAF patterns in `vaf_analysis.png`
4. Review detailed statistics in `analysis_report.txt`

## Customize Analysis

To change clustering parameters, edit `sv_analysis.py`:

```python
# Change number of clusters
analyzer.perform_clustering(method='ward', n_clusters=4)

# Try different linkage methods
analyzer.perform_clustering(method='average')  # or 'complete', 'single'
```
