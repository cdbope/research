# Structural Variant Analysis for GBM Samples

This pipeline analyzes structural variants (SVs) from Sniffles VCF files for Glioblastoma (GBM) samples, incorporating tumor purity and variant allele frequency (VAF) metrics.

## Overview

The analysis performs:
1. **SV Count Table** - Generates a summary table showing DEL, DUP, INV, BND, and INS counts per sample
2. **Clustering Analysis** - Hierarchical clustering optimized for GBM samples
3. **VAF Analysis** - Calculates and analyzes variant allele frequencies
4. **Visualizations** - Multiple plots including dendrograms, PCA, heatmaps, and distributions

## Input Files

- `sample.txt` - Tab-separated file with sample IDs and tumor purity values:
  ```
  Sample1	0.84
  Sample2	0.64
  ...
  ```

- VCF files - Gzipped Sniffles VCF files named as `{sample_id}.wf_sv.vcf.gz`

## Installation

### Dependencies

Install required Python packages:

```bash
pip install -r requirements.txt
```

Required packages:
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn

## Usage

### Run Complete Analysis

```bash
python sv_analysis.py
```

This will:
1. Parse all VCF files listed in `sample.txt`
2. Extract SV counts (DEL, DUP, INV, BND, INS)
3. Calculate VAF metrics for each sample
4. Perform hierarchical clustering
5. Generate all visualizations
6. Create comprehensive report

### Output Files

All results are saved to `results/` directory:

#### Tables
- `sv_summary_table.csv` - Full summary with all metrics including VAF
- `sv_summary_display.csv` - Display table (Sample, Purity, DEL, DUP, INV, BND, INS)
- `clustering_results.csv` - Sample assignments to clusters

#### Plots
- `clustering_dendrogram.png` - Hierarchical clustering dendrogram
- `pca_clustering.png` - PCA plot showing sample relationships
- `sv_heatmap.png` - Heatmap of normalized SV features
- `sv_distribution.png` - SV type distribution (stacked and grouped bars)
- `vaf_analysis.png` - VAF vs purity and VAF distribution plots

#### Report
- `analysis_report.txt` - Comprehensive text report with statistics

## Analysis Details

### Clustering Methodology

The clustering analysis is optimized for GBM samples and considers:

1. **Features used**:
   - Tumor purity
   - SV counts (DEL, DUP, INV, BND, INS)
   - Mean VAF (clonality indicator)
   - Median VAF

2. **Method**:
   - Hierarchical clustering with Ward linkage
   - Euclidean distance metric
   - Features are standardized (z-score normalization)

3. **Why this approach for GBM**:
   - **Purity**: Critical for GBM as tumor heterogeneity affects SV detection
   - **SV burden**: High genomic instability is characteristic of GBM
   - **VAF**: Indicates clonal vs subclonal events, important for understanding GBM evolution
   - **Ward linkage**: Minimizes within-cluster variance, good for identifying distinct genomic subtypes

### VAF Interpretation

- **Mean VAF**: Average allele frequency across all SVs
- **Median VAF**: Robust measure less affected by outliers
- **VAF vs Purity correlation**: Higher purity samples typically show higher VAF for clonal events

For GBM samples:
- High VAF (>0.5) suggests clonal/early events
- Low VAF (<0.3) suggests subclonal/late events
- VAF distribution indicates tumor heterogeneity

### Suggested Clustering Strategies

The default analysis uses 3 clusters, but you can modify based on:

1. **Genomic instability groups**:
   - Low SV burden (stable genomes)
   - Medium SV burden
   - High SV burden (chromothripsis/chromoplexy)

2. **Clonality groups**:
   - High clonal VAF (early driver events)
   - Mixed clonality (heterogeneous)
   - Low clonal VAF (subclonal diversity)

3. **Combined approach** (recommended):
   - Integrate SV burden, purity, and VAF
   - Identifies clinically relevant GBM subtypes

## Customization

### Modify Clustering Parameters

Edit `sv_analysis.py` in the `main()` function:

```python
# Change clustering method
analyzer.perform_clustering(method='average', n_clusters=4)

# Options:
# method: 'ward', 'average', 'complete', 'single'
# n_clusters: integer (default is 3)
```

### Add Custom Features

Extend the `SVAnalyzer` class to include additional metrics from VCF INFO fields.

## Expected Output Table Format

```
Sample      Purity  DEL    DUP   INV   BND   INS
Sample1    0.84    205    38    14    92    55
Sample2    0.64    134    22    10    55    28
...
```

## Notes

- VCF files must contain `SVTYPE` in INFO field
- `AF` (allele frequency) field is used for VAF calculation
- Samples without VCF files are skipped with a warning
- Missing AF values are handled gracefully

## GBM-Specific Considerations

Glioblastoma multiforme (GBM) is characterized by:
- High genomic instability
- Frequent copy number alterations
- Intratumoral heterogeneity
- Complex structural rearrangements

This pipeline accounts for these characteristics by:
1. Including tumor purity in clustering
2. Analyzing VAF to assess clonality
3. Using robust normalization for varying SV burdens
4. Providing multiple visualization methods

## Troubleshooting

**Issue**: Module not found errors
- **Solution**: Install dependencies with `pip install -r requirements.txt`

**Issue**: VCF file not found warning
- **Solution**: Ensure VCF files are named exactly as `{sample_id}.wf_sv.vcf.gz`

**Issue**: Empty VAF values
- **Solution**: Check that VCF files contain `AF=` in INFO field

## Citation

If you use this pipeline, please cite:
- Sniffles2 for SV calling
- Your GBM study details

## Contact

For questions or issues, please contact the analysis team.
