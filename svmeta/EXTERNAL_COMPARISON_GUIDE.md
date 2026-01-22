# External Dataset Comparison Guide

## Overview

The external dataset comparison module allows you to compare your GBM cohort's structural variant patterns against published reference datasets (TCGA-GBM and PCAWG). This helps:

1. **Validate findings** - Confirm your cohort shows known GBM SV patterns
2. **Identify novel SVs** - Find cohort-specific alterations not in reference datasets
3. **Contextualize results** - Understand if your findings are universal or cohort-specific
4. **Benchmark quality** - Verify your SV calling is detecting expected events

## Quick Start

### Step 1: Obtain External Data

You have three options:

#### Option A: Use Literature-Based Reference (Easiest)

The pipeline includes a curated reference file with well-characterized GBM SVs from TCGA publications:

```bash
# This file is created automatically by download_external_datasets.sh
./download_external_datasets.sh

# Use the literature reference file
ls external_datasets/known_gbm_svs_literature.csv
```

This file contains:
- EGFR amplifications (~45% frequency in TCGA)
- CDKN2A/B deletions (~52%)
- PTEN deletions (~41%)
- PDGFRA, CDK4, MDM2 amplifications
- RB1, NF1 deletions

**Update the script to use this file**:

Edit `04_external_dataset_comparison.py`:
```python
TCGA_SV_FILE = f"{EXTERNAL_DATA_DIR}/known_gbm_svs_literature.csv"
```

#### Option B: Download TCGA-GBM Data

1. Visit GDC Data Portal: https://portal.gdc.cancer.gov/
2. Filter: Project = TCGA-GBM, Data Category = Copy Number Variation or Structural Variation
3. Download aggregated data or individual VCFs
4. Process into summary format (see below)

#### Option C: Access PCAWG Data

1. Register at ICGC Data Portal: https://dcc.icgc.org/
2. Navigate to PCAWG datasets
3. Download glioblastoma SV calls
4. Process into summary format

### Step 2: Prepare Data Format

External datasets should be CSV files with these columns:

| Column | Type | Description |
|--------|------|-------------|
| sv_id | string | Unique identifier (e.g., "chr7:55M-56M_DUP") |
| chr | string | Chromosome |
| start | int | Start position |
| end | int | End position |
| svtype | string | SV type (DEL, DUP, INV, BND, INS) |
| frequency | float | Frequency in dataset (0-1) |
| num_samples | int | Number of samples with this SV |
| total_samples | int | Total samples in dataset |
| genes | string | Affected genes (optional) |

**Example**:
```csv
sv_id,chr,start,end,svtype,frequency,num_samples,total_samples,genes
chr7:55019032-55211628_DUP,chr7,55019032,55211628,DUP,0.45,150,333,EGFR
chr9:21967751-21994490_DEL,chr9,21967751,21994490,DEL,0.52,173,333,CDKN2A;CDKN2B
```

### Step 3: Run Comparison

```bash
# Make sure you've completed Steps 1-3 of main pipeline first
conda activate svmeta_env

# Run external comparison
python 04_external_dataset_comparison.py
```

## Output Files

### 1. Per-Dataset Comparisons

**tcga_comparison.csv** and **pcawg_comparison.csv**:

Columns:
- `cohort_sv`: Your SV ID
- `cohort_freq`: Frequency in your cohort
- `cohort_samples`: Number of samples in your cohort
- `external_match`: "Yes" if matched in external dataset
- `external_freq`: Frequency in external dataset
- `overlap_pct`: Percent reciprocal overlap
- `cohort_specific`: True if not found in external dataset
- `enriched_in_cohort`: True if >1.5x more frequent in cohort
- `enriched_in_external`: True if >1.5x more frequent externally

### 2. Universal Patterns

**universal_patterns.csv**:

Lists SV patterns shared across all datasets:
- Known GBM driver genes affected
- SV type distributions
- Common structural rearrangements

### 3. Visualizations

**frequency_comparison_TCGA-GBM.png**:
- Scatter plot: External frequency (x-axis) vs Cohort frequency (y-axis)
- Red points: Enriched in cohort
- Blue points: Enriched in external dataset
- Gray points: Similar frequency
- Diagonal line: Equal frequency (y=x)

**external_comparison_summary.png**:
- Bar charts showing:
  - Cohort-specific vs shared SV counts
  - Enrichment patterns across datasets

### 4. Text Report

**external_comparison_report.txt**:

Comprehensive summary including:
- Summary statistics (shared vs cohort-specific counts)
- Top cohort-specific SVs
- Universal patterns identified
- Gene-level comparisons

## Interpretation Guide

### What to Look For

#### 1. High Concordance (Expected)

If your cohort shows:
- EGFR amplifications in ~40-50% samples
- CDKN2A/B deletions in ~50-70% samples
- PTEN deletions in ~30-50% samples

**Interpretation**: Good! Your SV calling is detecting known GBM events.

#### 2. Cohort-Specific SVs

SVs found in your cohort but not in TCGA/PCAWG:

**Possible reasons**:
- Novel biological finding (exciting!)
- Technical differences (calling method, sequencing platform)
- Population-specific variations
- False positives (check manual if very frequent)

**Action**: Investigate top cohort-specific SVs:
- Are they in known cancer genes?
- Do they cluster in hotspots?
- Are they recurrent in your cohort?

#### 3. Missing Known SVs

Expected SVs not found in your cohort:

**Possible reasons**:
- Smaller sample size (some rare events not captured)
- Different SV calling sensitivity
- True biological differences in your cohort

**Action**: Check if major drivers (EGFR, CDKN2A, PTEN) are present at reasonable frequencies.

### Enrichment Analysis

**Enriched in cohort** (>1.5x higher frequency):
- May indicate selection in your cohort
- Could be technical artifacts if extreme
- Worth investigating if in driver genes

**Enriched in external** (>1.5x lower in cohort):
- May be under-called in your data
- Could be genuine biological difference
- Check SV calling sensitivity

## Advanced Usage

### Compare Specific Gene Sets

```python
import pandas as pd

# Load comparisons
tcga = pd.read_csv("results/external_comparison/tcga_comparison.csv")

# Filter to specific genes
driver_genes = ['EGFR', 'PTEN', 'CDKN2A', 'CDKN2B', 'PDGFRA']

# Extract gene-affected SVs (requires parsing sv_id or using gene annotations)
# ... custom filtering logic
```

### Calculate Concordance Metrics

```python
import pandas as pd
from scipy.stats import spearmanr

tcga = pd.read_csv("results/external_comparison/tcga_comparison.csv")
matched = tcga[tcga['external_match'] == 'Yes']

# Correlation of frequencies
corr, pval = spearmanr(matched['cohort_freq'], matched['external_freq'])
print(f"Frequency correlation: r={corr:.3f}, p={pval:.2e}")

# Percentage of shared SVs
pct_shared = (matched.shape[0] / tcga.shape[0]) * 100
print(f"Percentage of SVs shared with TCGA: {pct_shared:.1f}%")
```

### Identify High-Confidence Novel SVs

```python
import pandas as pd

tcga = pd.read_csv("results/external_comparison/tcga_comparison.csv")

# Cohort-specific + recurrent = high confidence novel
novel = tcga[
    (tcga['cohort_specific'] == True) &
    (tcga['cohort_samples'] >= 5)  # Present in 5+ samples
].sort_values('cohort_freq', ascending=False)

print("Top novel recurrent SVs:")
print(novel.head(10))
```

## Troubleshooting

### No External Data Available

If you can't access TCGA/PCAWG data:

1. Use the included literature reference file:
   ```bash
   ./download_external_datasets.sh  # Creates known_gbm_svs_literature.csv
   ```

2. Manually create reference from papers:
   - TCGA 2013 GBM paper (Brennan et al., Cell)
   - TCGA marker papers
   - Published GBM reviews

3. Focus on gene-level validation:
   - Check if known GBM genes are affected
   - Compare against gene-level recurrence literature

### Low Match Rate

If you see <20% of your SVs matching external datasets:

**Possible causes**:
1. Different SV size ranges (check size distributions)
2. Different calling sensitivity (common with long-read vs short-read)
3. Strict matching criteria (try reducing MIN_OVERLAP_BP in script)

**Solutions**:
- Relax matching parameters (reduce from 50% to 30% reciprocal overlap)
- Compare gene-level patterns instead of exact SV positions
- Focus on known driver genes

### Template Files Not Created

```bash
# Manually create external_datasets directory
mkdir -p external_datasets

# Run download script
./download_external_datasets.sh

# Check files
ls external_datasets/
```

## Known GBM SV Patterns (Reference)

Use this as a validation checklist:

### Core GBM SVs (Should see in most cohorts)

| Gene/Region | SV Type | Expected Frequency | Function |
|-------------|---------|-------------------|----------|
| EGFR (chr7:55-56M) | Amplification | 40-57% | RTK signaling |
| CDKN2A/B (chr9p21) | Deletion | 49-61% | Cell cycle |
| PTEN (chr10q23) | Deletion | 30-41% | PI3K pathway |
| TP53 (chr17p13) | Deletion | 25-35% | Tumor suppressor |
| RB1 (chr13q14) | Deletion | 7-11% | Cell cycle |

### Secondary GBM SVs (Common but variable)

| Gene/Region | SV Type | Expected Frequency | Function |
|-------------|---------|-------------------|----------|
| PDGFRA (chr4q12) | Amplification | 10-15% | RTK signaling |
| CDK4 (chr12q14) | Amplification | 15-18% | Cell cycle |
| MDM2/MDM4 (chr12q15) | Amplification | 10-15% | p53 regulation |
| MET (chr7q31) | Amplification | 4-8% | RTK signaling |
| NF1 (chr17q11) | Deletion | 10-18% | RAS pathway |

### Chromosomal Events

- Chromosome 7 gain: 50-70%
- Chromosome 10 loss: 60-80%
- Chromothripsis: 5-10% (complex rearrangements)

## References

1. **TCGA-GBM (2013)**: Brennan et al., "The Somatic Genomic Landscape of Glioblastoma", Cell
2. **TCGA Marker Paper (2008)**: Comprehensive genomic characterization, Nature
3. **PCAWG (2020)**: Pan-cancer structural variant analysis, Nature
4. **GBM Subtypes**: Verhaak et al. (2010), Cancer Cell

## Support

If you encounter issues:

1. Check input file formats match specifications
2. Verify Step 3 completed successfully (recurrent_svs.csv exists)
3. Review error messages in script output
4. Start with literature reference file before trying TCGA/PCAWG data

---

**Next Steps**: After comparison, use findings to:
1. Validate your cohort represents typical GBM
2. Identify novel therapeutic targets (cohort-specific SVs in druggable genes)
3. Contextualize findings in publications
4. Plan functional validation experiments
