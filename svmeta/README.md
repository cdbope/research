# SV Meta-Analysis Pipeline for GBM Cohort (svmeta)

Comprehensive pipeline for analyzing structural variants across 100+ GBM samples with focus on recurrence, hotspots, and patterns.

## Overview

This pipeline implements the complete meta-analysis workflow:

**Step 1**: Normalize & Merge SVs → cohort-level SV catalog
**Step 2**: Build sample × SV matrix → presence/absence matrix
**Step 3**: Identify recurrent SVs & hotspots → genomic patterns
**Step 4**: Gene-level analysis → driver gene identification
**Step 5**: SV pattern signatures → sample clustering
**Step 6**: Clinical integration → associations with outcomes

## Quick Start

### Installation (Automated Setup)

**Recommended**: Use the automated setup script to create the conda environment and install all dependencies:

```bash
cd /home/chbope/extension/script/svmeta

# Run setup script (installs conda env + SURVIVOR)
./00_setup_environment.sh

# Activate environment
conda activate svmeta_env
```

The setup script will:
- ✅ Create `svmeta_env` conda environment with all Python packages
- ✅ Install bcftools and samtools
- ✅ Download and build SURVIVOR
- ✅ Update script paths automatically

### Alternative: Manual Installation

If you prefer manual installation:

```bash
# 1. Create conda environment
cd /home/chbope/extension/script/svmeta
conda env create -f environment.yml
conda activate svmeta_env

# 2. Install SURVIVOR manually
git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug && make
export PATH=$PATH:$(pwd)
cd ../../
```

### Configure Paths

Edit the configuration sections in each script to match your system:

**01_normalize_vcfs.sh**:
```bash
REF_FASTA="/path/to/reference/genome.fa"
INPUT_VCF_DIR="/path/to/raw/vcf/files"
```

**03_build_matrix_and_analyze.py**:
```python
GENE_ANNOTATION = "/path/to/gencode_genes.bed"
```

### Run Pipeline

```bash
# Step 1: Normalize VCFs (5-10 min for 100 samples)
./01_normalize_vcfs.sh

# Step 2: Merge with SURVIVOR (10-30 min)
./02_merge_with_survivor.sh

# Step 3: Build matrix and analyze (5-15 min)
python 03_build_matrix_and_analyze.py

# Step 4 (Optional): Compare with external datasets
./download_external_datasets.sh  # Get external data (one-time)
python 04_external_dataset_comparison.py
```

## Files in This Directory

| File | Description |
|------|-------------|
| [00_setup_environment.sh](00_setup_environment.sh) | Automated installation of conda env and SURVIVOR |
| [01_prepare_vcfs.sh](01_prepare_vcfs.sh) | Prepare VCFs for SURVIVOR (decompress) |
| [02_merge_with_survivor.sh](02_merge_with_survivor.sh) | Merge SVs across samples using SURVIVOR |
| [03_build_matrix_and_analyze.py](03_build_matrix_and_analyze.py) | Build matrix, identify hotspots, analyze genes, cluster samples |
| [04_external_dataset_comparison.py](04_external_dataset_comparison.py) | Compare cohort with TCGA-GBM and PCAWG datasets |
| [download_external_datasets.sh](download_external_datasets.sh) | Guide for obtaining external reference datasets |
| [sv_meta_pipeline.py](sv_meta_pipeline.py) | Python framework (extended functionality) |
| [gbm_driver_genes.txt](gbm_driver_genes.txt) | Known GBM driver genes |
| [sample_metadata_template.txt](sample_metadata_template.txt) | Template for sample metadata |
| [environment.yml](environment.yml) | Conda environment definition |
| [requirements.txt](requirements.txt) | Python dependencies |
| README.md | This file |

## Output Structure

```
results/
├── prepared_vcfs/            # Step 1 output
│   ├── Sample001.wf_sv.vcf
│   └── ...
├── merged/                   # Step 2 output
│   ├── merged_SV.vcf.gz
│   └── vcf_list.txt
├── matrices/                 # Step 3: matrices
│   ├── sample_by_sv_matrix.csv
│   ├── sv_events.csv
│   └── sample_summary.csv
├── hotspots/                 # Step 3: recurrence
│   ├── recurrent_svs.csv
│   └── sv_hotspots.csv
├── genes/                    # Step 3: gene analysis
│   └── gene_sv_counts.csv
├── patterns/                 # Step 3: SV features
│   ├── sv_features.csv
│   └── sample_embeddings.csv
├── plots/                    # Step 3: visualizations
│   ├── dimensionality_reduction.png
│   └── clustering_dendrogram.png
└── external_comparison/      # Step 4: external datasets
    ├── tcga_comparison.csv
    ├── pcawg_comparison.csv
    ├── universal_patterns.csv
    ├── frequency_comparison_TCGA-GBM.png
    ├── frequency_comparison_PCAWG.png
    ├── external_comparison_summary.png
    └── external_comparison_report.txt
```

## Workflow Details

### Step 1: VCF Preparation

**Purpose**: Prepare VCF files for SURVIVOR merge

**Script**: `01_prepare_vcfs.sh`

**What it does**:
- Decompress VCF.gz files to plain VCF (SURVIVOR requires uncompressed)
- Create symbolic links to prepared VCF directory
- Validate files are non-empty

**Input**: Raw VCF files from Sniffles (*.vcf.gz)

**Output**: Prepared VCFs in `results/prepared_vcfs/`

**Key parameters**:
```bash
INPUT_VCF_DIR="/path/to/vcfs"            # Raw VCF directory
OUTPUT_VCF_DIR="results/prepared_vcfs"   # Output directory
```

### Step 2: SURVIVOR Merge

**Purpose**: Create cohort-level merged SV catalog

**Script**: `02_merge_with_survivor.sh`

**What it does**:
- Clusters similar SVs across samples into "meta-events"
- Each merged SV records which samples support it
- Creates unified catalog for cohort analysis

**Input**: Normalized VCFs from Step 1

**Output**: `results/merged/merged_SV.vcf.gz`

**Key parameters**:
```bash
MAX_DIST=1000          # Max distance to merge SVs (bp)
MIN_SUPPORT=1          # Min samples required
TYPE_MATCH=1           # Require same SV type (1=yes)
SIZE_SIMILARITY=0.3    # Size similarity threshold (0-1)
```

**Parameter tuning**:
- **Tight merging** (high confidence): MAX_DIST=500, SIZE_SIMILARITY=0.5
- **Loose merging** (high sensitivity): MAX_DIST=5000, SIZE_SIMILARITY=0.1

### Step 3: Matrix & Analysis

**Purpose**: Comprehensive SV meta-analysis

**Script**: `03_build_matrix_and_analyze.py`

**What it does**:

1. **Build sample × SV matrix**
   - Rows: samples
   - Columns: SV events
   - Values: 0/1 (presence/absence)

2. **Identify recurrent SVs**
   - SVs present in ≥3 samples (configurable)
   - Sorted by recurrence frequency

3. **Find genomic hotspots**
   - Sliding window approach (50 kb windows)
   - Regions with high SV density across samples

4. **Gene-level analysis**
   - Intersect SV breakpoints with genes
   - Count SVs affecting each gene
   - Flag known GBM driver genes

5. **Compute SV pattern features**
   - Size-binned SV counts (tiny/small/medium/large/xlarge)
   - By SV type (DEL, DUP, INV, BND, INS)
   - Feature vector per sample

6. **Cluster samples**
   - PCA, UMAP, t-SNE dimensionality reduction
   - Hierarchical clustering
   - Identify SV-pattern-based subtypes

**Input**:
- Merged VCF from Step 2
- Gene annotation (BED format)
- GBM driver genes list

**Output**: Multiple files (see Output Structure above)

**Key parameters**:
```python
MIN_RECURRENCE = 3           # Min samples for "recurrent"
HOTSPOT_WINDOW = 50000       # Hotspot window size (bp)
MIN_HOTSPOT_SAMPLES = 5      # Min samples per hotspot

SIZE_BINS = {
    'tiny': (50, 1000),
    'small': (1000, 10000),
    'medium': (10000, 100000),
    'large': (100000, 1000000),
    'xlarge': (1000000, float('inf'))
}
```

### Step 4: External Dataset Comparison (Optional)

**Purpose**: Compare cohort SV patterns against reference datasets

**Script**: `04_external_dataset_comparison.py`

**What it does**:

1. **Compare with TCGA-GBM**
   - Match SVs by genomic position and type
   - Compare frequencies
   - Identify cohort-specific vs shared SVs
   - Detect enriched/depleted SVs

2. **Compare with PCAWG**
   - Same analysis as TCGA
   - Broader pan-cancer perspective

3. **Identify universal patterns**
   - SVs common across all datasets
   - Known GBM driver alterations
   - Validate findings against literature

4. **Generate comparison reports**
   - Frequency scatter plots
   - Venn diagram summaries
   - Enrichment statistics
   - Text reports

**Setup**:
```bash
# Step 1: Run download guide
./download_external_datasets.sh

# Step 2: Place external datasets in external_datasets/
#   - tcga_gbm_sv_summary.csv
#   - pcawg_gbm_sv_summary.csv
# OR use provided templates and populate with literature data

# Step 3: Run comparison
python 04_external_dataset_comparison.py
```

**Input**:
- Cohort results from Step 3 (recurrent_svs.csv, gene_recurrence.csv)
- External dataset files (TCGA/PCAWG summaries)

**Output**:
- `external_comparison/tcga_comparison.csv`: Per-SV TCGA comparison
- `external_comparison/pcawg_comparison.csv`: Per-SV PCAWG comparison
- `external_comparison/universal_patterns.csv`: Shared patterns
- `external_comparison/frequency_comparison_*.png`: Scatter plots
- `external_comparison/external_comparison_summary.png`: Overview
- `external_comparison/external_comparison_report.txt`: Text report

**Key features**:
- **Reciprocal overlap matching**: SVs matched based on 50%+ reciprocal overlap
- **Frequency enrichment**: 1.5x fold-change threshold for enrichment
- **Literature validation**: Known GBM SVs from TCGA 2013 included
- **Visual comparisons**: Frequency scatter plots with enrichment zones

**Expected findings**:
- EGFR amplifications: Shared with TCGA (~40-50%)
- CDKN2A/B deletions: Shared with TCGA (~50-70%)
- PTEN deletions: Shared with TCGA (~30-50%)
- Cohort-specific SVs: Novel findings or technical differences
- Universal patterns: Core GBM SV signature

## Key Analyses

### 1. Recurrent SV Identification

**File**: `results/hotspots/recurrent_svs.csv`

Identifies SVs present in multiple samples:

```csv
sv_id,chr,pos,end,svtype,svlen,n_samples
SV0042,chr7,55000000,55100000,DEL,100000,45
SV0128,chr10,89000000,89500000,DUP,500000,38
```

**Interpretation**:
- **n_samples ≥ 10**: Highly recurrent, likely functional
- **n_samples 3-9**: Moderately recurrent
- Focus on recurrent SVs in known GBM genes

**Expected GBM recurrent regions**:
- **chr7:55-56 Mb**: EGFR amplifications (~40-60% samples)
- **chr9p21**: CDKN2A/B deletions (~50-70%)
- **chr10q23**: PTEN deletions (~30-50%)
- **chr12q13-15**: CDK4/MDM2 amplifications (~10-20%)

### 2. Genomic Hotspot Analysis

**File**: `results/hotspots/sv_hotspots.csv`

Identifies regions with clustered SVs across samples:

```csv
chr,start,end,n_svs,total_samples,mean_samples
chr7,55000000,55050000,15,125,8.3
chr9,21900000,21950000,12,98,8.2
```

**Interpretation**:
- High `total_samples`: Cohort-wide rearrangement
- Hotspots near known drivers are expected
- Novel hotspots may indicate new candidate genes

### 3. Gene-Level SV Burden

**File**: `results/genes/gene_sv_counts.csv`

Ranks genes by SV burden across cohort:

```csv
gene,n_svs,n_samples,DEL,DUP,INV,BND,INS,is_driver
EGFR,85,52,5,45,10,25,0,TRUE
CDKN2A,72,68,70,0,0,2,0,TRUE
PTEN,58,42,48,2,3,5,0,TRUE
```

**Expected top genes (GBM)**:
1. **EGFR**: Amplifications (DUP), rearrangements (BND)
2. **CDKN2A/CDKN2B**: Homozygous deletions (DEL)
3. **PTEN**: Deletions (DEL)
4. **NF1**: Deletions, truncations (DEL, BND)
5. **PDGFRA**: Amplifications (DUP)

### 4. SV Pattern Clustering

**Files**:
- `results/patterns/sample_embeddings.csv`: PCA/UMAP/t-SNE coordinates
- `results/plots/clustering_dendrogram.png`: Hierarchical clustering tree

**Interpretation**:
- **Tight clusters**: Samples with similar SV patterns
- **Outliers**: Unusual SV profiles (chromothripsis, high burden)
- **Gradients**: Continuous variation in SV load

**May correlate with**:
- GBM molecular subtypes (Classical, Mesenchymal, Proneural)
- IDH status, MGMT methylation
- Clinical outcomes

## Advanced Analysis Examples

### Compare SV Burden by IDH Status

```python
import pandas as pd
from scipy.stats import mannwhitneyu

# Load data
summary = pd.read_csv("results/matrices/sample_summary.csv")
metadata = pd.read_csv("sample_metadata.txt", sep='\t')

# Merge
df = summary.merge(metadata, left_on='sample', right_on='sample_id')

# Compare
idh_wt = df[df['IDH_status'] == 'WT']['total_svs']
idh_mut = df[df['IDH_status'] == 'mutant']['total_svs']

stat, pval = mannwhitneyu(idh_wt, idh_mut)
print(f"SV burden: IDH-WT vs IDH-mutant, p={pval:.4f}")
```

### Identify Samples with High Recurrent SV Load

```python
import pandas as pd

summary = pd.read_csv("results/matrices/sample_summary.csv")

# Top 10 samples by recurrent SVs
top_recurrent = summary.nlargest(10, 'recurrent_svs')
print(top_recurrent[['sample', 'total_svs', 'recurrent_svs']])
```

### Pathway Enrichment for SV-Affected Genes

```bash
# Extract top 50 genes
head -51 results/genes/gene_sv_counts.csv | tail -50 | cut -d',' -f1 > top_sv_genes.txt

# Use with g:Profiler, DAVID, Enrichr, etc.
```

**Expected GBM pathways**:
- RTK/RAS/PI3K signaling
- p53 pathway
- Cell cycle regulation
- DNA repair
- Chromatin remodeling

## Required Reference Files

### 1. Reference Genome

**What**: FASTA file used for SV calling

**Where**: Same reference used with Sniffles

**Index**:
```bash
samtools faidx genome.fa
```

### 2. Gene Annotation (BED format)

**What**: Gene coordinates for overlap analysis

**Format**:
```
chr1    11868   14409   DDX11L1    +
chr1    14361   29370   WASH7P     -
```

**How to create from GENCODE**:
```bash
# Download GENCODE GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip gencode.v38.annotation.gtf.gz

# Extract genes and convert to BED
grep -w gene gencode.v38.annotation.gtf | \
    awk 'OFS="\t" {print $1,$4,$5,$10,$7}' | \
    tr -d '";' > gencode_genes.bed
```

### 3. Sample Metadata (optional but recommended)

**What**: Clinical and molecular features per sample

**Format**: Tab-separated with header

**Required column**: `sample_id`

**Optional columns**: purity, age, sex, IDH_status, MGMT_status, survival_months, etc.

**Template**: See `sample_metadata_template.txt`

### 4. GBM Driver Genes

**What**: Known GBM driver genes for flagging

**Provided**: `gbm_driver_genes.txt` (based on your reference file)

**You can add more genes** to this file if needed.

## Expected Results (100 GBM Samples)

Based on typical GBM cohorts:

**SV Burden per Sample**:
- Mean: ~35,000 SVs
- Range: 20,000 - 60,000

**Recurrent SVs**:
- 500-1,500 SVs in ≥3 samples
- 100-300 SVs in ≥10 samples
- 20-50 SVs in ≥50 samples

**Top Recurrent Genes**:
1. EGFR: 40-60% samples
2. CDKN2A/B: 50-70%
3. PTEN: 30-50%
4. NF1: 20-30%
5. PDGFRA: 10-20%

**Genomic Hotspots**: 20-50 major hotspots

**Sample Clusters**: 3-5 distinct SV-pattern groups

## Troubleshooting

### bcftools norm fails

**Error**: `Failed to open reference`

**Solution**:
- Check REF_FASTA path
- Ensure reference is indexed: `samtools faidx genome.fa`

### SURVIVOR not found

**Error**: `SURVIVOR: command not found`

**Solution**:
```bash
export PATH=$PATH:/path/to/SURVIVOR/Debug
```

Or update `SURVIVOR_BIN` in script to full path.

### Gene annotation not found

**Warning**: `Gene annotation not found, skipping gene analysis`

**Solution**:
- Download and create BED file (see Required Files)
- Update `GENE_ANNOTATION` path in `03_build_matrix_and_analyze.py`

### Not enough memory

**Error**: Process killed during merge or clustering

**Solution**:
- Run on HPC/cluster with ≥32 GB RAM
- Or reduce cohort size for testing

## Parameter Tuning

### For High-Confidence Recurrent SVs

```bash
# In 02_merge_with_survivor.sh
MAX_DIST=500
SIZE_SIMILARITY=0.5

# In 03_build_matrix_and_analyze.py
MIN_RECURRENCE = 10
```

### For Sensitive Detection

```bash
# In 02_merge_with_survivor.sh
MAX_DIST=5000
SIZE_SIMILARITY=0.1

# In 03_build_matrix_and_analyze.py
MIN_RECURRENCE = 2
```

### For Fine-Grained Hotspots

```python
# In 03_build_matrix_and_analyze.py
HOTSPOT_WINDOW = 10000
MIN_HOTSPOT_SAMPLES = 3
```

## Citations

**SURVIVOR**:
- Jeffares DC et al. "Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast." Nat Commun. 2017.

**bcftools**:
- Danecek P et al. "Twelve years of SAMtools and BCFtools." GigaScience. 2021.

## Support

For issues with:
- **SURVIVOR**: https://github.com/fritzsedlazeck/SURVIVOR/issues
- **bcftools**: http://www.htslib.org/doc/bcftools.html
- **This pipeline**: Contact your bioinformatics team

## License

Research use only.

---

**Created for GBM SV meta-analysis | Last updated: 2025**
