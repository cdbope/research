# Quick Start Guide - SV Meta-Analysis Pipeline

**Location**: `/home/chbope/extension/script/svmeta`

## 3-Step Quick Start

### Step 1: Install Dependencies (One-Time Setup)

**Option A: Automated Setup (Recommended)**

```bash
cd /home/chbope/extension/script/svmeta

# Run automated setup
./00_setup_environment.sh

# Activate environment
conda activate svmeta_env
```

This will install everything: Python packages, bcftools, samtools, and SURVIVOR.

**Option B: Manual Setup**

```bash
cd /home/chbope/extension/script/svmeta

# Create conda environment
conda env create -f environment.yml
conda activate svmeta_env

# Install SURVIVOR
git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug && make
export PATH=$PATH:$(pwd)
cd ../../
```

### Step 2: Configure Paths

Edit these 3 paths in the scripts:

**In `01_normalize_vcfs.sh`**:
```bash
REF_FASTA="/path/to/reference/genome.fa"
INPUT_VCF_DIR="/path/to/your/vcf/files"
```

**In `03_build_matrix_and_analyze.py`** (line 22):
```python
GENE_ANNOTATION = "/path/to/gencode_genes.bed"
```

That's it! The rest is auto-configured.

### Step 3: Run Pipeline

**Make sure svmeta_env is activated first:**

```bash
# Activate environment (if not already active)
conda activate svmeta_env

# Run pipeline
# Step 1: Normalize VCFs (~5-10 min for 100 samples)
./01_normalize_vcfs.sh

# Step 2: Merge with SURVIVOR (~10-30 min)
./02_merge_with_survivor.sh

# Step 3: Analyze (~5-15 min)
python 03_build_matrix_and_analyze.py
```

**Tip**: You can also use the helper script:
```bash
source activate_svmeta.sh  # Activates environment and shows SURVIVOR path
```

## What You Get

After running, check `results/` directory:

```
results/
├── matrices/
│   ├── sample_by_sv_matrix.csv      # Sample × SV presence/absence matrix
│   ├── sv_events.csv                 # All SV events
│   └── sample_summary.csv            # Per-sample SV counts
│
├── hotspots/
│   ├── recurrent_svs.csv             # SVs in multiple samples
│   └── sv_hotspots.csv               # Genomic hotspot regions
│
├── genes/
│   └── gene_sv_counts.csv            # Gene-level SV burden
│
├── patterns/
│   ├── sv_features.csv               # SV feature vectors
│   └── sample_embeddings.csv         # PCA/UMAP/t-SNE coordinates
│
└── plots/
    ├── dimensionality_reduction.png  # PCA, UMAP, t-SNE plots
    └── clustering_dendrogram.png     # Hierarchical clustering
```

## Key Files to Check

### 1. Recurrent SVs
```bash
head results/hotspots/recurrent_svs.csv
```

Look for SVs with high `n_samples` (present in many samples).

### 2. Top Genes with SVs
```bash
head -20 results/genes/gene_sv_counts.csv
```

Expected top genes for GBM:
- EGFR (amplifications)
- CDKN2A, CDKN2B (deletions)
- PTEN (deletions)
- NF1, PDGFRA

### 3. Sample Clustering
```bash
open results/plots/clustering_dendrogram.png
```

Shows samples grouped by SV patterns.

## What If Things Don't Work?

### Missing bcftools
```bash
conda install -c bioconda bcftools
```

### Missing SURVIVOR
```bash
git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug && make
export PATH=$PATH:$(pwd)
```

### Missing gene annotation
Download GENCODE and convert to BED:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip gencode.v38.annotation.gtf.gz
grep -w gene gencode.v38.annotation.gtf | \
    awk 'OFS="\t" {print $1,$4,$5,$10,$7}' | \
    tr -d '";' > gencode_genes.bed
```

Then update path in `03_build_matrix_and_analyze.py`.

### Python import errors
```bash
pip install -r requirements.txt
```

## Next Steps

1. **View results**: Check all CSV files in `results/` subdirectories
2. **Explore clusters**: Look at `clustering_dendrogram.png`
3. **Identify drivers**: Check top genes in `gene_sv_counts.csv`
4. **Integrate clinical data**: Add your metadata and run associations

## Full Documentation

See [README.md](README.md) for complete documentation including:
- Detailed workflow explanation
- Parameter tuning guide
- Advanced analysis examples
- Troubleshooting
- Expected results for GBM

## File List

| File | Purpose |
|------|---------|
| `01_normalize_vcfs.sh` | Normalize VCFs with bcftools |
| `02_merge_with_survivor.sh` | Merge SVs across samples |
| `03_build_matrix_and_analyze.py` | Main analysis script |
| `sv_meta_pipeline.py` | Extended Python framework |
| `gbm_driver_genes.txt` | Known GBM driver genes |
| `sample_metadata_template.txt` | Template for metadata |
| `requirements.txt` | Python dependencies |
| `README.md` | Full documentation |
| `QUICKSTART.md` | This file |

## Support

For detailed information, see [README.md](README.md).

For SURVIVOR issues: https://github.com/fritzsedlazeck/SURVIVOR/issues

---

**Ready to go? Run the 3 steps above!**
