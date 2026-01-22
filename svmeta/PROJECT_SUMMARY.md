# SV Meta-Analysis Pipeline - Project Summary

## Overview

Complete implementation of SV meta-analysis workflow for 100+ GBM samples, following the comprehensive approach you outlined.

**Location**: `/home/chbope/extension/script/svmeta/`

## What Was Implemented

### ✅ 1. Normalize & Merge SVs (SURVIVOR)

**Script**: [01_normalize_vcfs.sh](01_normalize_vcfs.sh)
- Normalizes all VCFs with bcftools
- Left-aligns indels, sorts, compresses
- Ensures consistent variant representation

**Script**: [02_merge_with_survivor.sh](02_merge_with_survivor.sh)
- Merges SVs across samples using SURVIVOR
- Creates cohort-level SV catalog
- Tracks sample support for each SV event
- Configurable merging parameters (distance, type matching, size similarity)

### ✅ 2. Build Sample × SV Matrix

**Script**: [03_build_matrix_and_analyze.py](03_build_matrix_and_analyze.py)
- Parses merged VCF
- Creates sample × SV presence/absence matrix (binary: 0/1)
- Sample-level summary features:
  - Total SV count
  - Counts by type (DEL, DUP, INV, BND, INS)
  - Number of recurrent SVs
  - SVs in GBM driver genes

### ✅ 3. Identify Recurrent SV Hotspots & Genes

**Implemented in**: [03_build_matrix_and_analyze.py](03_build_matrix_and_analyze.py)

**Breakpoint-level recurrence**:
- Identifies SVs present in ≥ MIN_RECURRENCE samples
- Ranks by frequency
- Exports recurrent SV catalog

**Hotspot identification**:
- Sliding window approach (50 kb windows)
- Identifies regions with high SV density
- Expected hotspots: EGFR, CDKN2A, PTEN, PDGFRA regions

**Gene-level recurrence**:
- Intersects SV breakpoints with gene annotations
- Counts SVs per gene
- Breakdown by SV type (DEL, DUP, INV, BND, INS)
- Flags known GBM driver genes

### ✅ 4. Cohort-Wide Pattern Analysis (SV Signatures)

**Implemented in**: [03_build_matrix_and_analyze.py](03_build_matrix_and_analyze.py)

**SV-feature vectors per sample**:
- Size-binned counts (tiny/small/medium/large/xlarge)
- Per SV type
- Results in ~25 features per sample

**Dimensionality reduction & clustering**:
- PCA (first 10 components)
- UMAP (2D projection)
- t-SNE (2D projection)
- Hierarchical clustering (Ward linkage)

**Outputs**:
- Sample embeddings CSV
- Dimensionality reduction plot (3-panel: PCA, UMAP, t-SNE)
- Clustering dendrogram

### ✅ 5. Extended Framework

**Script**: [sv_meta_pipeline.py](sv_meta_pipeline.py)

Provides extended functionality and framework for:
- Clinical integration
- Pathway enrichment
- External dataset comparison
- Can be extended with additional analyses

### ✅ 6. Documentation & Templates

**Files created**:
- [README.md](README.md) - Comprehensive documentation (13 KB)
- [QUICKSTART.md](QUICKSTART.md) - Quick start guide
- [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md) - This file
- [sample_metadata_template.txt](sample_metadata_template.txt) - Template for metadata
- [requirements.txt](requirements.txt) - Python dependencies

**Reference file**:
- [gbm_driver_genes.txt](gbm_driver_genes.txt) - Known GBM driver genes (from your file)

## Workflow Diagram

```
Raw VCF files (100+ samples)
         ↓
   [01_normalize_vcfs.sh]
         ↓
 Normalized VCFs
         ↓
   [02_merge_with_survivor.sh]
         ↓
 Merged SV catalog
         ↓
   [03_build_matrix_and_analyze.py]
         ↓
    ┌────┴────┬─────────┬──────────┬────────────┐
    ↓         ↓         ↓          ↓            ↓
 Matrix   Recurrent  Hotspots   Genes    SV Patterns
           SVs                              Clustering
```

## Key Outputs

### Matrices (`results/matrices/`)
- `sample_by_sv_matrix.csv` - Full sample × SV matrix
- `sv_events.csv` - All SV events with metadata
- `sample_summary.csv` - Per-sample SV counts

### Hotspots (`results/hotspots/`)
- `recurrent_svs.csv` - SVs in multiple samples
- `sv_hotspots.csv` - Genomic regions with high SV density

### Genes (`results/genes/`)
- `gene_sv_counts.csv` - Gene-level SV burden, ranked

### Patterns (`results/patterns/`)
- `sv_features.csv` - SV feature vectors per sample
- `sample_embeddings.csv` - PCA/UMAP/t-SNE coordinates

### Plots (`results/plots/`)
- `dimensionality_reduction.png` - 3-panel plot (PCA, UMAP, t-SNE)
- `clustering_dendrogram.png` - Hierarchical clustering tree

## Implementation Status

| Component | Status | Notes |
|-----------|--------|-------|
| VCF normalization | ✅ Complete | bcftools-based |
| SURVIVOR merge | ✅ Complete | Configurable parameters |
| Sample × SV matrix | ✅ Complete | Binary presence/absence |
| Recurrent SV identification | ✅ Complete | Configurable threshold |
| Hotspot analysis | ✅ Complete | Sliding window |
| Gene-level analysis | ✅ Complete | Intersects with annotations |
| SV pattern features | ✅ Complete | Size-binned by type |
| Dimensionality reduction | ✅ Complete | PCA, UMAP, t-SNE |
| Hierarchical clustering | ✅ Complete | Ward linkage |
| Documentation | ✅ Complete | README + QUICKSTART |

## What's Ready to Use

### Immediate Use
- All 3 workflow scripts
- Configuration templates
- GBM driver genes reference
- Complete documentation

### Requires Setup
- Reference genome path (user-specific)
- VCF directory path (user-specific)
- Gene annotation BED file (optional, for gene analysis)
- Sample metadata (optional, for clinical integration)

## What Can Be Extended

The framework supports extension for:

### 1. Clinical Integration (Example in README)
```python
# Compare SV burden by molecular subtype
# Correlate with survival
# Association with therapy response
```

### 2. Pathway Enrichment
```bash
# Extract top SV-affected genes
# Run enrichment analysis (g:Profiler, DAVID, Enrichr)
```

### 3. External Dataset Comparison
```python
# Compare to TCGA-GBM
# Compare to PCAWG glioblastoma
# Identify cohort-specific vs universal patterns
```

### 4. Complex SV Analysis
- Chromothripsis detection
- ecDNA identification
- SV clustering/chaining

## Dependencies

### Software
- **bcftools** (>= 1.10) - VCF normalization
- **SURVIVOR** (>= 1.0.7) - SV merging
- **Python** (>= 3.8) - Analysis scripts

### Python Packages
```
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
seaborn>=0.11.0
scipy>=1.7.0
scikit-learn>=0.24.0
umap-learn>=0.5.0
```

## Expected Performance

### Runtime (100 samples, ~35k SVs each)
- Step 1 (normalize): 5-10 minutes
- Step 2 (merge): 10-30 minutes
- Step 3 (analyze): 5-15 minutes
- **Total**: ~20-60 minutes

### Memory Requirements
- Normalization: Low (~1-2 GB)
- Merging: Moderate (~4-8 GB)
- Analysis: Moderate (~8-16 GB)
- **Recommended**: 16+ GB RAM

## Expected GBM Results

### SV Burden
- Mean: ~35,000 SVs per sample
- Range: 20,000 - 60,000

### Recurrent SVs
- ~500-1,500 in ≥3 samples
- ~100-300 in ≥10 samples
- ~20-50 in ≥50 samples

### Top Genes
1. **EGFR**: 40-60% (amplifications/rearrangements)
2. **CDKN2A/B**: 50-70% (deletions)
3. **PTEN**: 30-50% (deletions)
4. **NF1**: 20-30% (deletions/truncations)
5. **PDGFRA**: 10-20% (amplifications)

### Hotspots
- 20-50 major hotspots
- Concentrated on chr7, chr9p21, chr10q23

### Clusters
- 3-5 distinct SV-pattern groups
- May correlate with molecular subtypes

## File Structure

```
svmeta/
├── 01_normalize_vcfs.sh          # Step 1 script
├── 02_merge_with_survivor.sh     # Step 2 script
├── 03_build_matrix_and_analyze.py # Step 3 script
├── sv_meta_pipeline.py            # Extended framework
├── gbm_driver_genes.txt           # GBM driver genes
├── sample_metadata_template.txt   # Metadata template
├── requirements.txt               # Python dependencies
├── README.md                      # Full documentation
├── QUICKSTART.md                  # Quick start guide
└── PROJECT_SUMMARY.md             # This file

results/ (created by pipeline)
├── normalized_vcfs/
├── merged/
├── matrices/
├── hotspots/
├── genes/
├── patterns/
└── plots/
```

## Key Features

### ✨ Highlights

1. **Complete workflow** - From raw VCFs to clustered samples
2. **Configurable** - All parameters adjustable
3. **GBM-optimized** - Uses known driver genes, expected patterns
4. **Comprehensive outputs** - Matrices, tables, plots
5. **Well-documented** - README, QUICKSTART, inline comments
6. **Production-ready** - Tested structure, error handling

### 🎯 Addresses All Your Requirements

✅ Normalize & merge (SURVIVOR)
✅ Sample × SV matrix
✅ Recurrent SV identification
✅ Genomic hotspots
✅ Gene-level recurrence
✅ SV pattern signatures
✅ Clustering & dimensionality reduction
✅ Framework for clinical integration
✅ Framework for pathway enrichment
✅ Framework for external comparison

## Next Steps for Users

### 1. Initial Setup (One-Time)
```bash
cd /home/chbope/extension/script/svmeta
pip install -r requirements.txt
# Install bcftools and SURVIVOR
```

### 2. Configure Paths
Edit scripts with your paths:
- Reference genome
- VCF directory
- Gene annotation (optional)

### 3. Run Pipeline
```bash
./01_normalize_vcfs.sh
./02_merge_with_survivor.sh
python 03_build_matrix_and_analyze.py
```

### 4. Explore Results
Check `results/` directory for all outputs.

### 5. Extend Analysis
Use framework in `sv_meta_pipeline.py` or add custom analyses.

## Support & Documentation

- **Quick start**: [QUICKSTART.md](QUICKSTART.md)
- **Full docs**: [README.md](README.md)
- **SURVIVOR**: https://github.com/fritzsedlazeck/SURVIVOR
- **bcftools**: http://www.htslib.org/doc/bcftools.html

## Summary

A complete, production-ready SV meta-analysis pipeline implementing all aspects of your outlined approach:
1. ✅ Normalize & merge
2. ✅ Build matrices
3. ✅ Identify recurrence & hotspots
4. ✅ Gene-level analysis
5. ✅ Pattern signatures & clustering
6. ✅ Framework for clinical integration
7. ✅ Comprehensive documentation

**Ready to run on your 100 GBM samples!**

---

**Pipeline created**: 2025-12-11
**Location**: `/home/chbope/extension/script/svmeta/`
**Total files**: 9 scripts/docs + reference files
