# SV Analysis Project - Complete Documentation Index

## 🚀 Quick Start

**New to this project?** Start here:
1. Read [PORTABILITY_UPDATE.md](PORTABILITY_UPDATE.md) - Latest changes for easy configuration
2. Follow [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md) - How to set up paths
3. Run: `python sv_analysis_v2.py`

## 📋 Project Overview

This project analyzes structural variants (SVs) from Sniffles VCF files for GBM samples, with focus on **significant VAF metrics** and GBM-optimized clustering.

## 📁 File Structure

### Analysis Scripts

| File | Description | Recommended |
|------|-------------|-------------|
| [sv_analysis.py](sv_analysis.py) | Version 1: Basic analysis with all VAF | ❌ |
| [sv_analysis_v2.py](sv_analysis_v2.py) | Version 2: Significant VAF analysis | ✅ **USE THIS** |
| [run_analysis.sh](run_analysis.sh) | Quick runner script | ✅ |

**Recommendation**: Use `sv_analysis_v2.py` for all new analyses.

### Configuration & Setup

| File | Purpose |
|------|---------|
| [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md) | **Complete guide to configuring paths** |
| [PORTABILITY_UPDATE.md](PORTABILITY_UPDATE.md) | Summary of portability improvements |
| [requirements.txt](requirements.txt) | Python dependencies |
| [sample.txt](sample.txt) | Sample IDs and purity values (input) |

### Documentation - User Guides

| File | When to Read |
|------|--------------|
| [README.md](README.md) | Version 1 documentation (legacy) |
| [README_V2.md](README_V2.md) | **Version 2 comprehensive documentation** |
| [QUICK_START.md](QUICK_START.md) | Quick reference for Version 1 |
| [VAF_ANALYSIS_GUIDE.md](VAF_ANALYSIS_GUIDE.md) | **Deep dive into VAF methodology** |
| [SV_HISTOGRAM_GUIDE.md](SV_HISTOGRAM_GUIDE.md) | Understanding the aggregate histogram plot |

### Documentation - Project Info

| File | Purpose |
|------|---------|
| [FINAL_SUMMARY.md](FINAL_SUMMARY.md) | Complete project summary and results |
| [CHANGES.md](CHANGES.md) | Change log (latest updates) |
| [INDEX.md](INDEX.md) | This file - documentation index |

## 🎯 Reading Path by Task

### Task: I want to run the analysis

1. **First time**: Read [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md)
2. **Configure**: Edit `sv_analysis_v2.py` lines 28-59
3. **Run**: `python sv_analysis_v2.py`

### Task: I want to understand the methodology

1. Read [VAF_ANALYSIS_GUIDE.md](VAF_ANALYSIS_GUIDE.md) - VAF concepts
2. Read [README_V2.md](README_V2.md) - Full methodology
3. Read [SV_HISTOGRAM_GUIDE.md](SV_HISTOGRAM_GUIDE.md) - Visualization details

### Task: I want to move the analysis to a different server

1. Read [PORTABILITY_UPDATE.md](PORTABILITY_UPDATE.md)
2. Follow [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md)
3. Update paths in `sv_analysis_v2.py`

### Task: I want to understand the results

1. Read [README_V2.md](README_V2.md) - Output files section
2. Read [FINAL_SUMMARY.md](FINAL_SUMMARY.md) - Interpretation
3. Read [SV_HISTOGRAM_GUIDE.md](SV_HISTOGRAM_GUIDE.md) - Histogram interpretation

### Task: I want to publish the results

1. Read [VAF_ANALYSIS_GUIDE.md](VAF_ANALYSIS_GUIDE.md) - Methods section
2. Read [FINAL_SUMMARY.md](FINAL_SUMMARY.md) - Publication recommendations
3. Use examples from [README_V2.md](README_V2.md)

### Task: I want to customize quality filters

1. Read [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md) - Advanced section
2. Edit `sv_analysis_v2.py` lines 48-55
3. See [VAF_ANALYSIS_GUIDE.md](VAF_ANALYSIS_GUIDE.md) for rationale

## 📊 Output Files

Analysis results are saved in `results_v2/` directory (configurable):

### Tables
- `sv_summary_display.csv` - Display table (Sample, Purity, DEL, DUP, INV, BND, INS)
- `sv_summary_table.csv` - Full table with all VAF metrics (24+ columns)
- `clustering_results.csv` - Cluster assignments

### Plots (13 files)
- **Clustering**: 3 dendrograms (high_quality, clonal, weighted)
- **PCA**: Sample relationships in 2D
- **Heatmap**: Normalized features across samples
- **SV Distribution**: Per-sample SV types
- **SV Stacked Histogram**: Aggregate SV type counts
- **VAF Analysis**: VAF vs purity, distributions
- **Comprehensive VAF**: 6-panel comparison

### Reports
- `analysis_report.txt` - Statistical summary with cluster characteristics

## 🔑 Key Concepts

### SV Types
- **DEL**: Deletions
- **DUP**: Duplications
- **INV**: Inversions
- **BND**: Breakends/Translocations
- **INS**: Insertions

### VAF Metrics (Version 2)
- **Mean_VAF_HQ**: Mean VAF of high-quality variants (RECOMMENDED)
- **Mean_VAF_Clonal**: Mean VAF of clonal variants (VAF ≥ 0.3)
- **Clonal_Fraction**: Proportion of variants that are clonal
- **Weighted_VAF**: VAF weighted by variant frequency
- **Top10_VAF**: Mean of top 10% highest VAF

### Quality Filters
- **PASS**: Variant passed quality filters
- **PRECISE**: Precise breakpoints (not IMPRECISE)
- **Support ≥ 10**: At least 10 supporting reads
- **VAF ≥ 0.1**: Minimum 10% variant allele frequency

## 🎨 Version Comparison

| Feature | Version 1 | Version 2 |
|---------|-----------|-----------|
| VAF calculation | All variants | High-quality only |
| Quality filtering | Basic | Multi-tier (PASS+PRECISE+Support) |
| Clonality analysis | No | Yes (≥0.3 threshold) |
| Clustering strategies | 1 | 3 (high_quality, clonal, weighted) |
| Weighted VAF | No | Yes |
| Top percentile VAF | No | Yes (10%, 25%) |
| Stacked histogram | No | Yes (aggregate SV types) |
| **Recommended for GBM** | ❌ | ✅ |

## ⚙️ Configuration Parameters

Located in `sv_analysis_v2.py` lines 28-59:

### Paths (Required)
- `SAMPLE_FILE`: Path to sample.txt
- `VCF_DIR`: Directory with VCF files
- `OUTPUT_DIR`: Results output directory

### Quality Filters (Optional)
- `MIN_SUPPORT`: Minimum read support (default: 10)
- `MIN_VAF`: Minimum VAF threshold (default: 0.1)
- `CLONAL_THRESHOLD`: Clonality threshold (default: 0.3)

## 🚨 Common Issues

### "VCF file not found"
- **Solution**: Check `VCF_DIR` and filename pattern `{sample_id}.wf_sv.vcf.gz`
- **See**: [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md) - Troubleshooting

### "No such file: sample.txt"
- **Solution**: Check `SAMPLE_FILE` path
- **See**: [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md) - File Requirements

### "Permission denied" (output)
- **Solution**: Choose writable `OUTPUT_DIR`
- **See**: [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md) - Troubleshooting

## 📖 Detailed Documentation

### For Methodology Details
- [VAF_ANALYSIS_GUIDE.md](VAF_ANALYSIS_GUIDE.md) - 342 lines of VAF methodology
- [README_V2.md](README_V2.md) - 650+ lines of comprehensive documentation

### For Quick Reference
- [QUICK_START.md](QUICK_START.md) - Version 1 quick start
- [PORTABILITY_UPDATE.md](PORTABILITY_UPDATE.md) - Latest changes summary

### For Visualization Understanding
- [SV_HISTOGRAM_GUIDE.md](SV_HISTOGRAM_GUIDE.md) - Complete histogram guide

## 🔬 GBM-Specific Features

Version 2 is optimized for GBM analysis:

1. ✅ **Tumor purity adjustment** - Critical for infiltrative GBM
2. ✅ **High SV burden handling** - GBM has 35k+ SVs per sample
3. ✅ **Clonality assessment** - GBM typically monoclonal (>80%)
4. ✅ **Homozygous deletion detection** - Common in GBM (CDKN2A, PTEN)
5. ✅ **Ward linkage clustering** - Best for genomic subtypes

## 📝 Citation & Acknowledgment

When publishing, cite:
- Significant VAF methodology (high-quality filtering)
- Quality criteria: PASS + PRECISE + Support ≥10
- Clonality threshold: VAF ≥ 0.3
- GBM-optimized clustering approach

## 🎓 Learning Path

**Beginner**:
1. [PORTABILITY_UPDATE.md](PORTABILITY_UPDATE.md)
2. [CONFIGURATION_GUIDE.md](CONFIGURATION_GUIDE.md)
3. Run the analysis

**Intermediate**:
1. [README_V2.md](README_V2.md)
2. [SV_HISTOGRAM_GUIDE.md](SV_HISTOGRAM_GUIDE.md)
3. Explore results

**Advanced**:
1. [VAF_ANALYSIS_GUIDE.md](VAF_ANALYSIS_GUIDE.md)
2. [FINAL_SUMMARY.md](FINAL_SUMMARY.md)
3. Customize parameters

## 📦 Dependencies

See [requirements.txt](requirements.txt):
- pandas ≥ 1.3.0
- numpy ≥ 1.21.0
- matplotlib ≥ 3.4.0
- seaborn ≥ 0.11.0
- scipy ≥ 1.7.0
- scikit-learn ≥ 0.24.0

## 🎯 Project Status

✅ **Version 2 Complete** - Recommended for all GBM analyses
✅ **Fully Documented** - 14 documentation files
✅ **Portable** - Easy configuration for any system
✅ **Production Ready** - Publication-quality results

---

**Last Updated**: 2025 (Portability improvements)

**Maintainer**: See project repository

**Support**: Read the relevant documentation file for your task above
