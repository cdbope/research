# ONT Time-Series CNV Monitor Pipeline

A comprehensive pipeline for analyzing Copy Number Variations (CNVs) from Oxford Nanopore Technology (ONT) sequencing data generated at intervals using **CNVpytor**.

## Quick Start

### 1. Setup (First time only)

```bash
cd /home/chbope/extension/script/cnvmonitor

# Create conda environment (takes 10-30 minutes)
conda env create -f cnvgen_environment.yml

# Activate environment
conda activate cnvgen

# Verify installation
cnvpytor --version
samtools --version
```

### 2. Run Pipeline

**Option A: Interactive Menu (Recommended for beginners)**
```bash
./ont_cnv_quickstart.sh
```

**Option B: Continuous Monitoring (Real-time analysis)**
```bash
./ont_timeseries_cnv_monitor.sh \
  -i /path/to/bam_directory \
  -o /path/to/output \
  -r /path/to/reference.fa \
  -s SAMPLE_ID \
  -c 300 \
  -t 4
```

**Option C: Batch Processing (Process existing BAMs once)**
```bash
./ont_timeseries_cnv_monitor.sh \
  -i /path/to/bam_directory \
  -o /path/to/output \
  -r /path/to/reference.fa \
  -s SAMPLE_ID \
  --no-continuous \
  -t 8
```

## Pipeline Files

- **cnvgen_environment.yml** - Conda environment configuration
- **ont_process_individual_bam_cnv.sh** - Process individual BAM files
- **ont_timeseries_cnv_visualization.R** - Generate time-series visualizations
- **ont_combined_cnv_analysis.R** - Perform combined CNV analysis
- **ont_timeseries_cnv_monitor.sh** - Main orchestration script
- **ont_cnv_quickstart.sh** - Interactive quick-start menu
- **ONT_TIMESERIES_CNV_README.md** - Comprehensive documentation
- **ONT_CNV_PIPELINE_SUMMARY.txt** - Quick reference guide

## Output Structure

```
output_directory/
├── individual_bams/     # Individual timepoint CNV analyses
│   └── SAMPLE_*_timestamp/
│       ├── *.pytor      # CNVpytor data files
│       ├── *.vcf        # VCF format CNV calls
│       ├── *.tsv        # Tab-separated CNV calls
│       ├── *.png        # Manhattan plots
│       └── *_summary.txt
├── timeseries_plots/    # Time-series visualizations
│   ├── *_cnv_count_timeseries.png
│   ├── *_genome_wide_cnv_heatmap.png
│   ├── *_per_chromosome_cnv.png
│   └── *_cnv_size_distribution.png
├── combined_analysis/   # Cross-timepoint CNV analysis
│   ├── *_cnv_classification.png
│   ├── *_persistent_cnv_map.png
│   ├── *_cnv_evolution_tracking.png
│   ├── *_cnv_classification.csv
│   └── *_persistent_cnvs.csv
└── logs/               # Processing logs and tracking
    ├── *_processing_log.txt
    └── *_processed_bams.txt
```

## Features

✓ **Fast CNV detection** using CNVpytor (optimized for ONT long reads)
✓ **Automatic monitoring** for new BAM files
✓ **Multi-resolution analysis** with multiple bin sizes (1kb, 10kb, 100kb)
✓ **Time-series tracking** of CNV evolution
✓ **CNV classification** into persistent, emerging, and transient categories
✓ **Rich visualizations** with multiple plot types
✓ **Comprehensive reports** in text and CSV formats

## CNV Classifications

- **Persistent CNVs**: Present in ALL time points (likely germline or early somatic)
- **Emerging CNVs**: Present in MULTIPLE time points (clonal evolution)
- **Transient CNVs**: Present in ONE time point only (technical artifacts or transient events)

## Documentation

For detailed documentation, see:
- [ONT_TIMESERIES_CNV_README.md](ONT_TIMESERIES_CNV_README.md) - Full documentation with examples
- [ONT_CNV_PIPELINE_SUMMARY.txt](ONT_CNV_PIPELINE_SUMMARY.txt) - Quick reference

## Requirements

- **Memory**: 8-16 GB RAM recommended
- **Disk Space**: ~50-100% of total BAM size
- **Processing Time**: 5-60 minutes per BAM (depends on size)

## Support

For issues or questions:
- Check the logs in `output_directory/logs/`
- Review the comprehensive documentation in `ONT_TIMESERIES_CNV_README.md`
- Verify conda environment with `conda activate cnvgen`

## Citation

If you use this pipeline, please cite:
- **CNVpytor**: Suvakov et al. (2021) CNVpytor: a tool for CNV/CNA detection and analysis from read depth and allele imbalance in whole-genome sequencing. Bioinformatics.

---

**Created**: 2025-10-14
**Location**: /home/chbope/extension/script/cnvmonitor/
