# ONT Time-Series CNV Analysis Pipeline

A comprehensive pipeline for analyzing Copy Number Variations (CNVs) from Oxford Nanopore Technology (ONT) sequencing data generated at intervals. This pipeline uses **CNVpytor** for fast and accurate CNV detection from long-read sequencing data.

## Overview

This pipeline automatically:
1. Monitors a directory for new BAM files generated from ONT basecalling
2. Processes each BAM file individually using CNVpytor for CNV detection
3. Generates time-series visualizations showing CNV evolution over time
4. Identifies persistent, emerging, and transient CNVs across all time points
5. Produces comprehensive reports and publication-quality figures

## Features

- **Fast CNV Detection**: Uses CNVpytor which is optimized for long-read data
- **Multi-resolution Analysis**: Supports multiple bin sizes (1kb, 10kb, 100kb)
- **Automatic Monitoring**: Continuously monitors for new BAM files
- **Time-Series Tracking**: Tracks CNV changes across multiple time points
- **CNV Classification**: Identifies persistent, emerging, and transient CNVs
- **Rich Visualizations**: Generates multiple plot types for comprehensive analysis
- **Modular Design**: Each component can be run independently

## Installation

### 1. Create Conda Environment

```bash
cd /home/chbope/extension/script
conda env create -f cnvgen_environment.yml
conda activate cnvgen
```

### 2. Verify Installation

```bash
# Check CNVpytor
cnvpytor --version

# Check R packages
Rscript -e "library(ggplot2); library(dplyr)"

# Check samtools
samtools --version
```

## Pipeline Components

### 1. Individual BAM Processing Script
**File**: `ont_process_individual_bam_cnv.sh`

Processes a single BAM file for CNV analysis using CNVpytor.

**Usage**:
```bash
./ont_process_individual_bam_cnv.sh \
  -b sample.bam \
  -o ./output_dir \
  -r /path/to/reference.fa \
  -s SAMPLE001 \
  --bin-sizes "1000 10000 100000" \
  -e cnvgen \
  -t 4
```

**Outputs**:
- `SAMPLE001.pytor`: CNVpytor data file
- `SAMPLE001_*bp.vcf`: VCF files for each bin size
- `SAMPLE001_cnv_calls.tsv`: Tab-separated CNV calls
- `SAMPLE001_manhattan_*bp.png`: Manhattan plots
- `SAMPLE001_summary.txt`: Summary statistics
- `SAMPLE001_metadata.json`: Metadata with timestamps

### 2. Time-Series Visualization Script
**File**: `ont_timeseries_cnv_visualization.R`

Creates time-series visualizations from multiple time-point BAM analyses.

**Usage**:
```bash
Rscript ont_timeseries_cnv_visualization.R \
  --input-dir ./output/individual_bams \
  --output-dir ./output/timeseries_plots \
  --sample-id SAMPLE001 \
  --bin-size 100000
```

**Outputs**:
- `SAMPLE001_cnv_count_timeseries.png`: CNV counts over time
- `SAMPLE001_genome_wide_cnv_heatmap.png`: Genome-wide heatmap
- `SAMPLE001_per_chromosome_cnv.png`: Per-chromosome profiles
- `SAMPLE001_cnv_size_distribution.png`: Size distribution plots
- `SAMPLE001_timeseries_summary.txt`: Summary report

### 3. Combined CNV Analysis Script
**File**: `ont_combined_cnv_analysis.R`

Performs comprehensive analysis across all time points to identify CNV patterns.

**Usage**:
```bash
Rscript ont_combined_cnv_analysis.R \
  --input-dir ./output/individual_bams \
  --output-dir ./output/combined_analysis \
  --sample-id SAMPLE001 \
  --overlap-threshold 0.5
```

**Outputs**:
- `SAMPLE001_cnv_classification.png`: CNV classification plot
- `SAMPLE001_persistent_cnv_map.png`: Persistent CNV regions
- `SAMPLE001_cnv_evolution_tracking.png`: CNV evolution heatmap
- `SAMPLE001_cnv_classification.csv`: Classification table
- `SAMPLE001_persistent_cnvs.csv`: Persistent CNV list
- `SAMPLE001_combined_analysis_summary.txt`: Summary report

### 4. Monitoring Script (Main Pipeline)
**File**: `ont_timeseries_cnv_monitor.sh`

Main orchestration script that monitors for new BAM files and runs the complete pipeline.

**Usage**:

**Continuous monitoring mode** (runs until stopped):
```bash
./ont_timeseries_cnv_monitor.sh \
  -i /path/to/bam_directory \
  -o /path/to/output \
  -r /path/to/reference.fa \
  -s SAMPLE001 \
  -c 300 \
  -e cnvgen \
  -t 4
```

**One-time processing mode** (process existing files and exit):
```bash
./ont_timeseries_cnv_monitor.sh \
  -i /path/to/bam_directory \
  -o /path/to/output \
  -r /path/to/reference.fa \
  -s SAMPLE001 \
  --no-continuous
```

**Parameters**:
- `-i, --input-dir`: Directory containing BAM files
- `-o, --output-dir`: Output directory for results
- `-r, --reference`: Reference genome FASTA file
- `-s, --sample-id`: Sample identifier
- `-c, --check-interval`: Check interval in seconds (default: 300)
- `-m, --min-size`: Minimum BAM size in bytes (default: 10000000)
- `-e, --conda-env`: Conda environment name (default: cnvgen)
- `--bin-sizes`: Bin sizes for CNVpytor (default: "1000 10000 100000")
- `-t, --threads`: Number of threads (default: 4)
- `--max-iterations`: Maximum iterations (0 = infinite)
- `--no-continuous`: Process existing files once and exit

## Workflow Example

### Scenario: Real-time CNV monitoring during ONT sequencing

```bash
# Step 1: Activate conda environment
conda activate cnvgen

# Step 2: Start monitoring (runs in continuous mode)
./ont_timeseries_cnv_monitor.sh \
  -i /data/ont_basecalling/bam_output \
  -o /data/cnv_analysis/PATIENT001 \
  -r /reference/hg38.fa \
  -s PATIENT001 \
  -c 600 \
  --bin-sizes "10000 100000" \
  -t 8

# The pipeline will:
# 1. Check for new BAM files every 600 seconds (10 minutes)
# 2. Process each new BAM file with CNVpytor
# 3. Update time-series visualizations after each new BAM
# 4. Update combined analysis showing CNV evolution
# 5. Continue monitoring until stopped (Ctrl+C)
```

### Scenario: Batch processing of existing BAM files

```bash
# Process all existing BAM files once
./ont_timeseries_cnv_monitor.sh \
  -i /data/historical_bams \
  -o /data/cnv_retrospective \
  -r /reference/hg38.fa \
  -s HISTORICAL_SAMPLE \
  --no-continuous \
  -t 16
```

## Output Directory Structure

```
output_directory/
├── individual_bams/
│   ├── SAMPLE001_bam1_20251014_120000/
│   │   ├── SAMPLE001_bam1_20251014_120000.pytor
│   │   ├── SAMPLE001_bam1_20251014_120000_1000bp.vcf
│   │   ├── SAMPLE001_bam1_20251014_120000_10000bp.vcf
│   │   ├── SAMPLE001_bam1_20251014_120000_100000bp.vcf
│   │   ├── SAMPLE001_bam1_20251014_120000_cnv_calls.tsv
│   │   ├── SAMPLE001_bam1_20251014_120000_manhattan_*.png
│   │   ├── SAMPLE001_bam1_20251014_120000_summary.txt
│   │   └── SAMPLE001_bam1_20251014_120000_metadata.json
│   └── SAMPLE001_bam2_20251014_130000/
│       └── ...
├── timeseries_plots/
│   ├── SAMPLE001_cnv_count_timeseries.png
│   ├── SAMPLE001_genome_wide_cnv_heatmap.png
│   ├── SAMPLE001_per_chromosome_cnv.png
│   ├── SAMPLE001_cnv_size_distribution.png
│   └── SAMPLE001_timeseries_summary.txt
├── combined_analysis/
│   ├── SAMPLE001_cnv_classification.png
│   ├── SAMPLE001_persistent_cnv_map.png
│   ├── SAMPLE001_cnv_evolution_tracking.png
│   ├── SAMPLE001_cnv_classification.csv
│   ├── SAMPLE001_persistent_cnvs.csv
│   └── SAMPLE001_combined_analysis_summary.txt
└── logs/
    ├── SAMPLE001_processing_log.txt
    └── SAMPLE001_processed_bams.txt
```

## CNV Classifications

The pipeline classifies CNVs into three categories:

1. **Persistent CNVs**: Present in ALL time points
   - Likely germline or early somatic events
   - Stable genomic alterations

2. **Emerging CNVs**: Present in MULTIPLE (but not all) time points
   - Potential clonal evolution events
   - Progressive genomic changes

3. **Transient CNVs**: Present in ONLY ONE time point
   - Possible technical artifacts
   - Transient biological events
   - Require validation

## Visualization Outputs

### 1. CNV Count Over Time
Line plot showing the number of deletions and duplications at each time point.

### 2. Genome-wide CNV Heatmap
Heatmap showing CNV profiles across all chromosomes over time, with color indicating copy number changes.

### 3. Per-Chromosome CNV Profiles
Faceted plots showing CNV distribution for each chromosome across all time points.

### 4. CNV Size Distribution
Histograms showing the size distribution of CNVs at each time point.

### 5. CNV Classification Plot
Bar plots showing counts of persistent, emerging, and transient CNVs.

### 6. Persistent CNV Map
Genomic map highlighting CNV regions present in all time points.

### 7. CNV Evolution Tracking
Heatmap tracking the top CNV regions across all time points.

## Tips and Best Practices

### For Real-time Monitoring

1. **Set appropriate check interval**: Balance between responsiveness and system load
   - Fast sequencing: `-c 300` (5 minutes)
   - Slow sequencing: `-c 1800` (30 minutes)

2. **Monitor disk space**: CNVpytor files can be large
   ```bash
   du -sh output_directory/
   ```

3. **Use screen/tmux**: Keep monitoring running in background
   ```bash
   screen -S cnv_monitor
   ./ont_timeseries_cnv_monitor.sh [options]
   # Detach: Ctrl+A, D
   # Reattach: screen -r cnv_monitor
   ```

### For Batch Processing

1. **Use more threads**: Process faster with `--threads 16` or higher

2. **Disable continuous mode**: Use `--no-continuous` for one-time processing

3. **Check logs**: Review `logs/SAMPLE_processing_log.txt` for any issues

### Bin Size Selection

- **1kb bins**: High resolution, large files, slower processing
  - Use for: Small regions, targeted analysis

- **10kb bins**: Balanced resolution and speed
  - Use for: Whole genome, most applications

- **100kb bins**: Fast, smaller files, lower resolution
  - Use for: Quick screening, large-scale events

## Troubleshooting

### Issue: "CNVpytor not found"
```bash
# Activate conda environment
conda activate cnvgen

# Verify installation
which cnvpytor
```

### Issue: "Failed to create BAM index"
```bash
# Manually index BAM
samtools index your_file.bam

# Check BAM file integrity
samtools quickcheck your_file.bam
```

### Issue: "R script fails with package error"
```bash
# Activate conda environment and test R packages
conda activate cnvgen
Rscript -e "library(ggplot2)"

# If missing, install within environment
conda install -c conda-forge r-ggplot2
```

### Issue: "No CNV calls generated"
- Check BAM file quality and coverage
- Try different bin sizes
- Review CNVpytor logs in output directory

### Issue: "Time-series plots not generating"
- Need at least 2 time points for time-series analysis
- Check that individual BAM processing completed successfully
- Review R script output in logs

## Performance Considerations

### Memory Usage
- CNVpytor: ~2-4 GB per BAM
- R scripts: ~1-2 GB
- Total: ~8-16 GB recommended

### Processing Time (per BAM)
- Small BAM (<5 GB): ~5-10 minutes
- Medium BAM (5-20 GB): ~15-30 minutes
- Large BAM (>20 GB): ~30-60 minutes

*Times vary based on coverage, bin sizes, and available CPUs*

### Disk Space
- CNVpytor .pytor files: ~10-50% of BAM size
- VCF files: ~1-10 MB each
- Plots: ~1-5 MB each
- Total: ~50-100% of total BAM size

## Citation

If you use this pipeline, please cite:

- **CNVpytor**: Suvakov et al. (2021) CNVpytor: a tool for CNV/CNA detection and analysis from read depth and allele imbalance in whole-genome sequencing. Bioinformatics.

- **Samtools**: Li et al. (2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics.

## Support

For issues, questions, or feature requests, please refer to:
- CNVpytor documentation: https://github.com/abyzovlab/CNVpytor
- Pipeline logs: `output_directory/logs/`

## License

This pipeline is provided as-is for research use.

## Authors

Created for ONT time-series CNV analysis
Date: 2025-10-14
