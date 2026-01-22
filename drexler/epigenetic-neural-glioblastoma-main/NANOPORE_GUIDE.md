# Nanopore Methylation to Drexler Neural Glioblastoma Classifier

## Overview

This guide explains how to use Oxford Nanopore methylation data (bedMethyl format from modkit) with the Drexler et al. neural glioblastoma classifier, which was originally designed for Illumina 450K/EPIC methylation arrays.

**Reference Paper:** Drexler, R., et al. "A prognostic neural epigenetic signature in high-grade glioma." Nature Medicine (2024).

## Background

### The Challenge

The Drexler classifier uses **1,289 specific CpG probes** from Illumina methylation arrays to classify glioblastoma samples into:
- **Low-Neural (Class 0):** Associated with worse prognosis
- **High-Neural (Class 1):** Associated with better prognosis

Nanopore sequencing provides methylation calls at genomic coordinates, while Illumina arrays use probe IDs (e.g., `cg20367788`). The conversion scripts bridge this gap by:

1. Mapping Illumina probe IDs to genomic coordinates (hg19)
2. Extracting methylation values from your bedMethyl file at those positions
3. Converting to beta values (0-1 scale) compatible with the classifier

### Data Flow

```
Nanopore BAM → modkit → bedMethyl → nanopore_to_drexler.py → betas.csv → classifier → prediction
```

## Prerequisites

### 1. Nanopore Data Requirements

Your bedMethyl file should be generated using **modkit** from Oxford Nanopore:

```bash
# Example modkit command to generate bedMethyl
modkit pileup input.bam output.bedmethyl.gz --ref reference.fasta --preset traditional
```

**Important:**
- The reference genome should be **hg19/GRCh37** (same as Illumina 450K annotation)
- If using hg38, you'll need to liftover coordinates first

### 2. bedMethyl File Format

The script expects standard modkit bedMethyl format:

| Column | Description | Example |
|--------|-------------|---------|
| 1 | Chromosome | chr1 |
| 2 | Start (0-based) | 10286 |
| 3 | End | 10287 |
| 4 | Modification code | m (5mC) or h (5hmC) |
| 5 | Score | 1 |
| 6 | Strand | + or - |
| 11 | Percent methylated | 100.00 |
| 12 | N_modified | 5 |
| 13 | N_canonical | 10 |

**Note:** Only `m` (5-methylcytosine) calls are used; `h` (5-hydroxymethylcytosine) is ignored.

### 3. Environment Setup

Ensure you have the drexler_env conda environment set up (see INSTALL.md):

```bash
conda activate drexler_env
```

## Usage

### Quick Start

```bash
cd /home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code

# Run the complete pipeline
python run_nanopore_classification.py /path/to/your/sample.bedmethyl.gz --output_dir ../results/
```

### Command-Line Options

```
python run_nanopore_classification.py <bedmethyl_file> [options]

Required:
  bedmethyl           Path to bedMethyl file (gzipped or plain text)

Options:
  --output_dir DIR    Output directory (default: current directory)
  --sample_name NAME  Sample name (default: derived from filename)
  --min_coverage N    Minimum read coverage per CpG site (default: 5)
  --genome {hg19,hg38}  Reference genome build (default: hg38)
  --probe_coords CSV  Pre-computed probe coordinates file (optional)
```

### Genome Build Selection

The `--genome` flag is **critical** for correct probe matching:

| Your Reference | Use |
|----------------|-----|
| GRCh38 / hg38 | `--genome hg38` (default) |
| GRCh37 / hg19 | `--genome hg19` |

On first run with a genome build, the pipeline will:
1. Load Illumina 450K probe coordinates (hg19)
2. If `--genome hg38`: perform liftOver to convert coordinates
3. Cache the result for future runs

### Example Commands

```bash
# Basic usage (hg38 reference - default)
python run_nanopore_classification.py ../T001.wf_mods.bedmethyl.gz

# Explicitly specify hg38
python run_nanopore_classification.py ../T001.wf_mods.bedmethyl.gz --genome hg38

# If your data uses hg19 reference
python run_nanopore_classification.py ../T001.wf_mods.bedmethyl.gz --genome hg19

# Specify output directory and sample name
python run_nanopore_classification.py ../T001.wf_mods.bedmethyl.gz \
    --genome hg38 \
    --output_dir ../results/ \
    --sample_name Patient_001

# Require higher coverage (more stringent)
python run_nanopore_classification.py ../T001.wf_mods.bedmethyl.gz \
    --genome hg38 \
    --min_coverage 10 \
    --output_dir ../results/
```

## Output Files

The pipeline generates two output files:

### 1. Beta Values File (`{sample}_betas.csv`)

Contains the methylation beta values (0-1 scale) for each CpG probe:

```csv
,Patient_001
cg20367788,0.713
cg01847620,0.970
cg08450501,0.026
...
```

### 2. Prediction File (`{sample}_prediction.csv`)

Contains the classification results:

```csv
,Prediction,Prediction_Score,Prob_Class_0,Prob_Class_1,Neural_Classification
Patient_001,0,0.992,0.992,0.008,Low-Neural
```

| Column | Description |
|--------|-------------|
| Prediction | 0 (Low-Neural) or 1 (High-Neural) |
| Prediction_Score | Confidence score (0-1) |
| Prob_Class_0 | Probability of Low-Neural |
| Prob_Class_1 | Probability of High-Neural |
| Neural_Classification | Human-readable classification |

## Step-by-Step Pipeline Details

### Step 1: Generate Probe Coordinates

The first run will generate `probe_coordinates_450k.csv` from the Illumina annotation:

```
Generating probe coordinates from Illumina450k annotation...
Saved probe coordinates
Loaded 485,512 probe coordinates
Classifier requires 1,289 probes
Found coordinates for 1,289 required probes
```

This file is cached for future runs.

### Step 2: Parse bedMethyl File

The script extracts 5mC methylation calls:

```
Parsing bedMethyl file: T001.wf_mods.bedmethyl.gz
  Processed 5,000,000 lines...
  Processed 10,000,000 lines...
  Loaded 28,456,789 5mC methylation calls
```

### Step 3: Map to Illumina Probes

Methylation values are matched to probe positions:

```
Mapping probes to Nanopore methylation values...
  Matched: 1,156 probes
  Low coverage (<5x): 45 probes
  Not found: 88 probes

Probe coverage: 1156/1289 (89.7%)
```

### Step 4: Impute Missing Values

Missing probes are filled using mean values from the training data:

```
Imputing missing values...
  Columns with missing values: 133
```

### Step 5: Classification

The logistic regression classifier predicts the neural signature:

```
CLASSIFICATION RESULTS
============================================================

Sample: Patient_001
  Classification: Low-Neural
  Confidence: 0.992
  Probabilities:
    Low-Neural (Class 0): 0.992
    High-Neural (Class 1): 0.008
```

## Interpreting Results

### Classification Meaning

| Classification | Description | Prognosis |
|----------------|-------------|-----------|
| **Low-Neural (0)** | Low expression of neural differentiation markers | Worse survival |
| **High-Neural (1)** | High expression of neural differentiation markers | Better survival |

### Confidence Score

- **> 0.9:** High confidence prediction
- **0.7 - 0.9:** Moderate confidence
- **< 0.7:** Low confidence, interpret with caution

### Probe Coverage Considerations

| Coverage | Interpretation |
|----------|----------------|
| **> 90%** | Excellent - reliable prediction |
| **80-90%** | Good - prediction should be reliable |
| **70-80%** | Moderate - some uncertainty |
| **< 70%** | Low - interpret with caution |

Low coverage may result from:
- Insufficient sequencing depth
- Different genome build (hg38 vs hg19)
- Regions not covered by sequencing

## Troubleshooting

### Issue: Low probe coverage

**Symptoms:** Less than 80% of probes have values

**Solutions:**
1. Increase sequencing depth
2. Check genome build (should be hg19)
3. Lower `--min_coverage` threshold (trade-off: less reliable individual values)

```bash
# Try lower coverage threshold
python run_nanopore_classification.py sample.bedmethyl.gz --min_coverage 3
```

### Issue: Wrong genome build

**Symptoms:** Very low probe matching (<50%)

**Solution:** Convert coordinates from hg38 to hg19 using liftOver:

```bash
# Download liftOver chain file
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz

# Convert bedMethyl (need to extract relevant columns first)
liftOver input_hg38.bed hg38ToHg19.over.chain.gz output_hg19.bed unmapped.bed
```

### Issue: Memory errors with large files

**Symptoms:** Out of memory when parsing large bedMethyl files

**Solution:** Filter the bedMethyl file to only include 5mC calls before processing:

```bash
# Pre-filter to only 5mC calls
zcat large_file.bedmethyl.gz | awk '$4=="m"' | gzip > filtered_5mC.bedmethyl.gz
```

### Issue: No 5mC calls found

**Symptoms:** "Loaded 0 5mC methylation calls"

**Check:**
1. Verify modkit was run with 5mC calling enabled
2. Check the 4th column of your bedMethyl file contains 'm'

```bash
# Check what modification types are present
zcat sample.bedmethyl.gz | cut -f4 | sort | uniq -c
```

## Advanced Usage

### Batch Processing Multiple Samples

```bash
#!/bin/bash
# batch_classify.sh

OUTPUT_DIR="../results"
mkdir -p $OUTPUT_DIR

for bedmethyl in /path/to/samples/*.bedmethyl.gz; do
    sample=$(basename $bedmethyl .wf_mods.bedmethyl.gz)
    echo "Processing: $sample"

    python run_nanopore_classification.py $bedmethyl \
        --output_dir $OUTPUT_DIR \
        --sample_name $sample \
        --min_coverage 5
done

# Combine all predictions
echo "Sample,Prediction,Score,Classification" > $OUTPUT_DIR/all_predictions.csv
for pred in $OUTPUT_DIR/*_prediction.csv; do
    tail -n +2 $pred >> $OUTPUT_DIR/all_predictions.csv
done
```

### Using Pre-computed Probe Coordinates

For faster repeated runs, save and reuse probe coordinates:

```bash
# First run generates probe_coordinates_450k.csv
python run_nanopore_classification.py sample1.bedmethyl.gz

# Subsequent runs can use the cached file
python run_nanopore_classification.py sample2.bedmethyl.gz \
    --probe_coords probe_coordinates_450k.csv
```

### Custom Analysis with Beta Values

After running the pipeline, you can use the beta values for custom analysis:

```python
import pandas as pd
import numpy as np

# Load beta values
betas = pd.read_csv('results/Patient_001_betas.csv', index_col=0)

# Check distribution
print(f"Mean beta: {betas.values.mean():.3f}")
print(f"Missing values: {betas.isna().sum().sum()}")

# Compare specific probes
signature_probes = ['cg20367788', 'cg01847620', 'cg08450501']
print(betas.loc[signature_probes])
```

## Technical Notes

### Coordinate System

- **Illumina 450K annotation:** 1-based coordinates (hg19)
- **bedMethyl format:** 0-based coordinates
- The script automatically converts between coordinate systems

### Strand Handling

CpG sites are palindromic, so methylation can be measured on either strand. The script:
1. First checks the strand specified in the Illumina annotation
2. If not found, checks the opposite strand
3. Uses the first available measurement

### Beta Value Calculation

Beta values are calculated as:

```
beta = percent_methylated / 100
     = N_modified / (N_modified + N_canonical)
```

Values range from 0 (unmethylated) to 1 (fully methylated).

## Citation

If you use this pipeline, please cite:

1. **Drexler et al.** A prognostic neural epigenetic signature in high-grade glioma. *Nature Medicine* (2024). https://doi.org/10.1038/s41591-024-02969-w

2. **modkit** (for methylation calling): https://github.com/nanoporetech/modkit

## Support

For issues with:
- **This conversion pipeline:** Open an issue describing your problem
- **The original Drexler classifier:** See the original paper's supplementary materials
- **modkit/Nanopore data:** Consult Oxford Nanopore documentation
