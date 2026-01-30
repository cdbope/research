# Cortical Neuron Probe Selection Methodology

## Overview

This document describes the methodology used to reduce the Moss reference atlas from **7,890 CpG probes** to **500 cortical neuron-specific probes** for improved classification in low-coverage Nanopore data.

## Rationale

The full Moss reference atlas contains ~7,890 CpG probes designed to distinguish 25 different cell types. For neural classification of glioblastoma samples, most of these probes are not informative because they:

1. Show similar methylation levels across all cell types
2. Are specific to non-neural cell types (e.g., hepatocytes, immune cells)
3. Add noise without improving cortical neuron detection

By selecting only probes where **cortical neurons have distinctly different methylation** compared to other cell types, we can:

- Reduce missing data impact in low-coverage samples
- Increase signal-to-noise ratio for neural classification
- Make each observed CpG more informative

## Methodology

### Input Data

- **Reference Atlas**: `reference_atlas.csv`
- **Dimensions**: 7,890 CpG probes × 25 cell types
- **Target Column**: `Cortical_neurons`

### Step 1: Calculate Probe Statistics

For each CpG probe, calculate:

| Metric | Formula | Description |
|--------|---------|-------------|
| `Neuron_beta` | Direct from atlas | Methylation fraction in cortical neurons (0-1) |
| `Other_mean` | Mean of 24 other cell types | Background methylation level |
| `Other_std` | Std of 24 other cell types | Variability across non-neural tissues |
| `Diff_vs_mean` | Neuron_beta - Other_mean | Direction and magnitude of difference |
| `Z_score` | (Neuron_beta - Other_mean) / Other_std | How many standard deviations from background |

### Step 2: Calculate Specificity Score

The **Specificity Score** combines Z-score magnitude with absolute difference:

```
Specificity_score = |Z_score| × |Diff_vs_mean|
```

This metric rewards probes that are:
- Many standard deviations from background (high |Z_score|)
- AND have large absolute difference in methylation (high |Diff|)

### Step 3: Classify Probe Types

Each probe is classified based on its methylation pattern:

| Probe Type | Criteria | Biological Meaning |
|------------|----------|-------------------|
| **Neuron-HIGH** | Neuron_beta > Other_mean + 0.15 AND Neuron_beta > max(Others) | Methylated specifically in neurons |
| **Neuron-LOW** | Neuron_beta < Other_mean - 0.15 AND Neuron_beta < min(Others) | Unmethylated specifically in neurons |
| **Neuron-DIFF** | \|Diff\| > 0.15 (but not extreme) | Different in neurons, but overlaps with some tissues |
| **Non-specific** | \|Diff\| ≤ 0.15 | Similar methylation across all cell types |

### Step 4: Select Top 500 Probes

1. Filter out "Non-specific" probes
2. Sort remaining probes by Specificity Score (descending)
3. Select top 500 probes

## Results

### Probe Type Distribution (Top 500)

Based on the reference atlas analysis:

- **Neuron-HIGH probes**: Probes with high methylation in neurons (e.g., Neuron=0.45, Others=0.01)
- **Neuron-LOW probes**: Probes with low methylation in neurons (e.g., Neuron=0.11, Others=0.96)

### Example Top Probes

| CpG | Neuron Beta | Others Mean | Z-score | Type |
|-----|-------------|-------------|---------|------|
| cg14149007 | 0.457 | 0.006 | +131.6 | Neuron-HIGH |
| cg07209034 | 0.112 | 0.959 | -67.6 | Neuron-LOW |
| cg18477204 | 0.774 | 0.025 | +71.0 | Neuron-HIGH |
| cg07804711 | 0.735 | 0.041 | +66.5 | Neuron-HIGH |
| cg02300356 | 0.521 | 0.977 | -99.6 | Neuron-LOW |

### Interpretation

- **cg14149007**: 46% methylated in neurons vs 0.6% in other tissues. Observing methylation at this site strongly indicates neural origin.
- **cg07209034**: 11% methylated in neurons vs 96% in other tissues. Observing unmethylation at this site strongly indicates neural origin.

## Output Files

| File | Description |
|------|-------------|
| `cortical_neuron_reference_atlas.csv` | Reduced atlas (500 probes) for use with deconvolve.py |
| `cortical_neuron_specific_probes.csv` | Full analysis with all metrics for each probe |
| `cortical_neuron_cpg_list.txt` | Simple list of 500 CpG IDs |

## Usage

### For Deconvolution

```bash
python deconvolve.py -a cortical_neuron_reference_atlas.csv input_betas.csv
```

### For Custom Probe Extraction

Modify parameters in `extract_cortical_neuron_probes.py`:

```python
MIN_DIFF = 0.15        # Minimum absolute difference threshold
TOP_N_PROBES = 500     # Number of probes to select
```

## Advantages for Low-Coverage Data

1. **Higher hit rate**: 500 targeted probes vs 7,890 means each read covering a probe is more likely to be informative
2. **Stronger priors**: Selected probes have extreme expected values (near 0 or near 1), making even single-read observations meaningful
3. **Reduced noise**: Non-informative probes don't dilute the signal
4. **Bayesian compatibility**: Prior distributions from training data are more distinct, improving classification with sparse observations

## References

- Moss et al. - Original reference atlas methodology
- Drexler et al. (Nature Medicine 2024) - Neural glioblastoma classification
- Script: `extract_cortical_neuron_probes.py`
