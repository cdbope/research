# MGMT Methylation Detection using Genome-Wide Methylation Data

A comprehensive machine learning pipeline for detecting MGMT methylation status using genome-wide methylation data from Oxford Nanopore Technologies' modkit package.

## Overview

This modular pipeline analyzes bedmethyl files to predict MGMT (O6-methylguanine-DNA methyltransferase) methylation status, which is crucial for glioblastoma treatment decisions. Unlike traditional approaches that focus only on the MGMT promoter region, this pipeline considers methylation patterns across the entire genome.

## Features

- **Modular Architecture**: Independent functions for data loading, feature engineering, model training, and evaluation
- **Genome-Wide Analysis**: Utilizes methylation patterns from all chromosomes, not just chr10
- **Multiple ML Algorithms**: Random Forest, Gradient Boosting, SVM, Neural Networks, and ensemble methods
- **Comprehensive Feature Engineering**: 
  - Basic methylation statistics
  - Chromosomal patterns
  - Regional analysis with multiple window sizes
  - MGMT-specific features
  - Statistical and distributional features
- **Advanced Evaluation**: ROC curves, precision-recall curves, calibration plots, and feature importance analysis

## Installation

1. Clone or download this repository
2. Install dependencies:

```bash
pip install -r requirements.txt
```

## Quick Start

### 1. Basic Usage (Single Sample)

```bash
# Create default configuration
python main_pipeline.py --input /path/to/sample.wf_mods.bedmethyl.gz --create-config

# Edit the created config.json file, then run:
python main_pipeline.py --config config.json
```

### 2. Multiple Samples

```bash
python main_pipeline.py \
  --input sample1.bedmethyl.gz sample2.bedmethyl.gz \
  --sample-names S1 S2 \
  --output-dir results \
  --metadata metadata.csv
```

### 3. Using the Example Data

```bash
# Run with the provided example
python main_pipeline.py \
  --input /home/chbope/extension/nWGS_manuscript_data/data/results/epi2me/modkit/T001.wf_mods.bedmethyl.gz \
  --output-dir T001_results \
  --min-coverage 3
```

## Input Data Format

The pipeline expects bedmethyl files from Oxford Nanopore's modkit package with the following columns:
- chromosome, start, end, modification_type, score, strand
- thick_start, thick_end, rgb, coverage, percent_modified
- n_modified, n_canonical, n_other_mod, n_delete, n_fail, n_diff

## Sample Metadata

For training models, you need a metadata CSV file with columns:
- `sample_id`: Sample identifier
- `mgmt_status`: MGMT methylation status ("methylated" or "unmethylated")

Example metadata.csv:
```csv
sample_id,mgmt_status,batch,tissue_type
T001,methylated,batch1,tumor
T002,unmethylated,batch1,tumor
```

## Module Description

### 1. `data_loader.py`
- Loads bedmethyl files (gzipped or plain text)
- Handles single or multiple samples
- Provides data validation and quality metrics
- Filters by chromosome, modification type, and coverage

### 2. `feature_engineering.py`
- Extracts comprehensive methylation features:
  - **Basic features**: Global methylation statistics, coverage metrics
  - **Chromosomal features**: Per-chromosome methylation patterns
  - **Regional features**: Window-based analysis, CpG island-like regions
  - **MGMT-specific features**: Promoter region analysis (chr10:129,466,683-129,467,448)
  - **Statistical features**: Distribution properties, spatial correlations

### 3. `model_training.py`
- Supports multiple ML algorithms:
  - Random Forest, Gradient Boosting, Extra Trees
  - Logistic Regression, SVM, Neural Networks
  - Ensemble methods (Voting Classifier)
- Automated hyperparameter tuning
- Feature selection (RFE, univariate, model-based)
- Cross-validation and model persistence

### 4. `model_evaluation.py`
- Comprehensive model evaluation:
  - Performance metrics (accuracy, precision, recall, F1, ROC-AUC)
  - Visualization (ROC curves, precision-recall curves, confusion matrices)
  - Feature importance analysis
  - Model calibration assessment
  - Learning curves

### 5. `main_pipeline.py`
- Orchestrates the complete workflow
- Configuration-based execution
- Command-line interface
- Results summary and reporting

## Configuration

The pipeline uses JSON configuration files. Key parameters:

```json
{
  "input_files": ["path/to/bedmethyl.gz"],
  "metadata_file": "metadata.csv",
  "output_dir": "results",
  "modification_types": ["m"],
  "min_coverage": 3,
  "window_sizes": [1000, 5000, 10000, 50000],
  "feature_selection": true,
  "n_features": 50,
  "hyperparameter_tuning": true,
  "create_ensemble": true
}
```

## Output

The pipeline generates:

```
results/
├── features.csv                    # Extracted features
├── feature_names.txt              # Feature names list
├── data_validation.json           # Data quality metrics
├── sample_metadata_template.csv   # Template for labels
├── models/                        # Trained models
│   ├── random_forest_model.joblib
│   ├── ensemble_model.joblib
│   └── feature_selector.joblib
├── evaluation/                    # Evaluation results
│   ├── roc_curves.png
│   ├── precision_recall_curves.png
│   ├── confusion_matrices.png
│   ├── feature_importance.png
│   └── evaluation_report.txt
├── predictions.csv               # Sample predictions
└── pipeline_summary.json        # Complete results summary
```

## Advanced Usage

### Custom Feature Engineering

```python
from feature_engineering import MethylationFeatureExtractor

# Initialize with custom window sizes
extractor = MethylationFeatureExtractor(window_sizes=[500, 2000, 10000])

# Extract specific feature types
basic_features = extractor.extract_basic_features(data)
mgmt_features = extractor.extract_mgmt_specific_features(data)
```

### Custom Model Training

```python
from model_training import MGMTModelTrainer

trainer = MGMTModelTrainer()
X_train, X_test, y_train, y_test, feature_names = trainer.prepare_data(features_df, metadata_df)

# Train specific models
results = trainer.train_single_model(X_train, y_train, 'random_forest')
```

### Batch Processing

For processing multiple samples in batches:

```python
from data_loader import BedMethylLoader

loader = BedMethylLoader()
file_list = ["sample1.bedmethyl.gz", "sample2.bedmethyl.gz", ...]
samples_data = loader.load_multiple_samples(file_list)
```

## Performance Considerations

- **Memory Usage**: Large bedmethyl files may require 8-16GB RAM
- **Processing Time**: Feature extraction ~5-10 minutes per sample
- **Model Training**: 10-30 minutes depending on hyperparameter tuning
- **Recommendations**: 
  - Use `min_coverage >= 3` to reduce noise
  - Consider chromosome filtering for faster processing
  - Use feature selection to reduce dimensionality

## Citation

If you use this pipeline in your research, please cite:

```
MGMT Methylation Detection Pipeline
Genome-wide methylation analysis for MGMT status prediction
[Your publication details]
```

## License

This project is licensed under the MIT License.

## Support

For issues and questions:
1. Check the generated log files in the output directory
2. Ensure input data format matches bedmethyl specifications
3. Verify sample metadata has correct MGMT status labels

## Examples

### Example 1: Single Sample Analysis

```bash
# Quick analysis of one sample
python main_pipeline.py \
  --input T001.wf_mods.bedmethyl.gz \
  --output-dir T001_analysis \
  --min-coverage 5
```

### Example 2: Multi-Sample Study

```bash
# Analyze multiple samples with metadata
python main_pipeline.py \
  --input *.bedmethyl.gz \
  --metadata study_metadata.csv \
  --output-dir cohort_study \
  --create-config

# Edit config.json for specific parameters, then:
python main_pipeline.py --config cohort_study/config.json
```

### Example 3: Chromosome-Specific Analysis

Create config.json:
```json
{
  "input_files": ["sample.bedmethyl.gz"],
  "chromosomes": ["chr10", "chr1", "chr2"],
  "window_sizes": [1000, 5000],
  "n_features": 30
}
```

## Key Features for MGMT Detection

The pipeline specifically extracts features known to be relevant for MGMT methylation:

1. **MGMT Promoter Region** (chr10:129,466,683-129,467,448)
   - Direct methylation levels in the canonical MGMT promoter
   - Extended promoter region analysis

2. **Chromosome 10 Features**
   - Overall chr10 methylation patterns
   - Relative methylation compared to genome-wide levels

3. **Genome-Wide Context**
   - Global methylation patterns that may correlate with MGMT status
   - Regional methylation signatures
   - Statistical distribution features

4. **Window-Based Analysis**
   - Multi-scale regional analysis (1kb to 50kb windows)
   - CpG island-like region identification
   - Spatial methylation correlations

This comprehensive approach captures both the canonical MGMT promoter methylation and the broader genomic context that may influence MGMT expression and therapeutic response.
