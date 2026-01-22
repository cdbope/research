# Installation Guide for Drexler et al. Epigenetic Neural Glioblastoma Code

This guide provides step-by-step instructions for setting up the environment to run the code from:

**"A prognostic neural epigenetic signature in high-grade glioma"** (Drexler et al., Nature Medicine 2024)

## Prerequisites

- Linux (tested on Ubuntu)
- Miniconda or Anaconda installed

## Step 1: Create Conda Environment with Python and R

```bash
# Create new conda environment with Python 3.10 and R 4.4
conda create -n drexler_env python=3.10 -y
conda activate drexler_env

# Upgrade R to version 4.4 (required for Bioconductor 3.20)
conda install -c conda-forge r-base=4.4 -y
```

## Step 2: Install Python Dependencies

```bash
conda activate drexler_env
pip install numpy==1.26.3 pandas==2.1.1 scikit-learn==1.2.2 joblib==1.3.2 scipy matplotlib seaborn lifelines jupyter
```

## Step 3: Install R Dependencies via Conda

Due to compiler compatibility issues between conda's gcc and system headers, it's recommended to install R packages via conda rather than from CRAN/Bioconductor directly.

```bash
conda activate drexler_env

# Install r-curl (critical - fails to compile from source)
conda install -c conda-forge r-curl -y

# Install core Bioconductor packages
conda install -c conda-forge -c bioconda r-httr bioconductor-genomeinfodb bioconductor-genomicranges bioconductor-biostrings -y

# Install GenomicFeatures and bumphunter
conda install -c conda-forge -c bioconda bioconductor-genomicfeatures bioconductor-bumphunter -y
```

## Step 4: Install Remaining Bioconductor Packages

After installing the conda packages, install the remaining packages via BiocManager in R:

```bash
conda activate drexler_env
Rscript -e '
# Install BiocManager if not present
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")

# Use Bioconductor 3.20 for R 4.4
BiocManager::install(version = "3.20", ask=FALSE)

# Install minfi and Illumina manifests
BiocManager::install(c(
    "minfi",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
), ask=FALSE, update=FALSE)
'
```

## Step 5: Verify Installation

### Test R packages:
```bash
conda activate drexler_env
Rscript -e '
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
cat("All R packages loaded successfully!\n")
'
```

### Test Python packages:
```bash
conda activate drexler_env
python -c '
import numpy as np
import pandas as pd
import sklearn
import joblib
print(f"numpy: {np.__version__}")
print(f"pandas: {pd.__version__}")
print(f"scikit-learn: {sklearn.__version__}")
print(f"joblib: {joblib.__version__}")
print("All Python packages loaded successfully!")
'
```

## Running the Example

Navigate to the example directory and run the notebook:

```bash
conda activate drexler_env
cd code/example_neural_classification
jupyter notebook run_example.ipynb
```

Or run it as a script:

```bash
conda activate drexler_env
cd code/example_neural_classification

# Step 1: Process IDAT files
Rscript ../DNAm_deconv/process_array.R sample_idats betas.csv ../DNAm_deconv/ref_sample.RData

# Step 2: Run prediction
python -c "
import pandas as pd
import numpy as np
import joblib

# Load beta values
df = pd.read_csv('betas.csv', index_col=0).T

# Load CpG sites
low_cgs = pd.read_csv('../neural_group_classification/low_manifest_all.csv', index_col=0).index.tolist()
high_cgs = pd.read_csv('../neural_group_classification/high_manifest_all.csv', index_col=0).index.tolist()
cgs = low_cgs + high_cgs

# Subset to required CpGs
df = df[cgs]

# Load model and predict
clf = joblib.load('../neural_group_classification/logregCV_allCpG.pkl')
preds = clf.predict(np.array(df))
preds_score = clf.predict_proba(np.array(df)).max(1)

# Save results
df_preds = pd.DataFrame([preds, preds_score], index=['Prediction', 'Prediction score'], columns=df.index).T
df_preds['Prediction'] = df_preds['Prediction'].astype(int)
df_preds.to_csv('prediction.csv')
print(df_preds)
"
```

## Troubleshooting

### Issue: R curl package fails to compile
**Error:** `error: unknown type name '__time64_t'`

**Solution:** Install r-curl via conda instead of from CRAN:
```bash
conda install -c conda-forge r-curl -y
```

### Issue: Bioconductor version mismatch
**Error:** `Bioconductor version '3.18' requires R version '4.3'`

**Solution:** Specify the correct Bioconductor version for your R version:
- R 4.3 → Bioconductor 3.18
- R 4.4 → Bioconductor 3.20

```R
BiocManager::install(version = "3.20", ask=FALSE)
```

### Issue: Matrix/MASS packages not available
**Error:** `package 'Matrix' is not available for this version of R`

**Solution:** Upgrade R to version 4.4:
```bash
conda install -c conda-forge r-base=4.4 -y
```

## Package Versions (Tested)

### Python
- numpy: 1.26.3
- pandas: 2.1.1
- scikit-learn: 1.2.2
- joblib: 1.3.2

### R/Bioconductor
- R: 4.4.x
- Bioconductor: 3.20
- minfi: latest
- IlluminaHumanMethylation450kmanifest: latest
- IlluminaHumanMethylationEPICmanifest: latest
- IlluminaHumanMethylation450kanno.ilmn12.hg19: latest

## Citation

If you use this code, please cite:

Drexler, R., Khatri, R., Sauvigny, T. et al. A prognostic neural epigenetic signature in high-grade glioma. Nat Med (2024). https://doi.org/10.1038/s41591-024-02969-w
