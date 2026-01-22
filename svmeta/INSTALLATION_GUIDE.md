# Installation Guide - svmeta_env Conda Environment

## Overview

This guide covers the installation and setup of the `svmeta_env` conda environment for the SV meta-analysis pipeline.

## Quick Installation (Automated)

The easiest way to set up everything:

```bash
cd /home/chbope/extension/script/svmeta
./00_setup_environment.sh
```

This script will:
1. ✅ Create `svmeta_env` conda environment
2. ✅ Install all Python packages (pandas, numpy, matplotlib, seaborn, scipy, scikit-learn, umap-learn)
3. ✅ Install bcftools and samtools
4. ✅ Download and build SURVIVOR
5. ✅ Update script paths automatically
6. ✅ Create activation helper script

**Time**: ~10-15 minutes (mostly downloading and compiling SURVIVOR)

## What Gets Installed

### Conda Environment: svmeta_env

**Python Version**: 3.9

**Bioinformatics Tools**:
- bcftools ≥ 1.10
- samtools ≥ 1.10
- htslib ≥ 1.10

**Python Packages**:
- pandas ≥ 1.3.0
- numpy ≥ 1.21.0
- matplotlib ≥ 3.4.0
- seaborn ≥ 0.11.0
- scipy ≥ 1.7.0
- scikit-learn ≥ 0.24.0
- umap-learn ≥ 0.5.0

**External Tool**:
- SURVIVOR (installed in `tools/SURVIVOR/`)

## Step-by-Step Manual Installation

If you prefer to install components manually:

### Step 1: Create Conda Environment

```bash
cd /home/chbope/extension/script/svmeta
conda env create -f environment.yml
```

This creates the `svmeta_env` environment with all dependencies.

### Step 2: Activate Environment

```bash
conda activate svmeta_env
```

### Step 3: Verify Installation

```bash
# Check Python packages
python -c "import pandas, numpy, matplotlib, seaborn, scipy, sklearn, umap; print('✓ All packages OK')"

# Check bcftools
bcftools --version

# Check samtools
samtools --version
```

### Step 4: Install SURVIVOR

```bash
# Create tools directory
mkdir -p tools
cd tools

# Clone and build SURVIVOR
git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug
make

# Test SURVIVOR
./SURVIVOR

cd ../../../
```

### Step 5: Update Script Paths

Edit `02_merge_with_survivor.sh` to point to your SURVIVOR installation:

```bash
SURVIVOR_BIN="/home/chbope/extension/script/svmeta/tools/SURVIVOR/Debug/SURVIVOR"
```

## Using the Environment

### Activate Environment

Every time you want to run the pipeline:

```bash
conda activate svmeta_env
```

**Tip**: Use the helper script:
```bash
source activate_svmeta.sh
```

### Deactivate Environment

When done:

```bash
conda deactivate
```

### Check What's Installed

```bash
conda list
```

## Directory Structure After Installation

```
svmeta/
├── 00_setup_environment.sh        # Setup script
├── 01_normalize_vcfs.sh
├── 02_merge_with_survivor.sh
├── 03_build_matrix_and_analyze.py
├── sv_meta_pipeline.py
├── environment.yml                # Conda environment definition
├── requirements.txt               # Alternative pip requirements
├── activate_svmeta.sh            # Helper activation script (created by setup)
├── gbm_driver_genes.txt
├── sample_metadata_template.txt
├── README.md
├── QUICKSTART.md
├── PROJECT_SUMMARY.md
├── INSTALLATION_GUIDE.md         # This file
└── tools/                        # Created by setup script
    └── SURVIVOR/
        └── Debug/
            └── SURVIVOR          # Compiled SURVIVOR binary
```

## Troubleshooting

### Conda environment already exists

If you get an error that `svmeta_env` already exists:

**Option 1**: Remove and recreate
```bash
conda env remove -n svmeta_env
./00_setup_environment.sh
```

**Option 2**: Update existing environment
```bash
conda activate svmeta_env
conda env update -f environment.yml
```

### SURVIVOR build fails

**Error**: `make: g++: Command not found`

**Solution**: Install build tools
```bash
# On Ubuntu/Debian
sudo apt-get install build-essential

# On CentOS/RHEL
sudo yum groupinstall "Development Tools"

# Then retry
cd tools/SURVIVOR/Debug
make
```

### bcftools not found after installation

**Solution**: Make sure environment is activated
```bash
conda activate svmeta_env
which bcftools  # Should show path in conda env
```

### Python package import errors

**Solution**: Reinstall packages
```bash
conda activate svmeta_env
pip install --force-reinstall umap-learn
```

## Verifying Your Installation

Run this verification script:

```bash
conda activate svmeta_env

# Test all components
python << 'EOF'
import sys
print("Testing Python packages...")
try:
    import pandas
    import numpy
    import matplotlib
    import seaborn
    import scipy
    import sklearn
    import umap
    print("✓ All Python packages OK")
except ImportError as e:
    print(f"✗ Error: {e}")
    sys.exit(1)
EOF

# Test bcftools
if command -v bcftools &> /dev/null; then
    echo "✓ bcftools found: $(bcftools --version | head -n1)"
else
    echo "✗ bcftools not found"
fi

# Test samtools
if command -v samtools &> /dev/null; then
    echo "✓ samtools found: $(samtools --version | head -n1)"
else
    echo "✗ samtools not found"
fi

# Test SURVIVOR
SURVIVOR_PATH="tools/SURVIVOR/Debug/SURVIVOR"
if [ -f "$SURVIVOR_PATH" ]; then
    echo "✓ SURVIVOR found at $SURVIVOR_PATH"
else
    echo "✗ SURVIVOR not found at $SURVIVOR_PATH"
fi

echo
echo "Installation verification complete!"
```

## Updating the Environment

If you need to add new packages:

### Option 1: Update environment.yml

1. Edit `environment.yml` to add new package
2. Update environment:
```bash
conda activate svmeta_env
conda env update -f environment.yml --prune
```

### Option 2: Install directly

```bash
conda activate svmeta_env
conda install new-package-name

# or with pip
pip install new-package-name
```

## Removing the Environment

If you need to completely remove the environment:

```bash
# Remove conda environment
conda env remove -n svmeta_env

# Remove SURVIVOR
rm -rf tools/SURVIVOR

# Remove helper script
rm -f activate_svmeta.sh
```

## Environment Export

To share your exact environment with others:

```bash
conda activate svmeta_env
conda env export > environment_exact.yml
```

Others can then create the exact same environment:

```bash
conda env create -f environment_exact.yml
```

## Alternative: Using Mamba

For faster installation, you can use `mamba` instead of `conda`:

```bash
# Install mamba (one-time)
conda install -n base -c conda-forge mamba

# Create environment with mamba
mamba env create -f environment.yml

# Much faster than conda!
```

## System Requirements

### Minimum
- **OS**: Linux (Ubuntu 18.04+, CentOS 7+, or similar)
- **RAM**: 16 GB
- **Disk**: 5 GB for software, 50+ GB for analysis
- **CPU**: 4+ cores

### Recommended
- **OS**: Linux (Ubuntu 20.04+)
- **RAM**: 32 GB
- **Disk**: 100+ GB SSD
- **CPU**: 8+ cores

## Support

### Installation Issues

1. Check conda is installed and updated:
```bash
conda --version
conda update -n base -c defaults conda
```

2. Try creating environment with verbose output:
```bash
conda env create -f environment.yml -v
```

3. If specific package fails, try installing without it first, then add it separately

### SURVIVOR Issues

- GitHub: https://github.com/fritzsedlazeck/SURVIVOR/issues
- Check build requirements: g++, make, zlib

### Python Package Issues

- Try using pip instead of conda for problematic packages
- Check Python version: `python --version` (should be 3.9)

## Next Steps

After installation is complete:

1. ✅ Verify installation (see verification script above)
2. ✅ Configure pipeline paths:
   - Edit `01_normalize_vcfs.sh` (REF_FASTA, INPUT_VCF_DIR)
   - Edit `03_build_matrix_and_analyze.py` (GENE_ANNOTATION)
3. ✅ Run test with small dataset
4. ✅ Run full analysis

See [QUICKSTART.md](QUICKSTART.md) for running the pipeline.

---

**Questions?** See [README.md](README.md) for full documentation.
