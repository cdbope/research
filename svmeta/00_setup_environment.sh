#!/bin/bash
###############################################################################
# Setup Script for SV Meta-Analysis Pipeline
#
# This script:
# 1. Creates the svmeta_env conda environment
# 2. Installs SURVIVOR
# 3. Verifies all dependencies
#
# Usage: ./00_setup_environment.sh
###############################################################################

set -e

echo "============================================================================"
echo "SV META-ANALYSIS PIPELINE - ENVIRONMENT SETUP"
echo "============================================================================"
echo

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

# ============================================================================
# Step 1: Create conda environment
# ============================================================================

echo "Step 1: Creating svmeta_env conda environment..."
echo "------------------------------------------------------------"

if conda env list | grep -q "^svmeta_env "; then
    echo "⚠ svmeta_env already exists"
    read -p "Do you want to remove and recreate it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing svmeta_env..."
        conda env remove -n svmeta_env -y
    else
        echo "Skipping conda environment creation"
        SKIP_CONDA=true
    fi
fi

if [ "$SKIP_CONDA" != "true" ]; then
    echo "Creating conda environment from environment.yml..."
    conda env create -f environment.yml
    echo "✓ Conda environment created"
fi

echo

# ============================================================================
# Step 2: Activate environment and verify
# ============================================================================

echo "Step 2: Activating environment and verifying dependencies..."
echo "------------------------------------------------------------"

# Source conda
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate environment
conda activate svmeta_env

# Check Python packages
echo "Checking Python packages..."
python -c "import pandas; import numpy; import matplotlib; import seaborn; import scipy; import sklearn" && \
    echo "✓ Core Python packages installed" || \
    echo "✗ Error: Some Python packages missing"

python -c "import umap" && \
    echo "✓ umap-learn installed" || \
    echo "✗ Error: umap-learn missing"

# Check bcftools
if command -v bcftools &> /dev/null; then
    BCFTOOLS_VERSION=$(bcftools --version | head -n1)
    echo "✓ bcftools installed: $BCFTOOLS_VERSION"
else
    echo "✗ Error: bcftools not found"
fi

# Check samtools
if command -v samtools &> /dev/null; then
    SAMTOOLS_VERSION=$(samtools --version | head -n1)
    echo "✓ samtools installed: $SAMTOOLS_VERSION"
else
    echo "✗ Error: samtools not found"
fi

echo

# ============================================================================
# Step 3: Install SURVIVOR
# ============================================================================

echo "Step 3: Installing SURVIVOR..."
echo "------------------------------------------------------------"

SURVIVOR_DIR="$SCRIPT_DIR/tools/SURVIVOR"

if [ -f "$SURVIVOR_DIR/Debug/SURVIVOR" ]; then
    echo "⚠ SURVIVOR already installed at $SURVIVOR_DIR"
    read -p "Do you want to reinstall it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing SURVIVOR..."
        rm -rf "$SURVIVOR_DIR"
    else
        echo "Skipping SURVIVOR installation"
        SKIP_SURVIVOR=true
    fi
fi

if [ "$SKIP_SURVIVOR" != "true" ]; then
    # Create tools directory
    mkdir -p "$SCRIPT_DIR/tools"
    cd "$SCRIPT_DIR/tools"

    # Clone SURVIVOR
    echo "Cloning SURVIVOR from GitHub..."
    git clone https://github.com/fritzsedlazeck/SURVIVOR.git
    cd SURVIVOR/Debug

    # Build SURVIVOR
    echo "Building SURVIVOR..."
    make

    if [ -f "SURVIVOR" ]; then
        echo "✓ SURVIVOR built successfully"
        cd "$SCRIPT_DIR"
    else
        echo "✗ Error: SURVIVOR build failed"
        cd "$SCRIPT_DIR"
        exit 1
    fi
fi

# Verify SURVIVOR
if [ -f "$SURVIVOR_DIR/Debug/SURVIVOR" ]; then
    SURVIVOR_VERSION=$("$SURVIVOR_DIR/Debug/SURVIVOR" 2>&1 | head -n1 || echo "Unknown")
    echo "✓ SURVIVOR installed: $SURVIVOR_VERSION"
else
    echo "✗ Error: SURVIVOR not found at $SURVIVOR_DIR/Debug/SURVIVOR"
fi

echo

# ============================================================================
# Step 4: Update script paths
# ============================================================================

echo "Step 4: Updating SURVIVOR path in scripts..."
echo "------------------------------------------------------------"

SURVIVOR_PATH="$SURVIVOR_DIR/Debug/SURVIVOR"

# Update 02_merge_with_survivor.sh
if [ -f "02_merge_with_survivor.sh" ]; then
    # Backup original
    cp 02_merge_with_survivor.sh 02_merge_with_survivor.sh.bak

    # Update SURVIVOR_BIN path
    sed -i "s|SURVIVOR_BIN=\"SURVIVOR\"|SURVIVOR_BIN=\"$SURVIVOR_PATH\"|g" 02_merge_with_survivor.sh

    echo "✓ Updated 02_merge_with_survivor.sh with SURVIVOR path"
fi

echo

# ============================================================================
# Step 5: Summary and next steps
# ============================================================================

echo "============================================================================"
echo "SETUP COMPLETE!"
echo "============================================================================"
echo
echo "Environment: svmeta_env"
echo "SURVIVOR location: $SURVIVOR_PATH"
echo
echo "To activate the environment:"
echo "  conda activate svmeta_env"
echo
echo "To deactivate:"
echo "  conda deactivate"
echo
echo "Next steps:"
echo "  1. Activate environment: conda activate svmeta_env"
echo "  2. Configure paths in scripts:"
echo "     - Edit 01_normalize_vcfs.sh (REF_FASTA, INPUT_VCF_DIR)"
echo "     - Edit 03_build_matrix_and_analyze.py (GENE_ANNOTATION)"
echo "  3. Run pipeline:"
echo "     ./01_normalize_vcfs.sh"
echo "     ./02_merge_with_survivor.sh"
echo "     python 03_build_matrix_and_analyze.py"
echo
echo "For more information, see README.md or QUICKSTART.md"
echo "============================================================================"

# Create activation helper script
cat > activate_svmeta.sh << 'EOF'
#!/bin/bash
# Helper script to activate svmeta_env

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate svmeta_env

echo "svmeta_env activated!"
echo "SURVIVOR location: $(dirname $(dirname $(readlink -f "$0")))/tools/SURVIVOR/Debug/SURVIVOR"
EOF

chmod +x activate_svmeta.sh
echo "✓ Created activate_svmeta.sh helper script"
echo

conda deactivate
