#!/bin/bash
###############################################################################
# Install g:Profiler for Pathway Enrichment Analysis
#
# This installs the gprofiler-official Python package for automated
# pathway enrichment analysis.
#
# Usage: ./install_gprofiler.sh
###############################################################################

echo "============================================================================"
echo "Installing g:Profiler for Pathway Enrichment"
echo "============================================================================"
echo

# Check if conda environment is active
if [[ "$CONDA_DEFAULT_ENV" != "svmeta_env" ]]; then
    echo "⚠ svmeta_env is not active"
    echo "Please activate the environment first:"
    echo "  conda activate svmeta_env"
    echo
    exit 1
fi

echo "Installing gprofiler-official..."
pip install gprofiler-official

echo
echo "============================================================================"
echo "Testing installation..."
echo "============================================================================"

python -c "from gprofiler import GProfiler; print('✓ g:Profiler installed successfully!')"

echo
echo "============================================================================"
echo "Installation complete!"
echo "============================================================================"
echo
echo "You can now run pathway enrichment:"
echo "  python 05_pathway_enrichment.py"
echo
