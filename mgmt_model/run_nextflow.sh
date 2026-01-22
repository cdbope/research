#!/bin/bash

# MGMT Methylation Detection - Nextflow Pipeline Runner
# 
# This script demonstrates different ways to run the Nextflow pipeline

echo "=========================================="
echo "MGMT Methylation Detection Pipeline"
echo "=========================================="

# Activate conda environment
echo "Activating mgmt_env environment..."
source /home/chbope/miniconda3/bin/activate mgmt_env

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
fi

echo "Nextflow version:"
nextflow -version

echo ""
echo "Available execution modes:"
echo "1. Single sample analysis"
echo "2. Multi-sample analysis with CSV input"
echo "3. Test mode (fast execution)"
echo "4. Show help"
echo ""

# Default: Single sample analysis
echo "Running single sample analysis..."
echo "Command: nextflow run main.nf --input /home/chbope/extension/nWGS_manuscript_data/data/results/epi2me/modkit/T001.wf_mods.bedmethyl.gz --metadata metadata.csv --outdir results_single"

# Uncomment the line below to run single sample analysis
# nextflow run main.nf --input /home/chbope/extension/nWGS_manuscript_data/data/results/epi2me/modkit/T001.wf_mods.bedmethyl.gz --metadata metadata.csv --outdir results_single

echo ""
echo "Multi-sample analysis with CSV:"
echo "Command: nextflow run main.nf --input samples.csv --metadata metadata.csv --outdir results_multi"

# Uncomment the line below to run multi-sample analysis
# nextflow run main.nf --input samples.csv --metadata metadata.csv --outdir results_multi

echo ""
echo "Test mode (quick execution):"
echo "Command: nextflow run main.nf --input samples.csv --metadata metadata.csv --outdir results_test -profile test"

# Uncomment the line below to run test mode
# nextflow run main.nf --input samples.csv --metadata metadata.csv --outdir results_test -profile test

echo ""
echo "To run with conda environment:"
echo "nextflow run main.nf --input samples.csv --metadata metadata.csv --outdir results_conda -profile conda"

echo ""
echo "To run on SLURM cluster:"
echo "nextflow run main.nf --input samples.csv --metadata metadata.csv --outdir results_slurm -profile slurm"

echo ""
echo "To show help:"
echo "nextflow run main.nf --help"

echo ""
echo "Choose your execution mode and uncomment the appropriate line in this script."
echo "Or run the commands manually."

