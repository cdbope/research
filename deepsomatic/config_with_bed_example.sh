#!/bin/bash
# DeepSomatic Pipeline Configuration File with BED Support
# Updated from your existing config.sh with BED file options added

# DeepSomatic Docker version
BIN_VERSION="1.9.0"

# Input/Output directories
INPUT_DIR="/home/chbope/extension/nWGS_manuscript_data/data/testdata/single_bam_folder/merge_bam"
OUTPUT_DIR="/home/chbope/extension/script/deepsomatic/output"
REF_DIR="/home/chbope/extension/nWGS_manuscript_data/data/reference"

# ANNOVAR paths
ANNOVAR_DIR="/home/chbope/Documents/nanopore/nWGS_manuscript/nWGS_pipeline_routine/bin"
HUMANDB_DIR="/home/chbope/extension/nWGS_manuscript_data/data/humandb"

# Reference genome file
REF_GENOME="GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa"

# Sample information (v2 - Single sample mode)
SAMPLE_ID="T25-256"
BAM_FILE="T25-256.occ.bam"

# DeepSomatic model type
# Options: ONT_TUMOR_ONLY, FFPE_WGS_TUMOR_ONLY, FFPE_WES_TUMOR_ONLY, etc.
MODEL_TYPE="ONT_TUMOR_ONLY"

# Use default PON filtering (true/false)
USE_PON_FILTERING=true

# Number of parallel shards (default: use all CPU cores)
NUM_SHARDS=$(nproc)

# ============================================================================
# BED FILE CONFIGURATION (NEW - FOR REGION-SPECIFIC VARIANT CALLING)
# ============================================================================
# Add BED file path here to limit variant calling to specific regions
# This makes DeepSomatic call variants only in the regions specified in the BED file
# (similar to ClairS-TO's --bed_fn parameter)

# Option 1: Use same BED as ClairS-TO for fair comparison
#           This is the protein coding BED file used in your benchmark
BED_FILE="/home/chbope/extension/nWGS_manuscript_data/data/reference/OCC.protein_coding.3col.bed"

# Option 2: No BED file - genome-wide calling
#           Uncomment the line below and comment out the line above
# BED_FILE=""

# Option 3: Custom cancer gene panel
#           Replace with your custom BED file path
# BED_FILE="/path/to/cancer_genes_panel.bed"

# ============================================================================
# INSTRUCTIONS:
# ============================================================================
# 1. To enable BED filtering: Set BED_FILE to a valid path (Option 1 is active now)
# 2. To disable BED filtering: Comment out BED_FILE or set to empty string
# 3. The BED file above matches your ClairS-TO analysis for fair comparison
#
# Expected results with BED file enabled:
#   - DeepSomatic: ~35-45 variants (BED regions only)
#   - ClairS-TO:   ~33 variants (BED regions only)
#   - Much closer comparison than genome-wide (66 vs 33)
#
# To use this config: Copy to config.sh or use with deepsomatic_v2_with_bed.sh
# ============================================================================
