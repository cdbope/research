# DeepSomatic Pipeline Configuration File
# Edit these paths according to your system setup

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
SAMPLE_ID="T21-140"
BAM_FILE="T21-140.occ.bam"

# DeepSomatic model type
# Options: ONT_TUMOR_ONLY, FFPE_WGS_TUMOR_ONLY, FFPE_WES_TUMOR_ONLY, etc.
MODEL_TYPE="ONT_TUMOR_ONLY"

# Use default PON filtering (true/false)
USE_PON_FILTERING=true

# Number of parallel shards (default: use all CPU cores)
NUM_SHARDS=$(nproc)
