# DeepSomatic Pipeline v3 Configuration File - Batch Processing
# Edit these paths according to your system setup

# DeepSomatic Docker version
BIN_VERSION="1.9.0"

# Input/Output directories
INPUT_DIR="/home/chbope/extension/nWGS_manuscript_data/data/testdata/single_bam_folder/merge_bam"
OUTPUT_DIR="/home/chbope/extension/data/200GMBs/pcgr/starsigndna/refit/tmp"
REF_DIR="/home/chbope/extension/nWGS_manuscript_data/data/reference"

# ANNOVAR paths
ANNOVAR_DIR="/home/chbope/Documents/nanopore/nWGS_manuscript/nWGS_pipeline_routine/bin"
HUMANDB_DIR="/home/chbope/extension/nWGS_manuscript_data/data/humandb"

# Reference genome file
REF_GENOME="GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa"

# Sample list file (one sample ID per line)
# Example content:
#   T25-152
#   T001
#   T002
SAMPLE_LIST_FILE="${INPUT_DIR}/sample_id.txt"

# BAM file suffix (will be appended to sample ID to form BAM filename)
# Example: If SAMPLE_ID="T25-152" and BAM_SUFFIX=".occ.bam", BAM file will be "T25-152.occ.bam"
BAM_SUFFIX=".occ.bam"

# DeepSomatic model type
# Options: ONT_TUMOR_ONLY, FFPE_WGS_TUMOR_ONLY, FFPE_WES_TUMOR_ONLY, etc.
MODEL_TYPE="ONT_TUMOR_ONLY"

# Use default PON filtering (true/false)
USE_PON_FILTERING=true

# Number of parallel shards (default: use all CPU cores)
NUM_SHARDS=$(nproc)
