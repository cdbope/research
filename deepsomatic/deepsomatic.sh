BIN_VERSION="1.9.0"
# Set up input and output directory data
INPUT_DIR="/home/chbope/extension/nWGS_manuscript_data/data/testdata/single_bam_folder/merge_bam"
OUTPUT_DIR="/home/chbope/extension/data/200GMBs/pcgr/starsigndna/refit/tmp"
REF_DIR="/home/chbope/extension/nWGS_manuscript_data/data/reference"

# Clean up any previous intermediate files to avoid corruption
echo "Cleaning up previous intermediate files..."
sudo rm -rf ${OUTPUT_DIR}/intermediate_results_dir/*

# Check if Docker image exists locally
echo "Checking DeepSomatic Docker image..."
if ! sudo docker image inspect google/deepsomatic:"${BIN_VERSION}" > /dev/null 2>&1; then
    echo "Image not found locally. Pulling DeepSomatic Docker image..."
    sudo docker pull google/deepsomatic:"${BIN_VERSION}"
else
    echo "DeepSomatic image ${BIN_VERSION} found locally. Skipping pull."
fi

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}":ro \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${REF_DIR}":"${REF_DIR}":ro \
google/deepsomatic:"${BIN_VERSION}" \
run_deepsomatic \
--model_type=ONT_TUMOR_ONLY \
--ref="${REF_DIR}/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa" \
--reads_tumor="${INPUT_DIR}/T25-152.occ.bam" \
--output_vcf="${OUTPUT_DIR}/T25-152_ont_deepsomatic_output.vcf.gz" \
--sample_name_tumor="T25-152" \
--num_shards=$(nproc) \
--logging_dir="${OUTPUT_DIR}/logs" \
--intermediate_results_dir="${OUTPUT_DIR}/intermediate_results_dir" \
--use_default_pon_filtering=true
#--regions=chr1
