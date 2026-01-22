BIN_VERSION="1.9.0"
INPUT_DIR="/home/chbope/extension/nWGS_manuscript_data/data/testdata/single_bam_folder/merge_bam"
OUTPUT_DIR="/home/chbope/extension/data/200GMBs/pcgr/starsigndna/refit/tmp"
REF_DIR="/home/chbope/extension/nWGS_manuscript_data/data/reference"

sudo docker pull google/deepsomatic:"${BIN_VERSION}"

sudo docker run \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
-v ${REF_DIR}:${REF_DIR} \
google/deepsomatic:"${BIN_VERSION}" \
run_deepsomatic \
--model_type=FFPE_WGS_TUMOR_ONLY \
--ref=/home/chbope/extension/nWGS_manuscript_data/data/reference/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa  \
--reads_tumor=/home/chbope/extension/data/t25-152/T25-152.occ.bam  \
--output_vcf=${OUTPUT_DIR}/T25-152_deepsomatic_output.vcf.gz \
--sample_name_tumor="T25-152" \
--num_shards=$(nproc) \
--logging_dir=${OUTPUT_DIR}/T25-152 \
--intermediate_results_dir=${OUTPUT_DIR}/T25-152 \
--use_default_pon_filtering=true \
#--regions=chr1
