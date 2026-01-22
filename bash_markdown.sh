#!/bin/bash

# Define the file containing the sample IDs
samples_file="/home/chbope/extension/trash/t23-057/t23-057.txt"

# RMarkdown template file path
rmd_template="/home/chbope/extension/script/nextflow_markdown_pipeline_update_final.Rmd"

# Loop through each sample ID in the samples file
while read -r sample_id tumor_content; do
    echo "Processing sample: $sample_id"

    # Build dynamic input paths for this sample
    craminoreport="/home/chbope/extension/trash/t23-057/${sample_id}_cramino_statistics.txt"
    sample_ids_file="/home/chbope/extension/trash/t23-057/t23-057.txt"
    nanodx="/home/chbope/extension/trash/t23-057/${sample_id}_nanodx_classifier.tsv"
    dictionaire="/media/chbope/Manga/images_reference/reference/nanoDx/static/Capper_et_al_dictionary.txt"
    logo="/home/chbope/extension/trash/t23-057/log_update.pdf"
    cnv_plot="/home/chbope/extension/trash/t23-057/${sample_id}_cnv_plot_full.pdf"
    tumor_number="/home/chbope/extension/trash/t23-057/${sample_id}_tumor_copy_number.txt"
    annotatecnv="/home/chbope/extension/trash/t23-057/${sample_id}_annotatedcnv_filter_header.csv"
    cnv_chr9="/home/chbope/extension/trash/t23-057/${sample_id}_cnv_chr9.pdf"
    cnv_chr7="/home/chbope/extension/trash/t23-057/${sample_id}_cnv_chr7.pdf"
    mgmt_results="/home/chbope/extension/trash/t23-057/${sample_id}_MGMT_results.csv"
    merge_results="/home/chbope/extension/trash/t23-057/${sample_id}_merge_annotation_filter_snvs_allcall_filter.csv"
    #fusion_events="/home/chbope/extension/trash/t23-057/${sample_id}_filter_fusion_event.tsv"
    fusion_events="/home/chbope/extension/trash/test2/test3/T22-024_test_filter_fusion_event23juneall.tsv"
    terphtml="/home/chbope/extension/trash/t23-057/${sample_id}_tertp_id1.html"
    #annotsvhtml="/data/pipeline/results/structure_variant/annotsv/${sample_id}_annotated_variants.html"
    svannahtml="/home/chbope/extension/trash/t23-057/${sample_id}_occ_svanna_annotation.html"
    egfr_coverage="/home/chbope/extension/trash/t23-057/${sample_id}_egfr_coverage.pdf"
    idh1_coverage="/home/chbope/extension/trash/t23-057/${sample_id}_idh1_coverage.pdf"
    tertp_coverage="/home/chbope/extension/trash/t23-057/${sample_id}_tertp_coverage.pdf"
    
    # Output PDF path
    output_file="/home/chbope/extension/trash/t23-057/${sample_id}_markdown_pipeline_report_final.pdf"

    # Now call the Rscript exactly as you want
    Rscript -e "rmarkdown::render('${rmd_template}', output_file=commandArgs(trailingOnly=TRUE)[20])" \
      "${sample_id}" \
      "${craminoreport}" \
      "${sample_ids_file}" \
      "${nanodx}" \
      "${dictionaire}" \
      "${logo}" \
      "${cnv_plot}" \
      "${tumor_number}" \
      "${annotatecnv}" \
      "${cnv_chr9}" \
      "${cnv_chr7}" \
      "${mgmt_results}" \
      "${merge_results}" \
      "${fusion_events}" \
      "${terphtml}" \
      "${svannahtml}" \
      "${egfr_coverage}" \
      "${idh1_coverage}" \
      "${tertp_coverage}" \
      "${output_file}"
    rm -rf /tmp/Rtmp*

    echo "Finished sample: $sample_id"

done < "${samples_file}"

