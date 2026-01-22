#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def start_time = new Date()
//params.mgmt_promoter_r_script = "mnt/scripts/MGMT_Prospective2.R"
// Define the base path as a parameter
// params.path = "/home/chbope/extension"
// params.input_path = "${params.path}"
// params.output_path = "${params.path}/results"
// params.annotate_dir = "/home/chbope/extension/data/annotations"
// params.config_dir = "/home/chbope/Documents/nanopore/packages/knotannotsv/knotAnnotSV"
// params.ref_dir ="/home/chbope/extension/data/reference"
// params.model_path="/home/chbope/extension/Data_for_Bope"
// params.clair3_dir="/home/chbope/Documents/nanopore/Data_for_Bope/results/sample_id1/callvariantclair3/clair3_output"
// params.humandb_dir="/home/chbope/extension/data/annovar/humandb"
// params.clairSTo_dir="/home/chbope/Documents/nanopore/Data_for_Bope/results/sample_id1/callvariantclairsto/clairsto_output"
// params.svanna_dir="/home/chbope/extension/data/svanna-cli-1.0.4/"
// params.bin_dir="/home/chbope/Documents/nanopore/nextflow/bin/"
// params.epi2me_dir="/home/chbope/Documents/nanopore/epi2me/wf-human-variation-master/"
// params.out_dir_epi2me ="/home/chbope/extension/out_dir_epi2me"

// // input reference files for all samples

// params.reference_genome = file("${params.ref_dir}/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa")
// params.reference_genome_bai = file("${params.ref_dir}/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa.fai")

// params.occ_fusions=file("${params.ref_dir}/OCC.fusions.bed")
// params.epicsites=file("${params.ref_dir}/EPIC_sites_NEW.bed")
// params.mgmt_cpg_island_hg38=file("${params.ref_dir}/MGMT_CpG_Island.hg38.bed")

// params.occ_snv_screening=file("${params.ref_dir}/OCC.SNV.screening.bed")
// params.tertp_variants=file("${params.ref_dir}/TERTp_variants.bed")
// params.ncbirefseq=file("${params.ref_dir}/ncbiRefSeq.txt.gz")
// params.refgene=file("${params.humandb_dir}/hg38_refGene.txt")
// params.hg38_refgenemrna=file("${params.humandb_dir}/hg38_refGeneMrna.fa")
// params.clinvar=file("${params.humandb_dir}/hg38_clinvar_20240611.txt")
// params.clinvarindex=file("${params.humandb_dir}/hg38_clinvar_20240611.txt.idx")
// params.hg38_cosmic100=file("${params.humandb_dir}/hg38_cosmic100coding2024.txt")
// params.hg38_cosmic100index=file("${params.humandb_dir}/hg38_cosmic100coding2024.txt.idx")
// params.knotannotsv_conf=file("${params.config_dir}/config_AnnotSV.yaml")
// params.occ_fusion_genes_list=file("${params.ref_dir}/occ_fusions_genes.txt")
// //params.bins_bed=file("/home/chbope/Documents/nanopore/Data_for_Bope/segs/T21-058_bins.bed")
// params.convertvariant_config = file("${params.ref_dir}/annotsv2.json")
// params.occ_genes = file("${params.ref_dir}/OCCgenes.rds")
// params.vcf2circos_json = file("${params.ref_dir}/options.json")
// params.sturgeon_model = file("${params.ref_dir}/general.zip")
// params.hg19_450model = file("${params.ref_dir}/nanoDx/static/hg19_450model.bed")
// params.nanodx_450model = file("${params.ref_dir}/nanoDx/static/Capper_et_al_NN.pkl")
// params.nanodx_dictinaire = file("${params.ref_dir}/nanoDx/static/Capper_et_al_dictionary.txt")
// params.nn_model = file("${params.ref_dir}/nanoDx/workflow/envs/NN_model.yaml")
// params.snakefile_nanodx = file("${params.ref_dir}/nanoDx/workflow/Snakefile")
// params.sub_snakefile = file("${params.ref_dir}/nanoDx/workflow/methylation.snakefile")
// params.tr_bed_file = file("${params.ref_dir}/human_GRCh38_trf.bed")
// params.epi2me_config_file = file("/home/chbope/Documents/nanopore/epi2me/wf-human-variation-master/nextflow.config")
// params.epi2me_base_config_file = file("/home/chbope/Documents/nanopore/epi2me/wf-human-variation-master/base.config")
   params.sample_id_file = file("/data/pipeline/data/200_gbm_remove_duplicate.txt")
// params.mardown_logo = file("/home/chbope/Downloads/log_update.pdf")
// //params.ma = file("/home/chbope/extension/Data_for_Bope/sample_ids.txt")
// //input individual samples files

// //input file folder
// merge_bam_folder='/home/chbope/extension/Data_for_Bope/merge_bams/'
// vcf_folder='/home/chbope/extension/Data_for_Bope/vcfs/'
// bams_folder='/home/chbope/extension/Data_for_Bope/bams/'
// bedmethyl_folder='/home/chbope/extension/Data_for_Bope/bedmethyl/'
// segsfromepi2me_folder='/home/chbope/extension/Data_for_Bope/segs/'
// sv_folder='/home/chbope/extension/Data_for_Bope/sv/'
// nanodx_workflow_dir='/home/chbope/Documents/nanopore/packages/nanoDx/workflow'


//params.mgmt_promoter_r_script = "mnt/scripts/MGMT_Prospective2.R"
// Define the base path as a parameter
params.path = "/data/pipeline"
params.input_path = "${params.path}"
//params.output_path = "/home/prom/africa/angola/200_GBMs_2025_analysis/results"
params.output_path = "/home/prom/africa/angola/P24_WGS_analysis/results"
params.annotate_dir = "/data/pipeline/data/annotations"
params.config_dir = "/data/pipeline/packages/knotannotsv/knotAnnotSV"
params.ref_dir ="/data/pipeline/data/reference"
params.model_path="/data/pipeline/Data_for_Bope"
params.clair3_dir="/data/pipeline/results/callvariantclair3/clair3_output"
params.humandb_dir="/data/pipeline/data/annovar/humandb"
params.clairSTo_dir="/data/pipeline/callvariantclairsto/clairsto_output"
params.svanna_dir="/data/pipeline/data/svanna-cli-1.0.4/"
params.bin_dir="/data/pipeline/nextflow/bin/"
params.out_dir_epi2me ="/data/pipeline/out_dir_epi2me"

// input reference files for all samples

params.reference_genome = file("${params.ref_dir}/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa")
params.reference_genome_bai = file("${params.ref_dir}/GCF_000001405.39_GRCh38.p13_genomic_chr_only_plus_mt.fa.fai")

params.occ_fusions=file("${params.ref_dir}/OCC.fusions.bed")
params.epicsites=file("${params.ref_dir}/EPIC_sites_NEW.bed")
params.mgmt_cpg_island_hg38=file("${params.ref_dir}/MGMT_CpG_Island.hg38.bed")

params.occ_snv_screening=file("${params.ref_dir}/OCC.SNV.screening.bed")
params.tertp_variants=file("${params.ref_dir}/TERTp_variants.bed")
params.ncbirefseq=file("${params.ref_dir}/ncbiRefSeq.txt.gz")
params.refgene=file("${params.humandb_dir}/hg38_refGene.txt")
params.hg38_refgenemrna=file("${params.humandb_dir}/hg38_refGeneMrna.fa")
params.clinvar=file("${params.humandb_dir}/hg38_clinvar_20240611.txt")
params.clinvarindex=file("${params.humandb_dir}/hg38_clinvar_20240611.txt.idx")
params.hg38_cosmic100=file("${params.humandb_dir}/hg38_cosmic100coding2024.txt")
params.hg38_cosmic100index=file("${params.humandb_dir}/hg38_cosmic100coding2024.txt.idx")
params.knotannotsv_conf=file("${params.config_dir}/config_AnnotSV.yaml")
params.occ_fusion_genes_list=file("${params.ref_dir}/occ_fusions_genes.txt")
//params.bins_bed=file("/home/chbope/Documents/nanopore/Data_for_Bope/segs/T21-058_bins.bed")
params.convertvariant_config = file("${params.ref_dir}/annotsv2.json")
params.occ_genes = file("${params.ref_dir}/OCCgenes.rds")
params.vcf2circos_json = file("${params.ref_dir}/options.json")
params.sturgeon_model = file("${params.ref_dir}/general.zip")
params.hg19_450model = file("${params.ref_dir}/nanoDx/static/hg19_450model.bed")
params.nanodx_450model = file("${params.ref_dir}/nanoDx/static/Capper_et_al_NN.pkl")
params.nn_model = file("${params.ref_dir}/nanoDx/workflow/envs/NN_model.yaml")
params.snakefile_nanodx = file("${params.ref_dir}/nanoDx/workflow/Snakefile")
//params.sub_snakefile = file("${params.ref_dir}/nanoDx/workflow/methylation.snakefile")
params.tr_bed_file = file("${params.ref_dir}/human_GRCh38_trf.bed")
params.epi2me_config_file = file("/home/chbope/Documents/nanopore/epi2me/wf-human-variation-master/nextflow.config")
params.epi2me_base_config_file = file("/home/chbope/Documents/nanopore/epi2me/wf-human-variation-master/base.config")
params.nanodx_dictinaire = file("${params.ref_dir}/nanoDx/static/Capper_et_al_dictionary.txt")
params.mardown_logo = file("/data/pipeline/data/log_update.pdf")
params.gviz_data = file("${params.ref_dir}/Gviz.RData")
params.cytoband_file = file("${params.ref_dir}/cytoBandIdeo.hg38.txt")
params.genecode_bed = file("${params.ref_dir}/gencode.v48.annotation.gff3")
//input individual samples files

//input file folder
//merge_bam_folder='/home/prom/africa/angola/200_GBMs_2025_analysis/epicnv/qdna_seq'
cnv_rds ='/home/prom/africa/angola/200_GBMs_2025_analysis/epicnv/qdna_seq'


//merge_bam_folder = '/home/prom/africa/angola/200_GBMs_merge_bam/repeat_double_cell/remove_duplicate/occ_bam'
merge_bam_folder = '/home/prom/africa/kenya/P24_WGS_merged_BAM'
//vcf_folder='/home/prom/africa/angola/200_GBMs_2025_analysis/episv/'
//vcf_folder='/data/pipeline/trash/data/sv/'
//vcf_folder = '/data/pipeline/trash/data/sv/'
bams_folder='/home/prom/africa/angola/P24_WGS_merged_BAM/P24_project_files/OCC.bam'
//bams_folder='/home/prom/africa/kenya/P24_WGS_merged_BAM/P24_project_files/OCC.bam'
///bams_folder = '/home/prom/africa/angola/OCC_ROI_BAM/'
bedmethyl_folder='/home/prom/africa/angola/P24_WGS_merged_BAM/P24_project_files/wf_mods.bedmethyl.gz'
//bedmethyl_folder='/home/prom/africa/angola/200_GBMs_2025_analysis/epimodkit'
//bedmethyl_folder = '/data/pipeline/trash/data/bedmethyl/'
///segsfromepi2me_folder='/home/prom/africa/angola/200_GBMs_2025_analysis/epicnv/qdna_seq/'
segsfromepi2me_folder='/home/prom/africa/angola/P24_WGS_merged_BAM/P24_project_files/cnv'
//segsfromepi2me_folder = '/data/pipeline/trash/epi2me/output/epicnv/qdna_seq'
//sv_folder='/home/prom/africa/angola/200_GBMs_2025_analysis/episv/'
sv_folder='/home/prom/africa/angola/P24_WGS_merged_BAM/P24_project_files/Sniffles2_SV'
///sv_folder = '/data/pipeline/trash/data/sv/'
nanodx_workflow_dir='/data/pipeline/packages/nanoDx/workflow'
nanodx_bed_folder = '/data/pipeline/results/classifier/nanodx'

// Define the available run modes
def mode_descriptions = [
    'methylation': 'Run methylation processes.',
    'annotsv': 'Run structure variant annotation process (AnnotSV and Svanna).',
    'cnv' : 'Run a copy number variation process.', 
    'occ' : 'Run OCC processes (Clair3 and Clais_to).', 
    'terp' : 'Run igv_tool process.',
    'rmd'   : 'Run the rmd report',
    'stat' : 'Generate Sequencing stats',
    'all': 'Run all the pipeline processes.'
    ]

// Define the run_mode parameter (default is "both")
params.run_mode = params.run_mode ?: 'all'


// Show available modes and descriptions before starting the pipeline
println "Available run modes with descriptions:"
mode_descriptions.each { mode, description ->
    println "- ${mode}: ${description}"
}
println "Selected run mode: ${params.run_mode}"

if (!mode_descriptions.keySet().contains(params.run_mode)) {
    println "ERROR: Invalid run mode '${params.run_mode}'. Please select from ${mode_descriptions.keySet().toList().join(', ')}."
    exit 1  // Exit with error status
}


//##############################3
//## Extract overlapping EPIC sites for methylation based classification 
//#############################

//# wf_mods_bedmethyl generated by modkit in Epi2me-labs on the P24
//# EPIC_sites_NEW.bed contains chromosomal coordinates and names of EPIC 850 probes
//# Uses Bedtools intersect, awk and sed


//epi2me pipeline to run methyl, cnv and sv



process extract_epic {
    cpus 2
    memory '2 GB'
    label 'epic'
    tag "${sample_id}"

    publishDir "${params.output_path}/methylation/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), file(bedmethyl), file(epicsites), file(mgmt_cpg_island_hg38)

    output:
   // path("${sample_id}_EpicSelect.bed")
    tuple val(sample_id), path("${sample_id}_EpicSelect_header.bed"), emit: epicselectnanodxinput
    path("${sample_id}_MGMT.bed")
    tuple val(sample_id), path("${sample_id}_MGMT_header.bed"), emit: MGMTheaderout
    tuple val(sample_id), path("${sample_id}_wf_mods.bedmethyl_intersect.bed"), emit: sturgeonbedinput

    script:
    """
    # Step 1: Run intersectBed and save the output to EpicSelect.bed
    which intersectBed 
    intersectBed -a $bedmethyl -b $epicsites -wb | \
    awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" }{print \$1,\$2,\$3,\$4,\$5,\$11,\$23}' > ${sample_id}_EpicSelect.bed

    intersectBed -a $bedmethyl  -b $epicsites   -wb | awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" } {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17, \$18}' > ${sample_id}_wf_mods.bedmethyl_intersect.bed

    # Step 2: Add header to EpicSelect.bed
    awk 'BEGIN {print "Chromosome\\tStart\\tEnd\\tmodBase\\tCoverage\\tMethylation_frequency\\tIllumina_ID"} 1' ${sample_id}_EpicSelect.bed > ${sample_id}_EpicSelect_header.bed

    # Step 3: Run intersectBed with MGMTp_hg38
    intersectBed -a $bedmethyl  -b $mgmt_cpg_island_hg38 | \
    awk -v OFS="\\t" '\$1=\$1' | awk -F'\\t' 'BEGIN{ OFS="\\t" }{print \$1,\$2,\$3,\$4,\$5,\$11,\$12,\$13,\$14,\$15,\$16}'  > ${sample_id}_MGMT.bed

    # Step 4: Add header to MGMT.bed
    awk 'BEGIN {print "Chrom\\tStart\\tEnd\\tmodBase\\tDepth\\tMethylation\\tNmod\\tNcanon\\tNother\\tNdelete\\tNfail"} 1' ${sample_id}_MGMT.bed > ${sample_id}_MGMT_header.bed
    """
}

//Sturgeon classifier
process sturgeon {
    cpus 2
    memory '2 GB'
    label 'epic'
    publishDir "${params.output_path}/classifier/sturgeon", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(sturgeon_bed), path(sturgeon_model)

    output:
    tuple path("${sample_id}_bedmethyl_sturgeon.bed"), path("${sample_id}_bedmethyl_sturgeon")


    """
    /sturgeon/venv/bin/sturgeon inputtobed -i $sturgeon_bed  -o ${sample_id}_bedmethyl_sturgeon.bed  -s modkit_pileup  --reference-genome hg38
   
    /sturgeon/venv/bin/sturgeon predict -i ${sample_id}_bedmethyl_sturgeon.bed   -o  ${sample_id}_bedmethyl_sturgeon --model-files $sturgeon_model  --plot-results

    """
}

//NanoDX classifier

process nanodx {
    cpus 2
    memory '16 GB'
    label 'epic'
    publishDir "${params.output_path}/classifier/nanodx", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(nanodx_bed), path(hg19_450model)

    output:
    tuple val(sample_id), path("${sample_id}_nanodx_bedmethyl.bed"), emit: nanodx450out

    script:
    """
    nanodx450intersectdataframe.py $hg19_450model $nanodx_bed ${sample_id}_output_cpg.bed ${sample_id}_nanodx_bedmethyl.bed ${sample_id}_nanodx_bedmethylfilter.bed
    """
}

//nanodx classifier for 45
process run_nn_classifier_desktop {
    label 'snakemake'
    publishDir "${params.output_path}/classifier/nanodx", mode: 'copy'
    
    cpus 2
    memory '12 GB'

    input:
    tuple val(sample_id), path(bed_file), path(model_file), path(nn_model) //, path(snakefile_nanodx)

    output:
    tuple val(sample_id), path("${sample_id}_nanodx_classifier.txt")
    tuple val(sample_id), path("${sample_id}_nanodx_classifier.tsv"), emit:rmdnanodx

    containerOptions = "-B /data/pipeline/trash/tmp:/data/pipeline/trash/tmp"

    script:
    """
    # Set TMPDIR to the provided path
    export TMPDIR="/data/pipeline/trash/tmp/"
    mkdir -p \$TMPDIR

    # Verify that the directory is writable
    echo "Testing write permissions in \$TMPDIR"
    touch \$TMPDIR/test_file
    if [ -f \$TMPDIR/test_file ]; then
        echo "Write permissions are working."
    else
        echo "Write permissions are NOT working."
    fi

    # Source the Conda environment
    source /opt/conda/etc/profile.d/conda.sh
    #conda info --envs
    #conda env list 
    conda activate

    echo \$CONDA_PREFIX
    which conda
    
    # Verify environment setup
    which python
    python -c "import pandas; print(pandas.__version__)"
    
    # Create a temporary Snakefile with the dynamic inputs
    cat << EOF > Snakefile
rule all:
    input:
        "${sample_id}_nanodx_classifier.txt",
        "${sample_id}_nanodx_classifier.tsv"

rule NN_classifier:
    input:
        bed = "${bed_file}",
        model = "${model_file}"
    output:
        txt = "${sample_id}_nanodx_classifier.txt",
        votes = "${sample_id}_nanodx_classifier.tsv"
    threads: 4
    resources: 
        mem_mb = 16384
    script: "${nanodx_workflow_dir}/scripts/classify_NN_bedMethyl.py"
EOF

    # Run snakemake within the container
    snakemake --use-conda \
        --conda-prefix \$TMPDIR/.snakemake/conda \
        --cores ${task.cpus} \
        --verbose \
        NN_classifier
    """
}


process run_nn_classifier {
    label 'snakemake'
    publishDir "${params.output_path}/classifier/nanodx", mode: 'copy'
    
    cpus 2
    memory '12 GB'

    input:
    tuple val(sample_id), path(bed_file), path(model_file), path(nn_model) //, path(snakefile_nanodx), 

    output:
    tuple val(sample_id), path("${sample_id}_nanodx_classifier.txt"), path("${sample_id}_nanodx_classifier.tsv")
    tuple val(sample_id), path("${sample_id}_nanodx_classifier.tsv"), emit:rmdnanodx

    //containerOptions = "-B /data/pipeline/trash/tmp:/data/pipeline/trash/tmp"

    script:
    """
    # Set up all necessary environment variables and directories
    export TMPDIR="/data/pipeline/trash/tmp"
    export CONDA_PKGS_DIRS="\$TMPDIR/conda_pkgs"
    export CONDA_ENVS_PATH="\$TMPDIR/conda_envs"
    export HOME="\$TMPDIR/home"
    export XDG_CACHE_HOME="\$TMPDIR/cache"
    export XDG_DATA_HOME="\$TMPDIR/local"

    # Create all necessary directories
    mkdir -p \$TMPDIR \$CONDA_PKGS_DIRS \$CONDA_ENVS_PATH \$HOME \$XDG_CACHE_HOME \$XDG_DATA_HOME

    # Verify that the directories are writable
    echo "Testing write permissions"
    touch \$TMPDIR/test_file
    if [ -f \$TMPDIR/test_file ]; then
        echo "Write permissions are working."
    else
        echo "Write permissions are NOT working."
    fi

   source /opt/conda/etc/profile.d/conda.sh
   # conda env list 
    conda activate nanodx_env2feb
    #conda activate nanodx_env


    echo \$CONDA_PREFIX
   # which conda
    
    # Verify environment setup
    #which python
    python -c "import pandas; print(pandas.__version__)"
    
    # Create a temporary Snakefile with the dynamic inputs
    cat << EOF > Snakefile
rule all:
    input:
        "${sample_id}_nanodx_classifier.txt",
        "${sample_id}_nanodx_classifier.tsv"

rule NN_classifier:
    input:
        bed = "${bed_file}",
        model = "${model_file}"
    output:
        txt = "${sample_id}_nanodx_classifier.txt",
        votes = "${sample_id}_nanodx_classifier.tsv"
    threads: 4
    resources: 
        mem_mb = 16384
    script: "${nanodx_workflow_dir}/scripts/classify_NN_bedMethyl.py"
EOF
    # Run snakemake within the container
    snakemake --use-conda \
        --conda-prefix \$TMPDIR/.snakemake/conda \
        --cores ${task.cpus} \
        --verbose \
        NN_classifier
    """
}
// Extract MGMT promoter Ccd pG-island #extracting MGMT prmoter sites
process mgmt_promoter {
    cpus 2
    memory '2 GB'
    label 'epic'
    publishDir "${params.output_path}/methylation/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(mgmt_bed)

    output:
    tuple val(sample_id), path("${sample_id}_MGMT_results.csv"), emit:mgmtresultsout
    
    script:
    MGMT_results = "${sample_id}_MGMT_results.csv"

  //# calculate methylation in MGMT promoter

  """

   MGMT_Prospective2.R ${mgmt_bed} ${sample_id}_MGMT_results.csv

  """
}

//####################################
//## Annotate structural variants
//####################

//# _wf_sv.vcf.gz greated by Sniffles2 via Epi2melabs on the P24
//# Fusion gene candidate regions in OCC.fusions.bed
//# Variants annotated with AnnotSV

process svannasv {

   label 'svannasv'
   cpus 2
   memory '2 GB'
   publishDir "${params.output_path}/structure_variant/svannasv/", mode: "copy", overwrite: true

   input:
   tuple val(sample_id), path(wf_sv), path(wf_sv_tbi),path(occ_fusions)

   output:
   //file("${sample_id}_OCC_SVs.vcf")
   tuple val(sample_id), path("${sample_id}_OCC_SVs.vcf"), emit: occsvannavcfout
   tuple val(sample_id), path("${sample_id}_occ_svanna_annotation.html"), emit:rmdsvannahtml 
   tuple val(sample_id), path("${sample_id}_occ_svanna_annotation.vcf.gz"), emit: occsvannaannotationannotationvcf

//   export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
//   export PATH=$JAVA_HOME/bin:$PATH
   script:
   """

   intersectBed -a $wf_sv  -b $occ_fusions  -header > ${sample_id}_OCC_SVs.vcf

   # Check if intersection file is empty (excluding header)
   if [ \$(grep -v '^#' ${sample_id}_OCC_SVs.vcf | wc -l) -eq 0 ]; then
      # If empty, use original wf_sv file
      INPUT_FILE=$wf_sv
   else
      # If not empty, use intersection file
      INPUT_FILE=${sample_id}_OCC_SVs.vcf
   fi
 
   /usr/lib/jvm/java-17-openjdk-amd64/bin/java -jar ${params.bin_dir}/svanna-cli-1.0.3.jar prioritize  \
   -d ${params.svanna_dir}/svanna-data  \
   --vcf \$INPUT_FILE \
   --phenotype-term HP:0100836 \
   --output-format html,vcf \
   --prefix ${sample_id}_occ_svanna_annotation

  # cp "${sample_id}_occ_svanna_annotation.html" "${params.output_path}/report/${sample_id}_svanna.html"

   """
}


process svannasv_fusion_events {
    cpus 4
    memory '2 GB'
    label 'svannasv'
    publishDir "${params.output_path}/structure_variant/svannasv/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(occ_svannavcf), path(genecode_bed), path(occ_fusions_genes)

    output:
    tuple val(sample_id), path("${sample_id}_filter_fusion_event.tsv"), emit: filterfusioneventout

    script:

    """
    breaking_point_bed_translocation.py --vcf $occ_svannavcf --out  ${sample_id}_breaking_bedpoints.bed

    awk 'BEGIN{OFS="\t"} {if (\$1 !~ /^chr/) \$1 = "chr"\$1; print}' ${sample_id}_breaking_bedpoints.bed > ${sample_id}_breaking_bedpoints_sort.bed

    intersectBed -a ${sample_id}_breaking_bedpoints_sort.bed  -b $genecode_bed  -wb  > ${sample_id}_breaking_bedpoints_genecode.bed

    #remove duplicate bed points

    remove_duplicate_bed_awk1.py --in  ${sample_id}_breaking_bedpoints_genecode.bed  \
            --formatted ${sample_id}_breaking_bedpoints_genecode_format.bed \
             --out ${sample_id}_breaking_bedpoints_genecode_clean.bed    \
             --paired ${sample_id}_breaking_bedpoints_genecode_clean_paired.bed  \
             --gene-list $occ_fusions_genes \
             --filtered ${sample_id}_filter_fusion_event.tsv

    """
}

process annotesvoriginal {
    cpus 8
    memory '20 GB'
   label 'annotation'
   publishDir "${params.output_path}/structure_variant/annotsv", mode: "copy", overwrite: true

   input:
        tuple val(sample_id), path(wf_sv), path(wf_sv_tbi), path(occ_fusions), path(occ_fusion_genes_list)

   output:
        path("${sample_id}_OCC_SVs.vcf")
        tuple val(sample_id), path("annotated_variants.tsv"), emit: annotatedvariantsout
        tuple val(sample_id), path("${sample_id}_annotSV_fusion_extraction.csv"), emit: annotsvfusion
        path("${sample_id}_annotated_variants.tsv")

   script:
   """
   export TMPDIR=\$PWD/tmp_dir
   mkdir -p \$TMPDIR

      # Check if intersection file is empty (excluding header)
   if [ \$(grep -v '^#' ${sample_id}_OCC_SVs.vcf | wc -l) -eq 0 ]; then
      # If empty, use original wf_sv file
      INPUT_FILE=$wf_sv
   else
      # If not empty, use intersection file
      INPUT_FILE=${sample_id}_OCC_SVs.vcf
   fi

   # Create intersection file
   intersectBed -a $wf_sv -b $occ_fusions -header > ${sample_id}_OCC_SVs.vcf



   # Run AnnotSV with selected input file
   AnnotSV -SVinputFile \$INPUT_FILE \
          -annotationsDir ${params.annotate_dir} \
          -vcf 1 \
          -genomeBuild GRCh38 \
          -outputFile ./annotated_variants.tsv
          -tmpDIR \$TMPDIR
           
   annotsv_fusion_filter.py ./annotated_variants.tsv $occ_fusion_genes_list ${sample_id}_annotSV_fusion_extraction.csv

   # Create the sample-specific copy using the correct path
   cp ./annotated_variants.tsv ${sample_id}_annotated_variants.tsv

   # Verify the copy operation
   if [ ! -f "${sample_id}_annotated_variants.tsv" ]; then
      echo "Error: Copy failed. Debugging information:"
      ls -la
      pwd
      exit 1
   fi

   #clean up temporary files and directory
   rm -rf \$TMPDIR
   rm -f *.tmp.tmp
   """
}


  process annotesvold {
    cpus 4
    memory '50 GB'
    label 'annotesv'
    publishDir "${params.output_path}/structure_variant/annotsv/", mode: "copy", overwrite: true

    input:
        tuple val(sample_id), path(sv_file), path(sv_file_tbi), path(occ_fusions), path(occ_fusion_genes_list)

    output:
        path("${sample_id}_OCC_SVs.vcf")
        tuple val(sample_id), path("annotated_variants.tsv"), emit: annotatedvariantsout
        tuple val(sample_id), path("${sample_id}_annotSV_fusion_extraction.csv"), emit: annotsvfusion
        path("${sample_id}_annotated_variants.tsv")

    script:
    """
    # Prepare input file
    if [[ "${sv_file}" == *.gz && "${sv_file_tbi}" != "NO_INDEX_NEEDED" ]]; then
        gunzip -c ${sv_file} > tmp.vcf
        SV_INPUT="tmp.vcf"
    else
        SV_INPUT="${sv_file}"
    fi
    
    intersectBed -a \$SV_INPUT -b $occ_fusions -header > ${sample_id}_OCC_SVs.vcf
    
    if [ \$(grep -v '^#' ${sample_id}_OCC_SVs.vcf | wc -l) -eq 0 ]; then
        INPUT_FILE=\$SV_INPUT
    else
        INPUT_FILE=${sample_id}_OCC_SVs.vcf
    fi

    # Create working directory
    mkdir -p workdir

    # Run AnnotSV with optimized settings
    AnnotSV \
        -SVinputFile \$INPUT_FILE \
        -annotationsDir ${params.annotate_dir} \
        -hpo "HPO:0100836" \
        -vcf 1 \
        -genomeBuild GRCh38 \
        -outputFile annotated_variants.tsv \
        -candidateGenesFile $occ_fusion_genes_list \
        -outputDir workdir \
        -overwrite 1 \
        -SVminSize 50 \
        -candidateGenesFiltering 1 \
        -annotationMode split
        -metrics us \
        -minTotalNumber 100 \
        -overlap 100

    annotsv_fusion_filter.py ./annotated_variants.tsv $occ_fusion_genes_list ${sample_id}_annotSV_fusion_extraction.csv
    cp ./annotated_variants.tsv ${sample_id}_annotated_variants.tsv

    if [ ! -f "${sample_id}_annotated_variants.tsv" ]; then
        echo "Error: Processing failed. Debugging information:"
        ls -la
        pwd
        exit 1
    fi

    # Cleanup
    rm -rf workdir
    """
}
  

process annotesv_good {
    cpus 4
    memory '4 GB'
    label 'annotesv'
    publishDir "${params.output_path}/structure_variant/annotsv/", mode: "copy", overwrite: true

    input:
        tuple val(sample_id), path(sv_file), path(sv_file_tbi), path(occ_fusions), path(occ_fusion_genes_list)

    output:
        path("${sample_id}_OCC_SVs.vcf")
        tuple val(sample_id), path("annotated_variants.tsv"), emit: annotatedvariantsout
        tuple val(sample_id), path("${sample_id}_annotSV_fusion_extraction.csv"), emit: annotsvfusion
        path("${sample_id}_annotated_variants.tsv")

    script:
    """
    # Prepare input file
    if [[ "${sv_file}" == *.gz && "${sv_file_tbi}" != "NO_INDEX_NEEDED" ]]; then
        gunzip -c ${sv_file} > tmp.vcf
        SV_INPUT="tmp.vcf"
    else
        SV_INPUT="${sv_file}"
    fi
    
    intersectBed -a \$SV_INPUT -b $occ_fusions -header > ${sample_id}_OCC_SVs.vcf
    
    if [ \$(grep -v '^#' ${sample_id}_OCC_SVs.vcf | wc -l) -eq 0 ]; then
        INPUT_FILE=\$SV_INPUT
    else
        INPUT_FILE=${sample_id}_OCC_SVs.vcf
    fi

    # Create directory for split files
    mkdir -p split_files

    # Extract header and variants separately
    grep '^#' \$INPUT_FILE > split_files/header.vcf
    grep -v '^#' \$INPUT_FILE > split_files/variants.vcf

    # Split variants into files of 5 lines each
    cd split_files
    split -l 2 variants.vcf split_variant_
    
    # Process each split file
    for split_file in split_variant_*; do
        # Add header to each split file
        cat header.vcf \$split_file > \${split_file}.vcf
        
        # Run AnnotSV on each split file
        AnnotSV \
            -SVinputFile \${split_file}.vcf \
            -annotationsDir ${params.annotate_dir} \
            #-hpo = "HPO:0100836" \
            -vcf 1 \
            -candidateGenesFile = ${occ_fusion_genes_list} \
            -genomeBuild GRCh38 \
            -candidateGenesFiltering 1 \
            -outputFile \${split_file}_annotated.tsv
    done
    cd ..

    # Merge results - keep header from first file only
    head -n 1 \$(ls split_files/*/split_variant_*.annotated.tsv | head -n 1) > annotated_variants.tsv
    for f in split_files/*/split_variant_*.annotated.tsv; do
        tail -n +2 \$f >> annotated_variants.tsv
    done


    # Run fusion filter on merged results
    annotsv_fusion_filter.py ./annotated_variants.tsv $occ_fusion_genes_list ${sample_id}_annotSV_fusion_extraction.csv
    cp ./annotated_variants.tsv ${sample_id}_annotated_variants.tsv

    # Cleanup
    rm -rf split_files

    if [ ! -f "${sample_id}_annotated_variants.tsv" ]; then
        echo "Error: Processing failed. Debugging information:"
        ls -la
        pwd
        exit 1
    fi
    """
}

  process annotesv {
    debug true
    cpus 4
    memory '4 GB'
    label 'annotesv'
    publishDir "${params.output_path}/structure_variant/annotsv/", mode: "copy", overwrite: true

    input:
        tuple val(sample_id), path(sv_file), path(sv_file_tbi), path(occ_fusions), path(occ_fusion_genes_list)

    output:
        path("${sample_id}_OCC_SVs.vcf")
        tuple val(sample_id), path("annotated_variants.tsv"), emit: annotatedvariantsout
        tuple val(sample_id), path("${sample_id}_annotSV_fusion_extraction.csv"), emit: annotsvfusion
        path("${sample_id}_annotated_variants.tsv")

    script:
    """
    #echo "Processing sample: ${sample_id}"
    #echo "SV file: ${sv_file}"
    #echo "SV index: ${sv_file_tbi}"
    # Prepare input file
    if [[ "${sv_file}" == *.gz && "${sv_file_tbi}" != "NO_INDEX_NEEDED" ]]; then
        gunzip -c ${sv_file} > tmp.vcf
        SV_INPUT="tmp.vcf"
    else
        SV_INPUT="${sv_file}"
    fi
    
    intersectBed -a \$SV_INPUT -b ${occ_fusions} -header > ${sample_id}_OCC_SVs.vcf
    
    if [ \$(grep -v '^#' ${sample_id}_OCC_SVs.vcf | wc -l) -eq 0 ]; then
        INPUT_FILE=\$SV_INPUT
    else
        INPUT_FILE=${sample_id}_OCC_SVs.vcf
    fi

    # Copy candidate genes file to local directory
    cp ${occ_fusion_genes_list} ./candidate_genes.txt

    # Create directory for split files
    mkdir -p split_files

    # Extract header and variants separately
    grep '^#' \$INPUT_FILE > split_files/header.vcf
    grep -v '^#' \$INPUT_FILE > split_files/variants.vcf

    # Calculate optimal split size based on file characteristics
    total_variants=\$(wc -l < split_files/variants.vcf)
    max_line_size=\$(awk '{print length}' split_files/variants.vcf | sort -nr | head -1)
    tcl_limit=2147483647
    
    # Calculate optimal number of variants per split
    # Using 80% of TCL limit to be safe
    safe_limit=\$(( tcl_limit * 8 / 10 ))
    optimal_split=\$(( safe_limit / max_line_size ))
    
    # Ensure split size is at least 1 and no more than 100
    if [ \$optimal_split -lt 1 ]; then
        optimal_split=1
    elif [ \$optimal_split -gt 5 ]; then
        optimal_split=5
    fi
    
    ###echo "Total variants: \$total_variants"
    ###echo "Max line size: \$max_line_size bytes"
    ###echo "Optimal split size: \$optimal_split variants"

    # Split variants using calculated optimal size
    cd split_files
    split -l \$optimal_split variants.vcf split_variant_
    
    # Process each split file
    for split_file in split_variant_*; do
        # Add header to each split file
        cat header.vcf \$split_file > \${split_file}.vcf
        
        # Run AnnotSV on each split file with error handling
        if ! AnnotSV \
            -SVinputFile \${split_file}.vcf \
            -annotationsDir ${params.annotate_dir} \
            -vcf 1 \
            -genomeBuild GRCh38 \
            -candidateGenesFile ../candidate_genes.txt \
            -candidateGenesFiltering 1 \
            -outputFile \${split_file}_annotated.tsv 2> \${split_file}.error; then
            
            # If failed, try with smaller split size
            echo "Failed with \$optimal_split variants, trying with half..."
            half_split=\$(( optimal_split / 2 ))
            split -l \$half_split \$split_file retry_split_
            for retry_file in retry_split_*; do
                cat header.vcf \$retry_file > \${retry_file}.vcf
                AnnotSV \
                    -SVinputFile \${retry_file}.vcf \
                    -annotationsDir ${params.annotate_dir} \
                    -vcf 1 \
                    -genomeBuild GRCh38 \
                    -candidateGenesFile ../candidate_genes.txt \
                    -candidateGenesFiltering 1 \
                    -outputFile \${retry_file}_annotated.tsv || true
            done
        fi
    done
    cd ..

    # Merge results - keep header from first file only
    head -n 1 \$(ls split_files/*/split_variant_*_annotated.tsv split_files/*/retry_split_*_annotated.tsv 2>/dev/null | head -n 1) > annotated_variants.tsv
    for f in split_files/*/split_variant_*_annotated.tsv split_files/*/retry_split_*_annotated.tsv; do
        if [ -f "\$f" ]; then
            tail -n +2 \$f >> annotated_variants.tsv
        fi
    done

    # Run fusion filter on merged results
    annotsv_fusion_filter.py ./annotated_variants.tsv ${occ_fusion_genes_list} ${sample_id}_annotSV_fusion_extraction.csv
    cp ./annotated_variants.tsv ${sample_id}_annotated_variants.tsv

    # Cleanup
    rm -rf split_files candidate_genes.txt

    if [ ! -f "${sample_id}_annotated_variants.tsv" ]; then
     ###   echo "Error: Processing failed. Debugging information:"
        ls -la
        pwd
        exit 1
    fi
    """
}

process annotesv_ori {
    cpus 4
    memory '4 GB'
    label 'annotesv'
    publishDir "${params.output_path}/structure_variant/annotsv/", mode: "copy", overwrite: true

    input:
        tuple val(sample_id), path(sv_file), path(sv_file_tbi), path(occ_fusions), path(occ_fusion_genes_list)

    output:
        path("${sample_id}_OCC_SVs.vcf")
        tuple val(sample_id), path("annotated_variants.tsv"), emit: annotatedvariantsout
        tuple val(sample_id), path("${sample_id}_annotSV_fusion_extraction.csv"), emit: annotsvfusion
        path("${sample_id}_annotated_variants.tsv")

    script:
    """
    # Prepare input file
    if [[ "${sv_file}" == *.gz && "${sv_file_tbi}" != "NO_INDEX_NEEDED" ]]; then
        gunzip -c ${sv_file} > tmp.vcf
        SV_INPUT="tmp.vcf"
    else
        SV_INPUT="${sv_file}"
    fi
    
    intersectBed -a \$SV_INPUT -b ${occ_fusions} -header > ${sample_id}_OCC_SVs.vcf
    
    if [ \$(grep -v '^#' ${sample_id}_OCC_SVs.vcf | wc -l) -eq 0 ]; then
        INPUT_FILE=\$SV_INPUT
    else
        INPUT_FILE=${sample_id}_OCC_SVs.vcf
    fi

    # Copy candidate genes file to local directory
    cp ${occ_fusion_genes_list} ./candidate_genes.txt

    # Create directory for split files
    mkdir -p split_files

    # Extract header and variants separately
    grep '^#' \$INPUT_FILE > split_files/header.vcf
    grep -v '^#' \$INPUT_FILE > split_files/variants.vcf

    # Calculate optimal split size based on file characteristics
    total_variants=\$(wc -l < split_files/variants.vcf)
    max_line_size=\$(awk '{print length}' split_files/variants.vcf | sort -nr | head -1)
    tcl_limit=2147483647
    
    # Calculate optimal number of variants per split
    # Using 80% of TCL limit to be safe
    safe_limit=\$(( tcl_limit * 8 / 10 ))
    optimal_split=\$(( safe_limit / max_line_size ))
    
    # Ensure split size is at least 1 and no more than 100
    if [ \$optimal_split -lt 1 ]; then
        optimal_split=1
    elif [ \$optimal_split -gt 5 ]; then
        optimal_split=5
    fi
    
    echo "Total variants: \$total_variants"
    echo "Max line size: \$max_line_size bytes"
    echo "Optimal split size: \$optimal_split variants"

    # Split variants using calculated optimal size
    cd split_files
    split -l \$optimal_split variants.vcf split_variant_
    
    # Process each split file
    for split_file in split_variant_*; do
        # Add header to each split file
        cat header.vcf \$split_file > \${split_file}.vcf
        
        # Run AnnotSV on each split file with error handling
        if ! AnnotSV \
            -SVinputFile \${split_file}.vcf \
            -annotationsDir ${params.annotate_dir} \
            -vcf 1 \
            -genomeBuild GRCh38 \
            -candidateGenesFile ../candidate_genes.txt \
            -candidateGenesFiltering 1 \
            -outputFile \${split_file}_annotated.tsv 2> \${split_file}.error; then
            
            # If failed, try with smaller split size
            echo "Failed with \$optimal_split variants, trying with half..."
            half_split=\$(( optimal_split / 2 ))
            split -l \$half_split \$split_file retry_split_
            for retry_file in retry_split_*; do
                cat header.vcf \$retry_file > \${retry_file}.vcf
                AnnotSV \
                    -SVinputFile \${retry_file}.vcf \
                    -annotationsDir ${params.annotate_dir} \
                    -vcf 1 \
                    -genomeBuild GRCh38 \
                    -candidateGenesFile ../candidate_genes.txt \
                    -candidateGenesFiltering 1 \
                    -outputFile \${retry_file}_annotated.tsv || true
            done
        fi
    done
    cd ..

    # Merge results - keep header from first file only
    head -n 1 \$(ls split_files/*/split_variant_*_annotated.tsv split_files/*/retry_split_*_annotated.tsv 2>/dev/null | head -n 1) > annotated_variants.tsv
    for f in split_files/*/split_variant_*_annotated.tsv split_files/*/retry_split_*_annotated.tsv; do
        if [ -f "\$f" ]; then
            tail -n +2 \$f >> annotated_variants.tsv
        fi
    done

    # Run fusion filter on merged results
    annotsv_fusion_filter.py ./annotated_variants.tsv ${occ_fusion_genes_list} ${sample_id}_annotSV_fusion_extraction.csv
    cp ./annotated_variants.tsv ${sample_id}_annotated_variants.tsv

    # Cleanup
    rm -rf split_files candidate_genes.txt

    """
}

//visualize AnnotSV result using knotAnnotSV
process knotannotsv {
    cpus 4
    memory '2 GB'
   label 'knotannotsv'
   publishDir "${params.output_path}/structure_variant/annotsv/", mode: "copy", overwrite: true
   
   input:
   tuple val(sample_id), path(annotsv_output), path(knotannov_conf)

    output:
    //tuple val(sample_id), path("annotated_variants.html"), path("annotated_variants.xlsm")
    tuple val(sample_id), path("${sample_id}_annotated_variants.html"), emit:rmdannotsvhtml
    path("${sample_id}_annotated_variants.xlsm")


    script:


    """
 
   knotAnnotSV.pl  --annotSVfile ${annotsv_output} --configFile ${knotannov_conf} --outDir 
   knotAnnotSV2XL.pl  --annotSVfile ${annotsv_output} --configFile ${knotannov_conf} --outDir
   cp annotated_variants.html ${sample_id}_annotated_variants.html
   cp annotated_variants.xlsm ${sample_id}_annotated_variants.xlsm

  # cp "annotated_variants.html" "${params.output_path}/report/${sample_id}_annotsv.html"
    """

}

//visualized Structure variant from AnnotSv using vfc2circos packages


process circosplot {
    cpus 2
    memory '2 GB'
   label 'circos'
   publishDir "${params.output_path}/structure_variant/svannasv/", mode: "copy", overwrite: true
   
   input:
   tuple val(sample_id), path(annotsv_output), path(vcf2circos_json)

   output:
   tuple val(sample_id), path("${sample_id}_vcf2circo.html"), optional: true, emit: circosout

   script:
   """
   # Check if file is empty (excluding header)
   if [ \$(grep -v '^#' ${annotsv_output} | wc -l) -eq 0 ]; then
      echo "Warning: ${annotsv_output} is empty. Skipping vcf2circos plot generation."
      touch ${sample_id}_vcf2circo.html
      exit 0
   else
      vcf2circos -i $annotsv_output -o ${sample_id}_vcf2circo.html -p $vcf2circos_json -a hg38
   fi
   """
}

//# Create an annotated CNV plot from the outputs of Epi2me qDNAseq
// # Uses _segs.vcf which is created by qDNAseq via Epi2melabs on P24
// # Custom Rscript to add annotation labels to CNV plot where gene amplifications or deletions are detected

// Process definition

process annotatecnv {
    cpus 1
    memory '2 GB'
   label 'annotatecnv'
   publishDir "${params.output_path}/copy_number_variation/t20_171_binsize50_ampli10", mode: "copy", overwrite: true

   input:
 //  tuple val(sample_id), path(segsfromepi2me), path(occ_fusions), path(bins_bed), val(threshold)
   tuple val(sample_id), path(segsfromepi2me), path(occ_fusions), path(calls_bed), path(seg_bed), val(threshold)

   output:
   output:
   tuple val(sample_id), path("${sample_id}_calls_fixed.vcf"), emit: callsfixedout
   path("${sample_id}_annotatedcnv.csv")
   path("${sample_id}_annotatedcnv_filter.csv")
   //path("${sample_id}_CNV_plot.pdf")
   path("${sample_id}_annotatedcnv_filter_header.csv"), emit:rmdannotatedcnvfilter
   //path("${sample_id}_CNV_plot.html")
   tuple path("${sample_id}_tumor_copy.txt"), path("${sample_id}_bins_filter.bed")
   //tuple path("${sample_id}_CNV_plot.pdf"), path("${sample_id}_annotatedcnv_filter.csv"), emit:cnvpdfandcsvout
   tuple val(sample_id), path("${sample_id}_cnv_plot_full.pdf"), path("${sample_id}_tumor_copy_number.txt"), path("${sample_id}_annotatedcnv_filter_header.csv"), path("${sample_id}_cnv_chr9.pdf"), path("${sample_id}_cnv_chr7.pdf"), emit: rmdcnvtumornumber
//   source /opt/conda/etc/profile.d/conda.sh
   script:
   """
   #export SSL_CERT_FILE0/etec/ssl/certs/ca-certificates.crt
   source /opt/conda/etc/profile.d/conda.sh
   #conda init
   conda activate annotatecnv_env
   awk 'OFS="\\t" {if (NR > 13) \$1="chr"\$1; print}' $segsfromepi2me  > ${sample_id}_calls_fixed.vcf

   intersectBed -a ${sample_id}_calls_fixed.vcf  -b $occ_fusions -wa -wb | \
   cut -f1,2,5,8,20 | awk '/protein_coding/'| awk -v OFS=";" '\$1=\$1'| awk 'BEGIN { FS=";"; OFS="\t"} {\$1=\$1; print}' | \
   cut -f1,2,3,5,6,8,9,13 > ${sample_id}_annotatedcnv.csv
   
   #cnv_html.R $calls_bed ${sample_id}_annotatedcnv.csv ${sample_id}_CNV_plot.pdf ${sample_id}_CNV_plot.html $sample_id
   CNV_function_new_update.R $calls_bed ${sample_id}_annotatedcnv.csv $seg_bed ${sample_id}_cnv_plot_full.pdf ${sample_id}_cnv_chr9.pdf ${sample_id}_cnv_chr7.pdf $sample_id 

   awk 'BEGIN { OFS="," } {
    gsub(/<[^>]+>/, substr(\$3, 2, length(\$3) - 2), \$3);
    print \$0
   }' ${sample_id}_annotatedcnv.csv > ${sample_id}_annotatedcnv_filter.csv

   awk 'BEGIN {print "Chrom,Start,Type,End,SVLEN,Score,LOG2CNT,Gene"} 1' ${sample_id}_annotatedcnv_filter.csv > ${sample_id}_annotatedcnv_filter_header.csv

   cnv_mapping_occfusion_update.py  $seg_bed $occ_fusions ${sample_id}_tumor_copy.txt ${sample_id}_bins_filter.bed  $threshold

   cnv_mapping_occfusion_update_nofilter.py  $seg_bed  ${sample_id}_tumor_copy_number.txt   $threshold
   """
}


process annotatecnv_old {
    cpus 4
    memory '2 GB'
   label 'annotatecnv'
   publishDir "${params.output_path}/copy_number_variation/", mode: "copy", overwrite: true

   input:
 //  tuple val(sample_id), path(segsfromepi2me), path(occ_fusions), path(bins_bed), val(threshold)
   tuple val(sample_id), path(segsfromepi2me), path(occ_fusions), path(calls_bed), path(seg_bed), val(threshold)

   output:
   tuple val(sample_id), path("${sample_id}_calls_fixed.vcf"), emit: callsfixedout
   path("${sample_id}_annotatedcnv.csv")
   path("${sample_id}_annotatedcnv_filter.csv")
   path("${sample_id}_CNV_plot.pdf")
   path("${sample_id}_annotatedcnv_filter_header.csv"), emit:rmdannotatedcnvfilter
   path("${sample_id}_CNV_plot.html")
   tuple path("${sample_id}_tumor_copy.txt"), path("${sample_id}_bins_filter.bed")
   tuple path("${sample_id}_CNV_plot.pdf"), path("${sample_id}_annotatedcnv_filter.csv"), emit:cnvpdfandcsvout
   tuple val(sample_id), path("${sample_id}_cnv_plot_full.pdf"), path("${sample_id}_tumor_copy_number.txt"), path("${sample_id}_annotatedcnv_filter_header.csv"), path("${sample_id}_cnv_chr9.pdf"), path("${sample_id}_cnv_chr7.pdf"), emit: rmdcnvtumornumber
//   source /opt/conda/etc/profile.d/conda.sh
   script:
   """
   source /opt/conda/etc/profile.d/conda.sh
   conda info --envs
   conda activate annotatecnv_env
   awk 'OFS="\\t" {if (NR > 13) \$1="chr"\$1; print}' $segsfromepi2me  > ${sample_id}_calls_fixed.vcf

   intersectBed -a ${sample_id}_calls_fixed.vcf  -b $occ_fusions -wa -wb | \
   cut -f1,2,5,8,20 | awk '/protein_coding/'| awk -v OFS=";" '\$1=\$1'| awk 'BEGIN { FS=";"; OFS="\t"} {\$1=\$1; print}' | \
   cut -f1,2,3,5,6,8,9,13 > ${sample_id}_annotatedcnv.csv
   
   cnv_html.R $calls_bed ${sample_id}_annotatedcnv.csv ${sample_id}_CNV_plot.pdf ${sample_id}_CNV_plot.html $sample_id
   CNV_function_new_update.R $calls_bed ${sample_id}_annotatedcnv.csv $seg_bed ${sample_id}_cnv_plot_full.pdf ${sample_id}_cnv_chr9.pdf ${sample_id}_cnv_chr7.pdf $sample_id 

   awk 'BEGIN { OFS="," } {
    gsub(/<[^>]+>/, substr(\$3, 2, length(\$3) - 2), \$3);
    print \$0
   }' ${sample_id}_annotatedcnv.csv > ${sample_id}_annotatedcnv_filter.csv

   awk 'BEGIN {print "Chrom,Start,Type,End,SVLEN,Score,LOG2CNT,Gene"} 1' ${sample_id}_annotatedcnv_filter.csv > ${sample_id}_annotatedcnv_filter_header.csv

   cnv_mapping_occfusion_update.py  $calls_bed $occ_fusions ${sample_id}_tumor_copy.txt ${sample_id}_bins_filter.bed  $threshold

   cnv_mapping_occfusion_update_nofilter.py $calls_bed ${sample_id}_tumor_copy_number.txt   $threshold
   """
}


//#######################
//# Call small variants in ROI.bam file with clair3
//###############
//# ROI.bam file created on P24 by cross-referencing merged bam file with list of OCC genes
//# ROI.bam file transfered to office PC for variant calling and annotation
//# clair3 installed from https://github.com/HKU-BAL/Clair3

process clair3 {
    cpus 2
    memory '5 GB'
    label 'clair3'
    publishDir "${params.output_path}/OCC1/$sample_id", mode: "copy", overwrite: true
   
    input:
    tuple val(sample_id), path(occ_bam), path(occ_bam_bai), path(reference_genome), path(reference_genome_bai),  path(refGene), path(hg38_refGeneMrna), path(clinvar), path(clinvarindex),path(hg38_cosmic100),path(hg38_cosmic100index)

    output:
    tuple val(sample_id), path('output_clair3/')
    tuple val(sample_id), path('occ_pileup_snvs_avinput')
    path("${sample_id}_occ_pileup_annotateandfilter.csv")
    tuple val(sample_id), path("${sample_id}_occ_pileup_annotateandfilter.csv"), emit:occpileupannotateandfilterout
    path('occ_merge_snv_avinpt')
    path('occ_merge.hg38_multianno.txt')
    tuple val(sample_id), path("${sample_id}_merge_annotateandfilter.csv"), emit:mergeannotateandfilterout
    tuple path("${sample_id}_occ_pileup_annotateandfilter.csv"), path("${sample_id}_merge_annotateandfilter.csv"), emit:clair3output 


    script:
   
   """ 

   run_clair3.sh \
    --bam_fn=$occ_bam \
    --ref_fn=$reference_genome  \
    --threads=8 \
    --var_pct_full=1 \
    --ref_pct_full=1 \
    --var_pct_phasing=1 \
    --platform="ont" \
    --no_phasing_for_fa \
    --model_path=${params.ref_dir}/r1041_e82_400bps_sup_v420 \
    --output=output_clair3
 
 convert2annovar.pl output_clair3/pileup.vcf.gz \
    --format vcf4 \
	--withfreq \
	--filter pass \
	--fraction 0.1 \
	--includeinfo \
	--outfile occ_pileup_snvs_avinput

   
   table_annovar.pl occ_pileup_snvs_avinput \
         -outfile occ_pileup \
         -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
         -operation g,f,f \
         ${params.humandb_dir} \
         -otherinfo
      
    awk '/exonic/ && /nonsynonymous/ && !/Benign/ && !/Likely_benign/|| /upstream/ || /Func.refGene/ || /splicing/ && !/Benign/ && !/Likely_benign/ || /frameshift/ && !/Benign/ && !/Likely_benign/ || /stopgain/ && !/Benign/ && !/Likely_benign/' occ_pileup.hg38_multianno.txt \
| awk '/exonic/ || /TERT/ || /Func.refGene/'  \
| awk '!/dist=166/' \
| cut -f1-16,26,28,29 > ${sample_id}_occ_pileup_annotateandfilter.csv

convert2annovar.pl \
    output_clair3/merge_output.vcf.gz \
    --format vcf4 \
    --withfreq \
    --filter pass \
    --fraction 0.1 \
    --includeinfo \
    --outfile occ_merge_snv_avinpt

table_annovar.pl occ_merge_snv_avinpt \
    -outfile occ_merge \
    -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
    -operation g,f,f \
    ${params.humandb_dir} \
    -otherinfo

    awk '/exonic/ && /nonsynonymous/ && !/Benign/ && !/Likely_benign/|| /upstream/ || /Func.refGene/ || /splicing/ && !/Benign/ && !/Likely_benign/ || /    frameshift/ && !/Benign/ && !/Likely_benign/ || /stopgain/ && !/Benign/ && !/Likely_benign/' \
        occ_merge.hg38_multianno.txt \
    | awk '/exonic/ || /TERT/ || /Func.refGene/'  \
    | awk '!/dist=166/' \
    | cut -f1-16,26,28,29 \
    > ${sample_id}_merge_annotateandfilter.csv 

    """
   }


//#################################
//#### ClairS-TO
//#################################
//# ClairS-TO is a recent development to specifically call somatic variants in Tumor-only samples
//# It's run separate from Clair3
// # installed from https://github.com/HKU-BAL/ClairS-TO via micromamba

process clairs_to {
    cpus 2
    memory '2 GB'
    label 'clairsto'
    publishDir "${params.output_path}/OCC1/$sample_id", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(occ_bam), path(occ_bam_bai), path(reference_genome), path(reference_genome_bai),  path(refGene), path(hg38_refGeneMrna), path(clinvar), path(clinvarindex),path(hg38_cosmic100),path(hg38_cosmic100index), path(occ_snv_screening)
    
    output:
    tuple val(sample_id), path('clairsto_output/')
    tuple val(sample_id), path('clairS_To_snv_avinput')
    path('ClairS_TO_snv.hg38_multianno.txt')
    tuple val(sample_id), path("${sample_id}_annotateandfilter_clairsto.csv"), emit:annotateandfilter_clairstoout
    path("${sample_id}_merge_snv_indel_claisto.vcf.gz")


    script:

   """

   run_clairs_to \
    --tumor_bam_fn=$occ_bam \
    --ref_fn=$reference_genome  \
    --threads=8 \
    --platform="ont_r10_dorado_4khz" \
    --output_dir=clairsto_output \
    --bed_fn=$occ_snv_screening \


    bcftools merge --force-samples clairsto_output/snv.vcf.gz clairsto_output/indel.vcf.gz -o ${sample_id}_merge_snv_indel_claisto.vcf.gz

    convert2annovar.pl ${sample_id}_merge_snv_indel_claisto.vcf.gz \
   --format vcf4 \
   --filter pass \
   --includeinfo \
   --outfile  clairS_To_snv_avinput


  table_annovar.pl clairS_To_snv_avinput \
   -outfile ClairS_TO_snv \
   -buildver hg38 -protocol refGene,clinvar_20240611,cosmic100coding2024\
   -operation g,f,f \
    ${params.humandb_dir} \
   -otherinfo  

   awk '/exonic/ && /nonsynonymous/ && !/Benign/ || /upstream/ || /Func.refGene/' \
   ClairS_TO_snv.hg38_multianno.txt \
   | awk '/exonic/ || /TERT/ || /Func.refGene/'  \
  | awk '!/dist=166/' \
  | cut -f1-16,25,26  > ${sample_id}_annotateandfilter_clairsto.csv


    """
   }


// Simple R-script to merge annotation calls

process merge_annotation {
    cpus 2
    memory '2 GB'
    label 'merge_annotation'
    publishDir "${params.output_path}/merge_annot_clair3andclairsto_update/", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(mergeannotateandfilterout), path(occpileupannotateandfilterout), path(annotateandfilter_clairstoout), path(occ_genes)
    
    output:
    tuple val(sample_id), path("${sample_id}_merge_annotation_filter_snvs_allcall.csv"), emit: occmergeout
    tuple val(sample_id), path("${sample_id}_merge_annotation_filter_snvs_allcall_filter.csv")

    """

    merge_annotations_prospective24june25.R $mergeannotateandfilterout $occpileupannotateandfilterout $annotateandfilter_clairstoout ${sample_id}_merge_annotation_filter_snvs_allcall.csv $occ_genes ${sample_id}_merge_annotation_filter_snvs_allcall_filter.csv 

    """
}
//# Uses igv_tools to extract 3 regions with recurring mutations.
process igv_tools {
    cpus 2
    memory '2 GB'
    label 'epic'


    publishDir "${params.output_path}/terp", mode: "copy", overwrite: true


    input:
    tuple val(sample_id),path(occ_bam), path(occ_bam_bai),  path(tertp_variants),path(ncbirefseq), path(reference_genome), path(reference_genome_bai)

    output:
    tuple val(sample_id), file("${sample_id}_tertp_id1.html"), emit: tertp_out_igv

    script:

    """
    export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
    create_report $tertp_variants \
        --fasta $reference_genome\
        --flanking 1000 \
        --tracks $tertp_variants $occ_bam  $ncbirefseq \
        --output ${sample_id}_tertp_id1.html





    """
}

    process cramino_report {
        cpus 2
        memory '2 GB'
        label 'epic'
        publishDir "${params.output_path}/cramino", mode: "copy", overwrite: true


    input:
    tuple val(sample_id),path(merge_bam), path(merge_bam_bai), path(reference_genome), path(reference_genome_bai)

    output:
    tuple val(sample_id), file("${sample_id}_cramino_statistics.txt"), emit:craminostatout
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda env list
    conda activate annotatecnv_env

    cramino $merge_bam --reference $reference_genome > ${sample_id}_cramino_statistics.txt
   """

    }

    process ace_tmc {
    label 'ace_tmc'
    cpus 1
    memory '2 GB'
    publishDir "${params.output_path}/ace/", mode: "copy", overwrite: true
    
    input:
    tuple val(sample_id), path(bam_dir)
    
    output:
    tuple val(sample_id), path("${sample_id}_ace_results"), emit: aceresults
    tuple val(sample_id), env(threshold_value), emit: threshold_value
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Create temporary directory for single sample
    mkdir -p temp_input_dir
    cp "${bam_dir}/${sample_id}_copyNumbersCalled.rds" temp_input_dir/
    
    # Check if we're in a container and use appropriate conda setup
    if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
        source /opt/conda/etc/profile.d/conda.sh
        conda activate ace_env
    else
        source activate ace_env
    fi
    
    # Debug info
    echo "Processing sample: ${sample_id}"
    echo "Input directory: temp_input_dir"
    ls -l temp_input_dir
    
    # Create output directory
    mkdir -p ${sample_id}_ace_results
    
    # Run ACE TMC analysis
    ace_tmc.R "temp_input_dir" "${sample_id}_ace_results" "${sample_id}"
    
    # Read and export the threshold value
    threshold_value=\$(cat "${sample_id}_ace_results/threshold_value.txt")
    echo "Threshold value for ${sample_id}: \$threshold_value"
    
    # Cleanup
    rm -rf temp_input_dir
    """
}

    process plot_genomic_regions {
        cpus 2
        memory '2 GB'
       // label 'plot_genomic_regions'
        publishDir "${params.output_path}/tertp", mode: "copy", overwrite: true

        input:
        tuple path(gviz_data), val(sample_id), path(occ_bam), path(occ_bam_bai), path(cytoband_file)

        output:
        tuple val(sample_id), path("${sample_id}_egfr_coverage.pdf"), path("${sample_id}_idh1_coverage.pdf"), path("${sample_id}_tertp_coverage.pdf"), emit:plot_genomic_regions_out
        
        script:
        """
        #!/bin/bash
        ls -l
       echo "Rscript path: \$(which Rscript)"
       echo "GViz file: ${gviz_data}"
        eval \"\$(micromamba shell hook --shell bash)\"
        micromamba activate gviz_env
        plot_genomic_regions.R \
            ${gviz_data} \
            ${sample_id} \
            ${occ_bam} \
            ${sample_id}_egfr_coverage.pdf \
            ${sample_id}_idh1_coverage.pdf \
            ${sample_id}_tertp_coverage.pdf \
            ${cytoband_file}

        """
    }

    

    process markdown_report_24april {
    cpus 2
    memory '2 GB'
    publishDir "${params.output_path}/report", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), 
          val(craminoreport),
          val(nanodx), 
          val(dictionaire), 
          val(logo), 
          val(cnv_plot), 
          val(tumor_number), 
          val(annotatecnv),
          val(cnv_chr9),
          val(cnv_chr7), 
          val(mgmt_results),
          val(merge_results),
          val(annotSV_fusion), 
          val(terphtml),
          val(svannahtml), 
          val(annotsvhtml),
          val(egfr_coverage),
          val(idh1_coverage),
          val(tertp_coverage)

    output:
    file("${sample_id}_markdown_pipeline_report.pdf")


    script:
    """
    Rscript -e "rmarkdown::render('/data/pipeline/nextflow/bin/nextflow_markdown_pipeline3.Rmd',output_file=commandArgs(trailingOnly=TRUE)[20])" "${sample_id}" "${craminoreport}" "${nanodx}" "${dictionaire}" "${logo}" "${cnv_plot}" "${tumor_number}" "${annotatecnv}" "${cnv_chr9}" "${cnv_chr7}" "${mgmt_results}" "${merge_results}" "${annotSV_fusion}"  "${terphtml}" "${svannahtml}" "${annotsvhtml}" "${egfr_coverage}" "${idh1_coverage}" "${tertp_coverage}" "\${PWD}/${sample_id}_markdown_pipeline_report.pdf"

    """
}
    
    process markdown_report {
    cpus 4
    memory '2 GB'
    publishDir "${params.output_path}/report", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), 
          val(craminoreport),
          val(nanodx), 
          val(dictionaire), 
          val(logo), 
          val(cnv_plot), 
          val(tumor_number), 
          val(annotatecnv),
          val(cnv_chr9),
          val(cnv_chr7), 
          val(mgmt_results),
          val(merge_results),
          val(annotSV_fusion), 
          val(terphtml),
          val(svannahtml), 
          val(annotsvhtml),
          val(egfr_coverage),
          val(idh1_coverage),
          val(tertp_coverage)

    output:
    file("${sample_id}_markdown_pipeline_report.pdf")

    script:
    """
    mkdir -p markdown_tmp_${sample_id}
    Rscript -e "rmarkdown::render('/data/pipeline/nextflow/bin/nextflow_markdown_pipeline3.Rmd', intermediates_dir='markdown_tmp_${sample_id}',knit_root_dir='.', clean=TRUE, envir=new.env(),output_file=commandArgs(trailingOnly=TRUE)[20])" \
"${sample_id}" \
"${craminoreport}" \
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
"${annotSV_fusion}" \
"${terphtml}" \
"${svannahtml}" \
"${annotsvhtml}" \
"${egfr_coverage}" \
"${idh1_coverage}" \
"${tertp_coverage}" \
"\${PWD}/${sample_id}_markdown_pipeline_report.pdf"
    """
}

    process markdown_report_old {
    cpus 4
    memory '2 GB'
    publishDir "${params.output_path}/report", mode: "copy", overwrite: true
    label 'markdown'

    input:
    tuple val(sample_id), 
          val(craminoreport),
          val(nanodx), 
          val(dictionaire), 
          val(logo), 
          val(cnv_plot), 
          val(tumor_number), 
          val(annotatecnv),
          val(cnv_chr9),
          val(cnv_chr7), 
          val(mgmt_results),
          val(merge_results),
          val(annotSV_fusion), 
          val(terphtml),
          val(svannahtml), 
          val(annotsvhtml),
          val(egfr_coverage),
          val(idh1_coverage),
          val(tertp_coverage)

    output:
    tuple val(sample_id), path("${sample_id}_markdown_pipeline_report.pdf"), emit: markdown_report

    script:
    """
    # Create header.tex file
    cat << 'EOT' > header.tex
    \\usepackage[utf8]{inputenc}
    \\usepackage{fontspec}
    \\setmainfont{Arial}
    \\usepackage{xcolor}
    \\usepackage{booktabs}
    \\usepackage{longtable}
    \\usepackage{array}
    \\usepackage{multirow}
    \\usepackage{float}
    \\usepackage{makecell}
    \\usepackage{graphicx}
    \\usepackage{caption}
    \\usepackage{placeins}
    EOT

    # Generate report with proper document structure
    Rscript -e '
    rmarkdown::render(
        "${workflow.projectDir}/bin/nextflow_markdown_pipeline3.Rmd",
        output_format = rmarkdown::pdf_document(
            latex_engine = "xelatex",
            includes = list(in_header = "header.tex"),
            keep_tex = TRUE
        ),
       # output_file = "${sample_id}_markdown_pipeline_report.pdf",
        params = list(
            sample_id = "${sample_id}",
            cramino_stat = "${craminoreport}",
            nanodx = "${nanodx}",
            dictionary_file = "${dictionaire}",
            logo_file = "${logo}",
            copy_number_plot_file = "${cnv_plot}",
            tumor_copy_number_file = "${tumor_number}",
            cnv_filter_file = "${annotatecnv}",
            cnv_chr9 = "${cnv_chr9}",
            cnv_chr7 = "${cnv_chr7}",
            mgmt_results_file = "${mgmt_results}",
            snv_results_file = "${merge_results}",
            structure_variant_file = "${annotSV_fusion}",
            terp_html = "${terphtml}",
            svanna_html = "${svannahtml}",
            annotsv_html = "${annotsvhtml}",
            egfr_plot_file = "${egfr_coverage}",
            idh1_plot_file = "${idh1_coverage}",
            tertp_plot_file = "${tertp_coverage}",
            output_pdf_file = "${sample_id}_markdown_pipeline_report.pdf"
        )
    )'
    """
    }





workflow {
    boosts_merge_annotation = []
    start_time = new Date()

    // Read sample_ids.txt and build a Map: [ sampleID : threshold_value, ... ]
    // def sample_thresholds = params.sample_id_file.readLines()
    //     .collectEntries { line ->
    //         def tokens = line.tokenize("\t").view()
    //         if (tokens.size() < 2) {
    //             throw new IllegalArgumentException("Invalid line format in sample_id_file: ${line}")
    //         }
    //         [(tokens[0]) : tokens[1]]
    // }
    // println "Sample thresholds: ${sample_thresholds}"

    def sample_thresholds = loadSampleThresholds()
    
    // Read BAM files and validate presence
    merge_bam_files = Channel.fromPath("${merge_bam_folder}/*.bam").map { bam -> 
        def sample_id = bam.getBaseName().split("\\.")[0]
        def bai = file(bam.toString() + '.bai')
        if (!bai.exists()) {
            throw new FileNotFoundException("Missing index file for ${bam}")
        }
        return tuple(sample_id, bam, bai)
    }
    .filter { it[0] in sample_thresholds }
    
    boosts_cramino = merge_bam_files.map { sample_id, bam, bai -> 
        return tuple(sample_id, bam, bai, params.reference_genome, params.reference_genome_bai)
    }
    
    // Load VCF files and check existence
    

    sv_files = Channel.fromPath("${sv_folder}/*.wf_sv.{vcf,vcf.gz}")
    boosts_annotsv_channel = sv_files.map { vcf -> 
        def sample_id = vcf.getBaseName().split("\\.")[0]
        def tbi = file(vcf.toString() + '.tbi')
        return tuple(sample_id, vcf, tbi, params.occ_fusions, params.occ_fusion_genes_list)
    }
    .filter { it[0] in sample_thresholds }
    
    boosts_svanna_channel = sv_files.map { vcf -> 
        def sample_id = vcf.getBaseName().split("\\.")[0]
        def tbi = file(vcf.toString() + '.tbi')
        return tuple(sample_id, vcf, tbi, params.occ_fusions)
    }
    .filter { it[0] in sample_thresholds }
    
    // Read and validate BED files
    bedmethyl_files = Channel.fromPath("${bedmethyl_folder}/*.wf_mods.bedmethyl.gz")
    boosts_epic_channel = bedmethyl_files.map { gz -> 
        def sample_id = gz.getBaseName().split("\\.")[0]
        return tuple(sample_id, gz, params.epicsites, params.mgmt_cpg_island_hg38)
    }
    .filter { it[0] in sample_thresholds }
    
    // Read tumor content ratio files
    segsfromepi2me_files = Channel.fromPath("${segsfromepi2me_folder}/*.{vcf,vcf.gz}").map { vcf ->
    def sample_id = vcf.getBaseName().split("_")[0]
    if (!vcf.exists()) {
        throw new FileNotFoundException("Missing segmentation VCF for ${sample_id}: ${vcf}")
    }
    return tuple(sample_id, vcf)  
}
    
   calls_bed_files = Channel.fromPath("${segsfromepi2me_folder}/*_bins.bed").filter { bed ->
        !bed.name.contains("raw")
   }
   .map {bed ->
        def sample_id = bed.getBaseName().split("_")[0]
        if (!bed.exists()) {
            throw new FileNotFoundException("Missing calls file for ${sample_id}: ${bed}")
    }
    return tuple(sample_id, bed)  
}

    // Channel for *_segs.bed files with existence check
    seg_bed_files = Channel.fromPath("${segsfromepi2me_folder}/*_segs.bed").map { bed ->
    def sample_id = bed.getBaseName().split("_")[0]
    if (!bed.exists()) {
        throw new FileNotFoundException("Missing segmentation file for ${sample_id}: ${bed}")
    }
    return tuple(sample_id, bed)  
}

   // Create the ace_tmc input channel with sample filtering
    ace_input = Channel
        .fromPath("${cnv_rds}/*_*.rds")
        .map { rds_file -> 
            def sample_id = rds_file.name.toString().split("_")[0]
            if (sample_thresholds.containsKey(sample_id)) {
                println "Found matching RDS file for sample ${sample_id}: ${rds_file}"
                tuple(sample_id, rds_file.parent)
            } else {
                //println "Skipping non-listed sample: ${sample_id}"
                null
            }
        }
        .filter { it != null }
        .filter { it[0] in sample_thresholds }

    // Run ACE analysis
    ace_tmc(ace_input)

    // Create channel for segsfromepi2me with thresholds from ace_tmc
    boosts_segsfromepi2me_channel = segsfromepi2me_files
        .filter { it[0] in sample_thresholds }
        .join(calls_bed_files)
        .join(seg_bed_files)
        .join(ace_tmc.out.threshold_value)  // Join with ace_tmc threshold output
        .map { sample_id, vcf, calls_bed, seg_bed, threshold_value ->
            println "ACE threshold for sample_id ${sample_id}: ${threshold_value}"
            tuple(
                sample_id, 
                vcf, 
                params.occ_fusions, 
                calls_bed, 
                seg_bed, 
                threshold_value
            )
        }
    
    // Read BAM files for variant calling
    bam_files = Channel.fromPath("${bams_folder}/*.bam").map { bam -> 
        def sample_id = bam.getBaseName().split("\\.")[0]
        def bai = file(bam.toString() + '.bai')
        if (!bai.exists()) {
            throw new FileNotFoundException("Missing index file for ${bam}")
        }
        return tuple(sample_id, bam, bai)
    }
    .filter { it[0] in sample_thresholds }
    
    boosts_clair3_channel = bam_files.map { sample_id, bam, bai -> 
        return tuple(sample_id, bam, bai, params.reference_genome, params.reference_genome_bai, params.refgene,
                     params.hg38_refgenemrna, params.clinvar, params.clinvarindex, params.hg38_cosmic100, params.hg38_cosmic100index)
    }
    
    boosts_clairSTo_channel = bam_files.map { sample_id, bam, bai -> 
        return tuple(sample_id, bam, bai, params.reference_genome, params.reference_genome_bai, params.refgene,
                     params.hg38_refgenemrna, params.clinvar, params.clinvarindex, params.hg38_cosmic100, 
                     params.hg38_cosmic100index, params.occ_snv_screening)
    }
    
    boosts_igv_channel = bam_files.map { sample_id, bam, bai -> 
        return tuple(sample_id, bam, bai, params.tertp_variants, params.ncbirefseq, params.reference_genome, params.reference_genome_bai)
    }
    
    boosts_plot_genomic_regions_channel = bam_files.map { sample_id, bam, bai -> 
                return tuple(params.gviz_data, sample_id, bam, bai, params.cytoband_file)
    }
    // Conditional logic to run the appropriate process based on the mode

    if (params.run_mode == 'methylation') {
        // Run only process_a
        extract_epic(boosts_epic_channel)
        MGMT_output = extract_epic.out.MGMTheaderout
        MGMT_sturgeon = extract_epic.out.sturgeonbedinput.map{sample_id, sturgeoninput -> [sample_id, sturgeoninput, params.sturgeon_model]}
        //sturgeon(MGMT_sturgeon)
        mgmt_promoter(MGMT_output)
        mgmt_promoter_out=mgmt_promoter.out.mgmtresultsout //.view()
        mgmt_nanodx = extract_epic.out.epicselectnanodxinput.map{sample_id, epicselectnanodxinput -> [sample_id, epicselectnanodxinput, params.hg19_450model]}
        nanodx(mgmt_nanodx).view()

        //def nanodx_out = nanodx(mgmt_nanodx).out.nanodx450out

        nanodx_out=nanodx.out.nanodx450out.map {sample_id, bed_file -> tuple(sample_id, bed_file, file(params.nanodx_450model), file(params.nn_model))}




// Pass to the classifier process
        run_nn_classifier(nanodx_out)
    }
        
    //    nanodx_beds_ch = nanodx.out.nanodx450out.map{ sample_id, bed_file ->
    //        tuple(sample_id, bed_file)}
    //    nanodx_input = sample_thresholds
    //        .map {sid -> tuple(sid)}.view()
    //        .join(nanodx_beds_ch, by: 0)
    //        .map{sample_id, bed_file -> tuple(sample_id, bed_file, params.nanodx_450model,params.snakefile_nanodx, params.nn_mode)}
        //nanodx_out=nanodx.out.nanodx450out.map{sample_id, nanodx450out -> [sample_id, nanodx450out, params.nanodx_450model,params.snakefile_nanodx, params.nn_model]}
    //    run_nn_classifier(nanodx_out)
    //}

    else if (params.run_mode == 'annotsv') {
        // Run only process_b
        svannasv(boosts_svanna_channel)
        svannasv_out = svannasv.out.occsvannaannotationannotationvcf.map{sample_id,svannavcfoutput -> [sample_id, svannavcfoutput, params.vcf2circos_json]}
        circosplot(svannasv_out)
        circosplot_out=circosplot.out.circosout //.view()
        svannaoutfusion_events= svannasv.out.occsvannavcfout.map{sample_id, occsvannavcfout -> [sample_id, occsvannavcfout, params.genecode_bed, params.occ_fusion_genes_list]}
        svannasv_fusion_events(svannaoutfusion_events)
     ///   annotesv(boosts_annotsv_channel)
     ///   annotsv_output = annotesv.out.annotatedvariantsout.map{sample_id,annotated_variants -> [sample_id, annotated_variants, params.knotannotsv_conf]}
    ///    knotannotsv(annotsv_output)
    }

    // run cnv pipeline 
    else if (params.run_mode == 'cnv') {

        annotatecnv(boosts_segsfromepi2me_channel)
    //    cramino_report(boosts_cramino)
//        plot_genomic_regions(boosts_plot_genomic_regions_channel)
     //   annotatecnv_report=annotatecnv.out.cnvpdfandcsvout //.view()
     //   boosts_cramino.out.craminostat.view()

    }


    else if(params.run_mode == 'occ') {

        plot_genomic_regions(boosts_plot_genomic_regions_channel)
        igv_tools(boosts_igv_channel)
        clair3(boosts_clair3_channel)
        clair3_out=clair3.out.clair3output
        clairs_to(boosts_clairSTo_channel)
        clairs_to_out=clairs_to.out.annotateandfilter_clairstoout
        combine_file = clair3_out.combine(clairs_to_out).map{def occ_pileup_annotateandfilter, def merge_annotateandfilter, def sample_id, def annotateandfilter_clairsto -> [sample_id, merge_annotateandfilter, occ_pileup_annotateandfilter, annotateandfilter_clairsto, params.occ_genes]}
        merge_annotation(combine_file)
        merge_annotation_out=merge_annotation.out.occmergeout //.view()
    }

    else if (params.run_mode == 'terp') {
        igv_tools(boosts_igv_channel)

    }

    else if (params.run_mode == 'stat') {
        cramino_report(boosts_cramino)
    }

    else {
        //run MGMT and classifier using sturgeon and nanodx
       extract_epic(boosts_epic_channel)
        MGMT_output = extract_epic.out.MGMTheaderout
        MGMT_sturgeon = extract_epic.out.sturgeonbedinput.map { sample_id, sturgeoninput -> [ sample_id, sturgeoninput, params.sturgeon_model ]
    }
       // sturgeon(MGMT_sturgeon)
        mgmt_promoter(MGMT_output)
        mgmt_promoter_out = mgmt_promoter.out.mgmtresultsout //.view()
        mgmt_nanodx = extract_epic.out.epicselectnanodxinput.map { sample_id, epicselectnanodxinput -> [ sample_id, epicselectnanodxinput, params.hg19_450model ]
    }
        nanodx(mgmt_nanodx).view()
        nanodx_out = nanodx.out.nanodx450out.map { sample_id, nanodx450out -> [ sample_id, nanodx450out, params.nanodx_450model, params.nn_model ]
    }
        run_nn_classifier(nanodx_out)
        rmd_nanodx_out = run_nn_classifier.out.rmdnanodx //.view()

        svannasv(boosts_svanna_channel)
        rmd_svanna_html = svannasv.out.rmdsvannahtml
        svannasv_out = svannasv.out.occsvannaannotationannotationvcf.map { sample_id, svannavcfoutput -> [ sample_id, svannavcfoutput, params.vcf2circos_json ]
    }
        circosplot(svannasv_out)
        circosplot_out = circosplot.out.circosout //.view()
        svannaoutfusion_events= svannasv.out.occsvannavcfout.map{sample_id, occsvannavcfout -> [sample_id, occsvannavcfout, params.genecode_bed, params.occ_fusion_genes_list]}
        svannasv_fusion_events(svannaoutfusion_events)

        ///annotesv(boosts_annotsv_channel)
        ///annotsv_output = annotesv.out.annotatedvariantsout.map { sample_id, annotated_variants -> [ sample_id, annotated_variants, params.knotannotsv_conf ]
    ///}
        ///rmd_annotsvfusion = annotesv.out.annotsvfusion
        ///knotannotsv(annotsv_output)
        ///rmd_knotannotsv_html = knotannotsv.out.rmdannotsvhtml

        annotatecnv(boosts_segsfromepi2me_channel)
//        def rmd_annotatecnv_report = annotatecnv.out.rmdcnvtumornumber 

        clair3(boosts_clair3_channel)
        clair3_out = clair3.out.clair3output
        clairs_to(boosts_clairSTo_channel)
        clairs_to_out = clairs_to.out.annotateandfilter_clairstoout
        combine_file = clair3_out.combine(clairs_to_out).map { occ_pileup_annotateandfilter, merge_annotateandfilter, sample_id, annotateandfilter_clairsto -> [ sample_id, merge_annotateandfilter, occ_pileup_annotateandfilter, annotateandfilter_clairsto, params.occ_genes ]
    }
        merge_annotation(combine_file)
        igv_tools(boosts_igv_channel)
        ///cramino_report(boosts_cramino)
        //occmergefile_out = merge_annotation.out.occmergefile
        // igv_tools(boosts_igv_channel)
        // cramino_report(boosts_cramino)
        plot_genomic_regions(boosts_plot_genomic_regions_channel)


        
      

    }
    
    }   

    workflow.onComplete {
    try {
        def end_time = new Date()
        def duration = (end_time.time - start_time.time) / 1000
        println "The script ran successfully."
        println "Total running time: ${duration} seconds"
    } catch (Exception e) {
        println "Error in onComplete block: ${e.message}"
        e.printStackTrace()
    }
}

// Function to load sample thresholds
def loadSampleThresholds() {
    def sampleFile
    
    // Check if sample_id_file is provided via command line
    if (params.sample_id_file instanceof String || params.sample_id_file instanceof GString) {
        sampleFile = file(params.sample_id_file)
    } else {
        sampleFile = params.sample_id_file
    }

    // Validate file existence
    if (!sampleFile.exists()) {
        error """
            Sample ID file not found: ${sampleFile}
            Please either:
            1. Provide the file path using --sample_id_file parameter
            2. Place the file at the default location: ${params.path}/testdata/sample_ids.txt
            """
    }

    def thresholds = [:]
    sampleFile.eachLine { line ->
        if (line.trim()) {  // Skip empty lines
            def parts = line.trim().split(/\s+/)  // Split on whitespace
            if (parts.size() >= 2) {
                thresholds[parts[0]] = parts[1]
                println "Loaded sample: ${parts[0]} with threshold: ${parts[1]}"
            } else {
                println "Warning: Skipping invalid line in sample file: ${line}"
            }
        }
    }

    if (thresholds.isEmpty()) {
        error "No valid sample entries found in ${sampleFile}"
    }

    println "Loaded ${thresholds.size()} samples from ${sampleFile}"
    return thresholds
}

// Update help message
def helpMessage() {
    log.info"""
    Usage:
        nextflow run pipeline.nf [--sample_id_file /path/to/samples.txt] [options]

    Arguments:
        --sample_id_file    Path to tab-delimited file containing sample IDs and thresholds
                           Format: <sample_id> <threshold>
                           Example: T20-170 0.63
                           Default: ${params.path}/testdata/sample_ids.txt

        --run_mode         Analysis mode to run (default: all)
                          Options: methylation, annotsv, cnv, occ, terp, rmd, stat, all

    Note: The sample_id_file can be specified either in the script or via command line flag.
          Command line flag takes precedence if provided.
    """
}

// Add parameter validation
if (params.help) {
    helpMessage()
    exit 0
}

