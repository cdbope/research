#!/usr/bin/env nextflow

/*
 * Nanopore Neural Glioblastoma Classification Pipeline (v2)
 * Based on Drexler et al. (Nature Medicine 2024)
 *
 * This pipeline classifies Nanopore methylation data as High-Neural or Low-Neural
 * using logistic regression with all 1289 CpG probes (original method from run_example.ipynb)
 */

nextflow.enable.dsl=2

// Default parameters
params.input_dir = null
params.sample_ids = null
params.output_dir = "./results"
params.min_coverage = 1
params.genome = "hg38"
params.script_dir = "${projectDir}/../code"
params.probe_coords = "${projectDir}/../code/probe_coordinates_450k_hg38.csv"

// Help message
def helpMessage() {
    log.info """
    =========================================
    Nanopore Neural Glioblastoma Classification Pipeline (v2)
    =========================================

    Usage:
        nextflow run main.nf --input_dir <DIR> --sample_ids <FILE> --output_dir <DIR>

    Required arguments:
        --input_dir       Directory containing *.wf_mods.bedmethyl.gz files
        --sample_ids      Text file with sample IDs (one per line)

    Optional arguments:
        --output_dir      Output directory (default: ./results)
        --min_coverage    Minimum read coverage per CpG (default: 5)
        --genome          Reference genome: hg19 or hg38 (default: hg38)
        --script_dir      Path to Drexler code directory

    Example:
        nextflow run main.nf \\
            --input_dir /path/to/bedmethyl_files \\
            --sample_ids sampleids.txt \\
            --output_dir results \\
            --genome hg38
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input_dir) {
    error "Please specify --input_dir"
}
if (!params.sample_ids) {
    error "Please specify --sample_ids"
}

// Log parameters
log.info """
=========================================
Nanopore Neural Glioblastoma Classification (v2)
=========================================
Input directory  : ${params.input_dir}
Sample IDs file  : ${params.sample_ids}
Output directory : ${params.output_dir}
Min coverage     : ${params.min_coverage}
Genome           : ${params.genome}
Script directory : ${params.script_dir}
Probe coords     : ${params.probe_coords}
=========================================
"""

/*
 * Read sample IDs from file and create channel
 */
process GET_SAMPLES {
    input:
    path sample_ids_file

    output:
    stdout

    script:
    """
    cat ${sample_ids_file}
    """
}

/*
 * Run classification for each sample
 */
process CLASSIFY_SAMPLE {
    tag "${sample_id}"
    publishDir "${params.output_dir}/individual", mode: 'copy'

    input:
    tuple val(sample_id), path(bedmethyl_file)
    path probe_coords

    output:
    tuple val(sample_id), path("${sample_id}_prediction.csv"), emit: predictions
    tuple val(sample_id), path("${sample_id}_betas.csv"), emit: betas

    script:
    """
    python ${params.script_dir}/run_nanopore_classification_v2.py \\
        ${bedmethyl_file} \\
        --sample_name ${sample_id} \\
        --output_dir . \\
        --min_coverage ${params.min_coverage} \\
        --genome ${params.genome} \\
        --probe_coords ${probe_coords}
    """
}

/*
 * Aggregate all predictions into a summary table
 */
process AGGREGATE_RESULTS {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path prediction_files

    output:
    path "neural_classification_summary.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    import os

    # Read all prediction files
    all_results = []
    for pred_file in glob.glob("*_prediction.csv"):
        df = pd.read_csv(pred_file, index_col=0)
        for sample_id, row in df.iterrows():
            all_results.append({
                'sample_id': sample_id,
                'Prediction': int(row['Prediction']),
                'Prediction_score': round(row['Prediction_score'], 4),
                'Neural_Classification': row['Neural_Classification']
            })

    # Create summary DataFrame
    summary_df = pd.DataFrame(all_results)
    summary_df = summary_df.sort_values('sample_id')

    # Save as TSV
    summary_df.to_csv('neural_classification_summary.tsv', sep='\\t', index=False)

    print(f"Aggregated {len(summary_df)} samples")
    """
}

/*
 * Main workflow
 */
workflow {
    // Pre-computed probe coordinates file (shared across all processes)
    probe_coords_ch = Channel.value(file(params.probe_coords))

    // Read sample IDs from file
    sample_ids_ch = Channel
        .fromPath(params.sample_ids)
        .splitText()
        .map { it.trim() }
        .filter { it }  // Remove empty lines

    // Create channel with sample ID and corresponding bedmethyl file
    samples_ch = sample_ids_ch
        .map { sample_id ->
            def bedmethyl = file("${params.input_dir}/${sample_id}.wf_mods.bedmethyl.gz")
            if (!bedmethyl.exists()) {
                log.warn "WARNING: File not found for sample ${sample_id}: ${bedmethyl}"
                return null
            }
            return tuple(sample_id, bedmethyl)
        }
        .filter { it != null }

    // Run classification with pre-computed probe coordinates
    CLASSIFY_SAMPLE(samples_ch, probe_coords_ch)

    // Collect all predictions and aggregate
    all_predictions = CLASSIFY_SAMPLE.out.predictions
        .map { sample_id, pred_file -> pred_file }
        .collect()

    AGGREGATE_RESULTS(all_predictions)
}

workflow.onComplete {
    log.info """
    =========================================
    Pipeline completed!
    =========================================
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration  : ${workflow.duration}
    Output    : ${params.output_dir}
    Summary   : ${params.output_dir}/neural_classification_summary.tsv
    =========================================
    """
}
