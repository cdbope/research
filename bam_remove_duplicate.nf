#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters
params {
    merge_bam_folder = "/home/prom/africa/angola/200_GBMs_merge_bam/repeat_double_cell/"  // Directory containing BAM files
    output_dir = "/home/prom/africa/angola/200_GBMs_merge_bam/repeat_double_cell"        // Output directory for deduplicated BAM files
    occ_bam_dir = "/home/prom/africa/angola/200_GBMs_merge_bam/repeat_double_cell/occ_bam"     // Output directory for ROI BAM files
    threads = 12                                    // Number of threads to use
    roi_bed = "/path/to/OCC.protein_coding.bed"    // Path to ROI BED file
    sample_id_file = "/path/to/sample_id.txt"      // File containing sample IDs
}

// Process configuration
process {
    errorStrategy = 'retry'
    maxRetries = 3
    
    withLabel: bam_processing {
        cpus = params.threads
        memory = '16GB'
    }
}

// Process to remove duplicates from BAM files
process dedup_bams {
    tag "${sample_id}"
    label 'bam_processing'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), 
          path("${sample_id}.dedup.bam"), 
          path("${sample_id}.dedup.bam.bai"), 
          path("${sample_id}.dedup.metrics.txt"), emit: dedup_bam

    script:
    """
    # Sort by coordinates for duplicate marking
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    
    # Mark and remove duplicates
    samtools markdup -@ ${task.cpus} -r -s \
        -f ${sample_id}.dedup.metrics.txt \
        ${sample_id}.sorted.bam \
        ${sample_id}.dedup.bam
    
    # Index the final BAM file
    samtools index -@ ${task.cpus} ${sample_id}.dedup.bam
    
    # Clean up intermediate files
    rm ${sample_id}.sorted.bam
    """
}

// Process to extract regions of interest
process extract_roi {
    tag "${sample_id}"
    label 'bam_processing'
    publishDir "${params.occ_bam_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(metrics), path(roi_bed)

    output:
    tuple val(sample_id), 
          path("${sample_id}.occ.bam"), 
          path("${sample_id}.occ.bam.bai")

    script:
    """
    # Extract regions of interest
    intersectBed -a ${bam} -b ${roi_bed} > ${sample_id}.occ.bam
    
    # Index the ROI BAM file
    samtools index -@ ${task.cpus} ${sample_id}.occ.bam
    """
}

// Main workflow
workflow {
    // Read sample IDs from file
    sample_ids = file(params.sample_id_file).readLines().collect { it.trim() }
    log.info "Found ${sample_ids.size()} sample IDs"

    // Create channel for BAM files
    bam_files = Channel
        .fromPath("${params.merge_bam_folder}/*.bam")
        .map { bam -> 
            def sample_id = bam.simpleName.replaceAll(/\..*/, '')
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                error "Missing index file for ${bam}"
            }
            if (!(sample_id in sample_ids)) {
                return null
            }
            tuple(sample_id, bam, bai)
        }
        .filter { it != null }

    // Create channel for ROI BED file
    roi_bed = Channel.fromPath(params.roi_bed)

    // Run deduplication process
    dedup_results = dedup_bams(bam_files)

    // Prepare input for ROI extraction
    roi_input = dedup_results.dedup_bam
        .combine(roi_bed)
        .map { sample_id, bam, bai, metrics, bed -> 
            tuple(sample_id, bam, bai, metrics, bed)
        }

    // Run ROI extraction process
    extract_roi(roi_input)
}

// Completion handler
workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Duration           : ${workflow.duration}
    Success           : ${workflow.success}
    Exit status       : ${workflow.exitStatus}
    """
    
    if (!workflow.success) {
        log.error "Pipeline failed! Check the logs for details."
    }
}

// Error handler
workflow.onError {
    log.error """
    Pipeline execution failed!
    Error message: ${workflow.errorMessage}
    """
}
