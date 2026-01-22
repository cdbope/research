#!/bin/bash
# Job name:
#SBATCH --job-name=Nanopore
# Project:
#SBATCH --account=p1753_tsd
# Wall clock limit:
#SBATCH --time=200:00:00
# Max memory usage:
#SBATCH --mem-per-cpu=8G
# Number of cores:
#SBATCH --cpus-per-task=16

# module load
module purge
module load Java/11.0.20
export _JAVA_OPTIONS='-XX:ParallelGCThreads=1'
export NXF_OFFLINE='TRUE'
set -o errexit

#### start job ###
cd $SCRATCH


### merge multiple small bam to a single bam
/ess/p1753/cluster/Software/nWGS_pipeline/nextflow-23.10.2-all run /ess/p1753/cluster/Software/nWGS_pipeline/main.nf \
   -c /ess/p1753/cluster/Software/nWGS_pipeline/conf/mergebam.config \
   --input_dir /ess/p1753/data/durable/Projects/Nanopore/IPD0091-DXX-p01-F08/merged_dir \
   -with-apptainer \
   --max_cpus 16 --max_memory '128.GB' --max_time '166.h' \
   --run_mode_mergebam &&

### run epi analysis
/ess/p1753/cluster/Software/nWGS_pipeline/nextflow-23.10.2-all run /ess/p1753/cluster/Software/nWGS_pipeline/main.nf \
   -c /ess/p1753/cluster/Software/nWGS_pipeline/conf/epi2me.config \
   --input_dir /ess/p1753/data/durable/Projects/Nanopore/IPD0091-DXX-p01-F08/merged_dir \
   --outdir /ess/p1753/data/durable/Projects/Nanopore/IPD0091-DXX-p01-F08/analysis_output \
   -with-apptainer \
   --max_cpus 16 --max_memory '128.GB' --max_time '166.h' \
   --run_mode_epi2me all &&

### run postprocessing analysis
/ess/p1753/cluster/Software/nWGS_pipeline/nextflow-23.10.2-all run /ess/p1753/cluster/Software/nWGS_pipeline/main.nf \
   -c /ess/p1753/cluster/Software/nWGS_pipeline/conf/analysis.config \
   --input_dir /ess/p1753/data/durable/Projects/Nanopore/IPD0091-DXX-p01-F08/merged_dir \
   --outdir /ess/p1753/data/durable/Projects/Nanopore/IPD0091-DXX-p01-F08/analysis_output \
   -with-apptainer \
   -with-report "report.html" \
   -with-timeline "timeline.html" \
   -with-trace "trace.txt" \
   --max_cpus 16 --max_memory '128.GB' --max_time '166.h' \
   --run_mode_analysis all

