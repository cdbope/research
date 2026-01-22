#!/usr/bin/env Rscript

################################################################################
# Generate Synthetic CNV Test Data
#
# Purpose: Create synthetic CNV VCF files for testing the visualization pipeline
#          when real data has insufficient coverage
################################################################################

suppressPackageStartupMessages({
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-o", "--output-dir"), type="character", default=NULL,
              help="Output directory for synthetic VCF files", metavar="character"),
  make_option(c("-s", "--sample-id"), type="character", default="T001",
              help="Sample identifier", metavar="character"),
  make_option(c("-b", "--bin-size"), type="integer", default=100000,
              help="Bin size (default: 100000)", metavar="integer"),
  make_option(c("-n", "--num-timepoints"), type="integer", default=3,
              help="Number of time points to generate (default: 3)", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$`output-dir`)) {
  print_help(opt_parser)
  stop("Output directory is required", call.=FALSE)
}

OUTPUT_DIR <- opt$`output-dir`
SAMPLE_ID <- opt$`sample-id`
BIN_SIZE <- opt$`bin-size`
NUM_TIMEPOINTS <- opt$`num-timepoints`

# Create output directory
dir.create(OUTPUT_DIR, showWarnings=FALSE, recursive=TRUE)

cat(sprintf("[INFO] Generating synthetic CNV test data\n"))
cat(sprintf("[INFO] Sample ID: %s\n", SAMPLE_ID))
cat(sprintf("[INFO] Bin size: %d bp\n", BIN_SIZE))
cat(sprintf("[INFO] Number of time points: %d\n", NUM_TIMEPOINTS))

# Chromosome sizes (hg38)
chrom_sizes <- data.frame(
  chrom = paste0("chr", c(1:22, "X", "Y")),
  size = c(248956422, 242193529, 198295559, 190214555, 181538259,
           170805979, 159345973, 145138636, 138394717, 133797422,
           135086622, 133275309, 114364328, 107043718, 101991189,
           90338345, 83257441, 80373285, 58617616, 64444167,
           46709983, 50818468, 156040895, 57227415),
  stringsAsFactors = FALSE
)

# Generate synthetic CNVs
generate_synthetic_cnvs <- function(chrom_sizes, n_cnvs = 50) {
  cnvs <- data.frame(
    chrom = character(n_cnvs),
    start = integer(n_cnvs),
    end = integer(n_cnvs),
    type = character(n_cnvs),
    qual = numeric(n_cnvs),
    stringsAsFactors = FALSE
  )

  for (i in 1:n_cnvs) {
    # Random chromosome
    chrom_idx <- sample(1:nrow(chrom_sizes), 1, prob = chrom_sizes$size / sum(chrom_sizes$size))
    chrom <- chrom_sizes$chrom[chrom_idx]
    chrom_size <- chrom_sizes$size[chrom_idx]

    # Random CNV size (100kb - 10Mb)
    cnv_size <- sample(100000:10000000, 1)

    # Random start position
    max_start <- max(1, chrom_size - cnv_size)
    start <- sample(1:max_start, 1)
    end <- start + cnv_size

    # Random type (60% deletion, 40% duplication to simulate cancer)
    type <- sample(c("deletion", "duplication"), 1, prob = c(0.6, 0.4))

    # Random quality score
    qual <- round(runif(1, 10, 100), 2)

    cnvs[i, ] <- list(chrom, start, end, type, qual)
  }

  return(cnvs)
}

# Write CNVpytor-style VCF
write_cnvpytor_vcf <- function(cnvs, vcf_file, sample_id) {
  # VCF header
  header <- c(
    "##fileformat=VCFv4.2",
    "##source=CNVpytor",
    sprintf("##fileDate=%s", format(Sys.time(), "%Y/%m/%d %H:%M:%S")),
    "##ALT=<ID=DEL,Description=\"Deletion\">",
    "##ALT=<ID=DUP,Description=\"Duplication\">",
    sprintf("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of CNV\">"),
    sprintf("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"),
    sprintf("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">"),
    sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  )

  writeLines(header, vcf_file)

  # Write CNV records
  if (nrow(cnvs) > 0) {
    for (i in 1:nrow(cnvs)) {
      cnv <- cnvs[i, ]
      svtype <- toupper(substring(cnv$type, 1, 3))
      svlen <- cnv$end - cnv$start

      info <- sprintf("END=%d;SVTYPE=%s;SVLEN=%d", cnv$end, cnv$type, svlen)

      vcf_line <- sprintf("%s\t%d\t.\tN\t<%s>\t%.2f\tPASS\t%s",
                          cnv$chrom, cnv$start, svtype, cnv$qual, info)

      write(vcf_line, file = vcf_file, append = TRUE)
    }
  }
}

# Generate time points
for (tp in 1:NUM_TIMEPOINTS) {
  # Create subdirectory for this time point
  timepoint_id <- sprintf("%s_timepoint_%d", SAMPLE_ID, tp)
  timepoint_dir <- file.path(OUTPUT_DIR, timepoint_id)
  dir.create(timepoint_dir, showWarnings=FALSE, recursive=TRUE)

  # Generate synthetic CNVs (increase CNVs over time to simulate progression)
  n_cnvs <- 20 + (tp - 1) * 15
  cnvs <- generate_synthetic_cnvs(chrom_sizes, n_cnvs = n_cnvs)

  # Write VCF file
  vcf_file <- file.path(timepoint_dir, sprintf("%s_%dbp.vcf", timepoint_id, BIN_SIZE))
  write_cnvpytor_vcf(cnvs, vcf_file, timepoint_id)

  # Create metadata JSON
  metadata_file <- file.path(timepoint_dir, sprintf("%s_metadata.json", timepoint_id))

  # Generate timestamps 1 day apart
  timestamp <- Sys.time() + (tp - 1) * 86400

  metadata <- sprintf('{
  "sample_id": "%s",
  "timepoint": %d,
  "timestamp": "%s",
  "bin_size": %d,
  "n_cnvs": %d,
  "synthetic": true
}', timepoint_id, tp, format(timestamp, "%Y-%m-%dT%H:%M:%S%z"), BIN_SIZE, n_cnvs)

  writeLines(metadata, metadata_file)

  cat(sprintf("[INFO] Generated time point %d: %d CNVs\n", tp, n_cnvs))
  cat(sprintf("       VCF: %s\n", vcf_file))
  cat(sprintf("       Metadata: %s\n", metadata_file))
}

cat(sprintf("\n[SUCCESS] Synthetic test data generated in: %s\n", OUTPUT_DIR))
cat(sprintf("[INFO] You can now run the visualization script:\n"))
cat(sprintf("       Rscript ont_timeseries_cnv_visualization.R \\\n"))
cat(sprintf("         --input-dir %s \\\n", OUTPUT_DIR))
cat(sprintf("         --output-dir %s/visualization \\\n", OUTPUT_DIR))
cat(sprintf("         --sample-id %s \\\n", SAMPLE_ID))
cat(sprintf("         --bin-size %d\n", BIN_SIZE))
