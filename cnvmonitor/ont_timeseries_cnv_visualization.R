#!/usr/bin/env Rscript

################################################################################
# ONT Time-Series CNV Visualization Script
#
# Purpose: Create time-series visualizations of CNV data from multiple
#          time-point BAM files processed with CNVpytor
#
# Usage: Rscript ont_timeseries_cnv_visualization.R --input-dir <dir> --output-dir <dir> --sample-id <id>
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(cowplot)
  library(viridis)
  library(scales)
  library(gridExtra)
  library(ggpubr)
})

# Define command line options
option_list <- list(
  make_option(c("-i", "--input-dir"), type="character", default=NULL,
              help="Input directory containing individual BAM CNV results", metavar="character"),
  make_option(c("-o", "--output-dir"), type="character", default=NULL,
              help="Output directory for time-series plots", metavar="character"),
  make_option(c("-s", "--sample-id"), type="character", default=NULL,
              help="Sample identifier", metavar="character"),
  make_option(c("-y", "--cytoband"), type="character", default=NULL,
              help="Cytoband file (optional)", metavar="character"),
  make_option(c("-b", "--bin-size"), type="integer", default=100000,
              help="Bin size to use for visualization (default: 100000)", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$`input-dir`) || is.null(opt$`output-dir`) || is.null(opt$`sample-id`)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

INPUT_DIR <- opt$`input-dir`
OUTPUT_DIR <- opt$`output-dir`
SAMPLE_ID <- opt$`sample-id`
CYTOBAND_FILE <- opt$cytoband
BIN_SIZE <- opt$`bin-size`

# Create output directory
dir.create(OUTPUT_DIR, showWarnings=FALSE, recursive=TRUE)

cat(sprintf("[INFO] %s - Starting time-series CNV visualization\n", Sys.time()))
cat(sprintf("[INFO] Input directory: %s\n", INPUT_DIR))
cat(sprintf("[INFO] Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("[INFO] Sample ID: %s\n", SAMPLE_ID))

################################################################################
# Function: Parse CNVpytor VCF output
################################################################################
parse_cnvpytor_vcf <- function(vcf_file) {
  if (!file.exists(vcf_file)) {
    warning(sprintf("VCF file not found: %s", vcf_file))
    return(NULL)
  }

  # Read VCF file
  lines <- readLines(vcf_file)

  # Skip header lines
  data_lines <- lines[!grepl("^#", lines)]

  if (length(data_lines) == 0) {
    warning(sprintf("No CNV calls found in: %s", vcf_file))
    return(NULL)
  }

  # Parse VCF fields
  cnv_data <- lapply(data_lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]

    if (length(fields) < 8) return(NULL)

    # Extract INFO field
    info <- fields[8]

    # Parse INFO field
    end_match <- regmatches(info, regexpr("END=([0-9]+)", info))
    end <- as.numeric(sub("END=", "", end_match))

    cnvtype_match <- regmatches(info, regexpr("SVTYPE=([^;]+)", info))
    cnvtype <- sub("SVTYPE=", "", cnvtype_match)

    # Calculate copy number from GT field if available
    gt_field <- if (length(fields) >= 10) fields[10] else "."

    data.frame(
      chrom = fields[1],
      start = as.numeric(fields[2]),
      end = end,
      type = cnvtype,
      qual = as.numeric(fields[6]),
      info = info,
      stringsAsFactors = FALSE
    )
  })

  # Combine into data frame
  cnv_df <- do.call(rbind, cnv_data[!sapply(cnv_data, is.null)])

  return(cnv_df)
}

################################################################################
# Function: Extract metadata from JSON files
################################################################################
extract_metadata <- function(metadata_file) {
  if (!file.exists(metadata_file)) {
    return(NULL)
  }

  # Simple JSON parsing for our specific format
  json_lines <- readLines(metadata_file)
  timestamp_line <- json_lines[grepl("timestamp", json_lines)]

  if (length(timestamp_line) > 0) {
    timestamp <- sub('.*"timestamp":\\s*"([^"]+)".*', "\\1", timestamp_line)
    return(as.POSIXct(timestamp))
  }

  # Fallback to file modification time
  return(file.info(metadata_file)$mtime)
}

################################################################################
# Main Analysis
################################################################################

cat("[INFO] Searching for processed BAM directories...\n")

# Find all subdirectories with CNV results
sample_dirs <- list.dirs(INPUT_DIR, recursive=FALSE)
sample_dirs <- sample_dirs[grepl(SAMPLE_ID, basename(sample_dirs))]

if (length(sample_dirs) == 0) {
  stop(sprintf("No sample directories found for sample ID: %s", SAMPLE_ID))
}

cat(sprintf("[INFO] Found %d time point(s)\n", length(sample_dirs)))

# Collect CNV data from all time points
all_cnv_data <- list()
time_points <- c()

for (sample_dir in sample_dirs) {
  cat(sprintf("[INFO] Processing: %s\n", basename(sample_dir)))

  # Find VCF file for the specified bin size
  vcf_pattern <- sprintf("*_%dbp.vcf", BIN_SIZE)
  vcf_files <- list.files(sample_dir, pattern=glob2rx(vcf_pattern), full.names=TRUE)

  if (length(vcf_files) == 0) {
    warning(sprintf("No VCF found in %s", sample_dir))
    next
  }

  vcf_file <- vcf_files[1]

  # Parse CNV data
  cnv_data <- parse_cnvpytor_vcf(vcf_file)

  if (is.null(cnv_data) || nrow(cnv_data) == 0) {
    warning(sprintf("No CNV data extracted from %s", vcf_file))
    next
  }

  # Extract timestamp
  metadata_file <- list.files(sample_dir, pattern="*_metadata.json", full.names=TRUE)
  if (length(metadata_file) > 0) {
    timestamp <- extract_metadata(metadata_file[1])
  } else {
    timestamp <- file.info(sample_dir)$mtime
  }

  cnv_data$timepoint <- timestamp
  cnv_data$timepoint_label <- basename(sample_dir)

  all_cnv_data[[length(all_cnv_data) + 1]] <- cnv_data
  time_points <- c(time_points, as.character(timestamp))
}

if (length(all_cnv_data) == 0) {
  stop("No CNV data could be extracted from any time point")
}

# Combine all CNV data
combined_cnv <- do.call(rbind, all_cnv_data)

# Sort by time
combined_cnv <- combined_cnv %>%
  arrange(timepoint)

# Create time index
combined_cnv$time_index <- as.numeric(factor(combined_cnv$timepoint))

cat(sprintf("[INFO] Total CNV events across all time points: %d\n", nrow(combined_cnv)))
cat(sprintf("[INFO] Time points: %d\n", length(unique(combined_cnv$timepoint))))

################################################################################
# Visualization 1: CNV count over time
################################################################################
cat("[INFO] Generating CNV count over time plot...\n")

cnv_counts <- combined_cnv %>%
  group_by(timepoint, time_index, type) %>%
  summarise(count = n(), .groups = "drop")

p1 <- ggplot(cnv_counts, aes(x = time_index, y = count, color = type, group = type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("deletion" = "#d62728", "duplication" = "#2ca02c"),
                     labels = c("Deletions", "Duplications")) +
  labs(title = sprintf("CNV Events Over Time - %s", SAMPLE_ID),
       subtitle = sprintf("Bin size: %d bp", BIN_SIZE),
       x = "Time Point",
       y = "Number of CNV Events",
       color = "CNV Type") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, sprintf("%s_cnv_count_timeseries.png", SAMPLE_ID)),
       plot = p1, width = 12, height = 6, dpi = 300)

################################################################################
# Visualization 2: Genome-wide CNV heatmap over time
################################################################################
cat("[INFO] Generating genome-wide CNV heatmap...\n")

# Create genomic bins for visualization
combined_cnv$genomic_position <- combined_cnv$start

# Chromosome order
chrom_order <- paste0("chr", c(1:22, "X", "Y"))
combined_cnv$chrom <- factor(combined_cnv$chrom, levels = chrom_order)

# Assign copy number values based on type
combined_cnv$copy_number <- ifelse(combined_cnv$type == "duplication", 1,
                                   ifelse(combined_cnv$type == "deletion", -1, 0))

p2 <- ggplot(combined_cnv, aes(x = genomic_position / 1e6, y = time_index, fill = copy_number)) +
  geom_tile() +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0,
                       name = "Copy Number\nChange") +
  labs(title = sprintf("Genome-wide CNV Profile Over Time - %s", SAMPLE_ID),
       x = "Genomic Position (Mb)",
       y = "Time Point") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    strip.text = element_text(size = 8),
    legend.position = "right"
  )

ggsave(file.path(OUTPUT_DIR, sprintf("%s_genome_wide_cnv_heatmap.png", SAMPLE_ID)),
       plot = p2, width = 20, height = 8, dpi = 300)

################################################################################
# Visualization 3: Per-chromosome CNV profiles
################################################################################
cat("[INFO] Generating per-chromosome CNV profiles...\n")

# Focus on autosomes and sex chromosomes
main_chroms <- paste0("chr", c(1:22, "X"))
chrom_data <- combined_cnv %>%
  filter(chrom %in% main_chroms)

if (nrow(chrom_data) > 0) {
  p3 <- ggplot(chrom_data, aes(x = start / 1e6, y = copy_number, color = as.factor(time_index))) +
    geom_point(alpha = 0.6, size = 1.5) +
    facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
    scale_color_viridis_d(name = "Time Point") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = sprintf("Per-Chromosome CNV Profiles - %s", SAMPLE_ID),
         x = "Genomic Position (Mb)",
         y = "Copy Number Change") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      strip.background = element_rect(fill = "lightblue")
    )

  ggsave(file.path(OUTPUT_DIR, sprintf("%s_per_chromosome_cnv.png", SAMPLE_ID)),
         plot = p3, width = 16, height = 12, dpi = 300)
}

################################################################################
# Visualization 4: CNV size distribution over time
################################################################################
cat("[INFO] Generating CNV size distribution plot...\n")

combined_cnv$cnv_size <- (combined_cnv$end - combined_cnv$start) / 1e6  # Size in Mb

p4 <- ggplot(combined_cnv, aes(x = cnv_size, fill = type)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~ time_index, scales = "free_y") +
  scale_fill_manual(values = c("deletion" = "#d62728", "duplication" = "#2ca02c"),
                    labels = c("Deletions", "Duplications")) +
  scale_x_log10(labels = comma) +
  labs(title = sprintf("CNV Size Distribution Over Time - %s", SAMPLE_ID),
       x = "CNV Size (Mb, log scale)",
       y = "Count",
       fill = "CNV Type") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, sprintf("%s_cnv_size_distribution.png", SAMPLE_ID)),
       plot = p4, width = 14, height = 10, dpi = 300)

################################################################################
# Generate summary report
################################################################################
cat("[INFO] Generating summary report...\n")

summary_file <- file.path(OUTPUT_DIR, sprintf("%s_timeseries_summary.txt", SAMPLE_ID))

sink(summary_file)
cat("========================================\n")
cat("Time-Series CNV Analysis Summary\n")
cat("========================================\n")
cat(sprintf("Sample ID: %s\n", SAMPLE_ID))
cat(sprintf("Bin Size: %d bp\n", BIN_SIZE))
cat(sprintf("Analysis Date: %s\n", Sys.time()))
cat(sprintf("Number of Time Points: %d\n", length(unique(combined_cnv$timepoint))))
cat("========================================\n\n")

# Summary by time point
cat("CNV Counts by Time Point:\n")
cat("----------------------------------------\n")
summary_by_time <- combined_cnv %>%
  group_by(timepoint_label, type) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = count, values_fill = 0)

print(as.data.frame(summary_by_time))
cat("\n")

# Overall statistics
cat("Overall CNV Statistics:\n")
cat("----------------------------------------\n")
cat(sprintf("Total CNV events: %d\n", nrow(combined_cnv)))
cat(sprintf("Total deletions: %d\n", sum(combined_cnv$type == "deletion")))
cat(sprintf("Total duplications: %d\n", sum(combined_cnv$type == "duplication")))
cat(sprintf("Mean CNV size: %.2f Mb\n", mean(combined_cnv$cnv_size)))
cat(sprintf("Median CNV size: %.2f Mb\n", median(combined_cnv$cnv_size)))
cat("\n")

sink()

cat(sprintf("[SUCCESS] Summary report saved: %s\n", summary_file))
cat("[SUCCESS] Time-series CNV visualization completed!\n")
cat(sprintf("[INFO] Output files saved to: %s\n", OUTPUT_DIR))
