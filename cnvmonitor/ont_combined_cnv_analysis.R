#!/usr/bin/env Rscript

################################################################################
# ONT Combined CNV Analysis Script
#
# Purpose: Generate comprehensive CNV analysis combining all time points
#          Identifies persistent CNVs, emerging CNVs, and transient CNVs
#
# Usage: Rscript ont_combined_cnv_analysis.R --input-dir <dir> --output-dir <dir> --sample-id <id>
################################################################################

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

# Command line options
option_list <- list(
  make_option(c("-i", "--input-dir"), type="character", default=NULL,
              help="Input directory containing individual BAM CNV results", metavar="character"),
  make_option(c("-o", "--output-dir"), type="character", default=NULL,
              help="Output directory for combined analysis", metavar="character"),
  make_option(c("-s", "--sample-id"), type="character", default=NULL,
              help="Sample identifier", metavar="character"),
  make_option(c("-y", "--cytoband"), type="character", default=NULL,
              help="Cytoband file (optional)", metavar="character"),
  make_option(c("-b", "--bin-size"), type="integer", default=100000,
              help="Bin size to use (default: 100000)", metavar="integer"),
  make_option(c("--overlap-threshold"), type="numeric", default=0.5,
              help="Overlap threshold for CNV matching (default: 0.5)", metavar="numeric")
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
BIN_SIZE <- opt$`bin-size`
OVERLAP_THRESH <- opt$`overlap-threshold`

dir.create(OUTPUT_DIR, showWarnings=FALSE, recursive=TRUE)

cat(sprintf("[INFO] %s - Starting combined CNV analysis\n", Sys.time()))

################################################################################
# Helper Functions
################################################################################

parse_cnvpytor_vcf <- function(vcf_file) {
  if (!file.exists(vcf_file)) {
    return(NULL)
  }

  lines <- readLines(vcf_file)
  data_lines <- lines[!grepl("^#", lines)]

  if (length(data_lines) == 0) {
    return(NULL)
  }

  cnv_data <- lapply(data_lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 8) return(NULL)

    info <- fields[8]
    end_match <- regmatches(info, regexpr("END=([0-9]+)", info))
    end <- as.numeric(sub("END=", "", end_match))

    cnvtype_match <- regmatches(info, regexpr("SVTYPE=([^;]+)", info))
    cnvtype <- sub("SVTYPE=", "", cnvtype_match)

    data.frame(
      chrom = fields[1],
      start = as.numeric(fields[2]),
      end = end,
      type = cnvtype,
      qual = as.numeric(fields[6]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, cnv_data[!sapply(cnv_data, is.null)])
}

extract_timestamp <- function(metadata_file) {
  if (!file.exists(metadata_file)) {
    return(NULL)
  }

  json_lines <- readLines(metadata_file)
  timestamp_line <- json_lines[grepl("timestamp", json_lines)]

  if (length(timestamp_line) > 0) {
    timestamp <- sub('.*"timestamp":\\s*"([^"]+)".*', "\\1", timestamp_line)
    return(as.POSIXct(timestamp))
  }

  return(file.info(metadata_file)$mtime)
}

# Calculate overlap between two CNV regions
calculate_overlap <- function(start1, end1, start2, end2) {
  overlap_start <- max(start1, start2)
  overlap_end <- min(end1, end2)
  overlap_length <- max(0, overlap_end - overlap_start)

  region1_length <- end1 - start1
  region2_length <- end2 - start2

  if (region1_length == 0 || region2_length == 0) {
    return(0)
  }

  # Return reciprocal overlap
  min(overlap_length / region1_length, overlap_length / region2_length)
}

################################################################################
# Load all CNV data
################################################################################

cat("[INFO] Loading CNV data from all time points...\n")

sample_dirs <- list.dirs(INPUT_DIR, recursive=FALSE)
sample_dirs <- sample_dirs[grepl(SAMPLE_ID, basename(sample_dirs))]

if (length(sample_dirs) == 0) {
  stop(sprintf("No sample directories found for sample ID: %s", SAMPLE_ID))
}

cat(sprintf("[INFO] Found %d time point(s)\n", length(sample_dirs)))

all_cnv_list <- list()

for (i in seq_along(sample_dirs)) {
  sample_dir <- sample_dirs[i]
  cat(sprintf("[INFO] Loading: %s\n", basename(sample_dir)))

  vcf_pattern <- sprintf("*_%dbp.vcf", BIN_SIZE)
  vcf_files <- list.files(sample_dir, pattern=glob2rx(vcf_pattern), full.names=TRUE)

  if (length(vcf_files) == 0) {
    warning(sprintf("No VCF found in %s", sample_dir))
    next
  }

  cnv_data <- parse_cnvpytor_vcf(vcf_files[1])

  if (is.null(cnv_data) || nrow(cnv_data) == 0) {
    warning(sprintf("No CNV data from %s", vcf_files[1]))
    next
  }

  metadata_file <- list.files(sample_dir, pattern="*_metadata.json", full.names=TRUE)
  if (length(metadata_file) > 0) {
    timestamp <- extract_timestamp(metadata_file[1])
  } else {
    timestamp <- file.info(sample_dir)$mtime
  }

  cnv_data$timepoint <- timestamp
  cnv_data$timepoint_id <- i
  cnv_data$timepoint_label <- basename(sample_dir)
  cnv_data$cnv_id <- paste0("CNV_", i, "_", 1:nrow(cnv_data))

  all_cnv_list[[i]] <- cnv_data
}

if (length(all_cnv_list) == 0) {
  stop("No CNV data could be loaded")
}

combined_cnv <- do.call(rbind, all_cnv_list)
combined_cnv <- combined_cnv %>% arrange(timepoint)

cat(sprintf("[INFO] Total CNVs loaded: %d\n", nrow(combined_cnv)))
cat(sprintf("[INFO] Time points: %d\n", length(unique(combined_cnv$timepoint))))

################################################################################
# CNV Classification: Persistent vs Emerging vs Transient
################################################################################

cat("[INFO] Classifying CNVs across time points...\n")

# Function to find matching CNVs across time points
find_matching_cnvs <- function(cnv_data, overlap_threshold = 0.5) {
  cnv_groups <- list()
  group_id <- 1
  assigned <- rep(FALSE, nrow(cnv_data))

  for (i in 1:nrow(cnv_data)) {
    if (assigned[i]) next

    current_cnv <- cnv_data[i, ]
    group <- c(i)
    assigned[i] <- TRUE

    # Look for matching CNVs in other time points
    for (j in (i+1):nrow(cnv_data)) {
      if (j > nrow(cnv_data) || assigned[j]) next

      other_cnv <- cnv_data[j, ]

      # Must be same chromosome and type
      if (current_cnv$chrom != other_cnv$chrom ||
          current_cnv$type != other_cnv$type) {
        next
      }

      # Calculate overlap
      overlap <- calculate_overlap(
        current_cnv$start, current_cnv$end,
        other_cnv$start, other_cnv$end
      )

      if (overlap >= overlap_threshold) {
        group <- c(group, j)
        assigned[j] <- TRUE
      }
    }

    cnv_groups[[group_id]] <- group
    group_id <- group_id + 1
  }

  return(cnv_groups)
}

cnv_groups <- find_matching_cnvs(combined_cnv, OVERLAP_THRESH)

# Assign group IDs to CNVs
combined_cnv$cnv_group <- NA
for (i in seq_along(cnv_groups)) {
  combined_cnv$cnv_group[cnv_groups[[i]]] <- i
}

# Classify CNVs based on their persistence
cnv_classification <- combined_cnv %>%
  group_by(cnv_group) %>%
  summarise(
    num_timepoints = n_distinct(timepoint),
    first_timepoint = min(timepoint_id),
    last_timepoint = max(timepoint_id),
    chrom = first(chrom),
    start = min(start),
    end = max(end),
    type = first(type),
    .groups = "drop"
  )

total_timepoints <- length(unique(combined_cnv$timepoint))

cnv_classification$classification <- with(cnv_classification, {
  case_when(
    num_timepoints == total_timepoints ~ "Persistent",
    num_timepoints == 1 ~ "Transient",
    TRUE ~ "Emerging"
  )
})

cat("[INFO] CNV Classification Summary:\n")
print(table(cnv_classification$classification))

# Add classification back to combined data
combined_cnv <- combined_cnv %>%
  left_join(cnv_classification %>% select(cnv_group, classification),
            by = "cnv_group")

################################################################################
# Visualization 1: CNV Classification Overview
################################################################################

cat("[INFO] Generating CNV classification plot...\n")

class_summary <- cnv_classification %>%
  group_by(classification, type) %>%
  summarise(count = n(), .groups = "drop")

p1 <- ggplot(class_summary, aes(x = classification, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("deletion" = "#d62728", "duplication" = "#2ca02c"),
                    labels = c("Deletions", "Duplications")) +
  labs(title = sprintf("CNV Classification - %s", SAMPLE_ID),
       subtitle = sprintf("%d time points analyzed", total_timepoints),
       x = "CNV Classification",
       y = "Number of CNV Regions",
       fill = "CNV Type") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, sprintf("%s_cnv_classification.png", SAMPLE_ID)),
       plot = p1, width = 10, height = 6, dpi = 300)

################################################################################
# Visualization 2: Persistent CNV Regions
################################################################################

cat("[INFO] Generating persistent CNV map...\n")

persistent_cnvs <- cnv_classification %>%
  filter(classification == "Persistent")

if (nrow(persistent_cnvs) > 0) {
  chrom_order <- paste0("chr", c(1:22, "X", "Y"))
  persistent_cnvs$chrom <- factor(persistent_cnvs$chrom, levels = chrom_order)

  persistent_cnvs$copy_number <- ifelse(persistent_cnvs$type == "duplication", 1, -1)

  p2 <- ggplot(persistent_cnvs, aes(x = start / 1e6, xend = end / 1e6,
                                     y = 0, yend = 0, color = type)) +
    geom_segment(size = 4) +
    facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
    scale_color_manual(values = c("deletion" = "#d62728", "duplication" = "#2ca02c"),
                       labels = c("Deletions", "Duplications")) +
    labs(title = sprintf("Persistent CNV Regions - %s", SAMPLE_ID),
         subtitle = sprintf("CNVs present in all %d time points", total_timepoints),
         x = "Genomic Position (Mb)",
         y = "",
         color = "CNV Type") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "lightblue")
    )

  ggsave(file.path(OUTPUT_DIR, sprintf("%s_persistent_cnv_map.png", SAMPLE_ID)),
         plot = p2, width = 16, height = 12, dpi = 300)
}

################################################################################
# Visualization 3: CNV Evolution Tracking
################################################################################

cat("[INFO] Generating CNV evolution plot...\n")

# Select top CNV groups by size for tracking
top_cnvs <- cnv_classification %>%
  mutate(size = end - start) %>%
  arrange(desc(size)) %>%
  head(20)

if (nrow(top_cnvs) > 0) {
  tracking_data <- combined_cnv %>%
    filter(cnv_group %in% top_cnvs$cnv_group)

  tracking_data$region_label <- with(tracking_data,
    paste0(chrom, ":", round(start/1e6, 1), "-", round(end/1e6, 1), "Mb"))

  p3 <- ggplot(tracking_data, aes(x = timepoint_id, y = region_label, fill = type)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = c("deletion" = "#d62728", "duplication" = "#2ca02c"),
                      labels = c("Deletions", "Duplications")) +
    labs(title = sprintf("CNV Evolution - Top 20 Regions by Size - %s", SAMPLE_ID),
         x = "Time Point",
         y = "Genomic Region",
         fill = "CNV Type") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )

  ggsave(file.path(OUTPUT_DIR, sprintf("%s_cnv_evolution_tracking.png", SAMPLE_ID)),
         plot = p3, width = 12, height = 10, dpi = 300)
}

################################################################################
# Export Results
################################################################################

cat("[INFO] Exporting results...\n")

# Export CNV classification table
class_output <- cnv_classification %>%
  mutate(size_mb = (end - start) / 1e6) %>%
  arrange(chrom, start)

write.csv(class_output,
          file.path(OUTPUT_DIR, sprintf("%s_cnv_classification.csv", SAMPLE_ID)),
          row.names = FALSE)

# Export persistent CNVs
if (nrow(persistent_cnvs) > 0) {
  persistent_output <- persistent_cnvs %>%
    mutate(size_mb = (end - start) / 1e6) %>%
    arrange(chrom, start)

  write.csv(persistent_output,
            file.path(OUTPUT_DIR, sprintf("%s_persistent_cnvs.csv", SAMPLE_ID)),
            row.names = FALSE)
}

# Generate summary report
summary_file <- file.path(OUTPUT_DIR, sprintf("%s_combined_analysis_summary.txt", SAMPLE_ID))

sink(summary_file)
cat("========================================\n")
cat("Combined CNV Analysis Summary\n")
cat("========================================\n")
cat(sprintf("Sample ID: %s\n", SAMPLE_ID))
cat(sprintf("Bin Size: %d bp\n", BIN_SIZE))
cat(sprintf("Analysis Date: %s\n", Sys.time()))
cat(sprintf("Number of Time Points: %d\n", total_timepoints))
cat(sprintf("Overlap Threshold: %.2f\n", OVERLAP_THRESH))
cat("========================================\n\n")

cat("CNV Classification:\n")
cat("----------------------------------------\n")
print(table(cnv_classification$classification))
cat("\n")

cat("CNV Types by Classification:\n")
cat("----------------------------------------\n")
print(table(cnv_classification$classification, cnv_classification$type))
cat("\n")

if (nrow(persistent_cnvs) > 0) {
  cat("Persistent CNV Statistics:\n")
  cat("----------------------------------------\n")
  cat(sprintf("Total persistent CNV regions: %d\n", nrow(persistent_cnvs)))
  cat(sprintf("Persistent deletions: %d\n", sum(persistent_cnvs$type == "deletion")))
  cat(sprintf("Persistent duplications: %d\n", sum(persistent_cnvs$type == "duplication")))
  cat(sprintf("Total genomic span: %.2f Mb\n", sum(persistent_cnvs$end - persistent_cnvs$start) / 1e6))
  cat("\n")
}

sink()

cat(sprintf("[SUCCESS] Summary report saved: %s\n", summary_file))
cat("[SUCCESS] Combined CNV analysis completed!\n")
cat(sprintf("[INFO] Output files saved to: %s\n", OUTPUT_DIR))
