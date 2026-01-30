#!/usr/bin/env Rscript
# Compare beta values between 5mC only and 5mC+5hmC merged data for sample T19-138
# Creates curve plots showing beta value distributions with probes on x-axis

library(ggplot2)
library(dplyr)
library(tidyr)

# File paths
file_5mc <- "/home/chbope/extension/script/drexler/results/200gbm/deconv/T19-138_deconv_betas.csv"
file_5mc_5hmc <- "/home/chbope/extension/script/drexler/results/200gbm/deconv/mergemh/T19-138_deconv_betas_merged.csv"
output_dir <- "/home/chbope/extension/script/drexler/epigenetic-neural-glioblastoma-main/code/DNAm_deconv/meth_atlas-master"

# Read data
cat("Loading 5mC data...\n")
data_5mc <- read.csv(file_5mc, row.names = 1)
colnames(data_5mc) <- "beta_5mC"

cat("Loading 5mC+5hmC data...\n")
data_5mc_5hmc <- read.csv(file_5mc_5hmc, row.names = 1)
colnames(data_5mc_5hmc) <- "beta_5mC_5hmC"

# Add probe IDs as column
data_5mc$CpG <- rownames(data_5mc)
data_5mc_5hmc$CpG <- rownames(data_5mc_5hmc)

cat(sprintf("5mC probes: %d\n", nrow(data_5mc)))
cat(sprintf("5mC+5hmC probes: %d\n", nrow(data_5mc_5hmc)))

# Merge on common probes
merged_data <- merge(data_5mc, data_5mc_5hmc, by = "CpG")
cat(sprintf("Common probes: %d\n", nrow(merged_data)))

# Sort by 5mC beta values for better visualization
merged_data <- merged_data[order(merged_data$beta_5mC), ]
merged_data$probe_index <- 1:nrow(merged_data)

# Convert to long format for ggplot
long_data <- merged_data %>%
  pivot_longer(
    cols = c(beta_5mC, beta_5mC_5hmC),
    names_to = "Method",
    values_to = "Beta"
  ) %>%
  mutate(Method = recode(Method,
    "beta_5mC" = "5mC only",
    "beta_5mC_5hmC" = "5mC + 5hmC"
  ))

# Create curve plot with probes on x-axis (sorted by 5mC value)
p1 <- ggplot(long_data, aes(x = probe_index, y = Beta, color = Method)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("5mC only" = "#E74C3C", "5mC + 5hmC" = "#3498DB")) +
  labs(
    title = "T19-138: Beta Value Comparison (5mC vs 5mC+5hmC)",
    subtitle = sprintf("Common probes: %d (sorted by 5mC beta value)", nrow(merged_data)),
    x = "Probe Index (sorted by 5mC beta)",
    y = "Beta Value",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "top"
  )

# Save plot
output_file <- file.path(output_dir, "T19-138_5mC_vs_5mC5hmC_curve_plot.png")
ggsave(output_file, p1, width = 12, height = 6, dpi = 300)
cat(sprintf("Saved curve plot: %s\n", output_file))

# Create scatter plot for direct comparison
p2 <- ggplot(merged_data, aes(x = beta_5mC, y = beta_5mC_5hmC)) +
  geom_point(alpha = 0.3, color = "#2C3E50", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = "T19-138: 5mC vs 5mC+5hmC Beta Value Correlation",
    subtitle = sprintf("N = %d probes, Correlation = %.3f",
                       nrow(merged_data),
                       cor(merged_data$beta_5mC, merged_data$beta_5mC_5hmC)),
    x = "5mC Beta Value",
    y = "5mC + 5hmC Beta Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

# Save scatter plot
output_file2 <- file.path(output_dir, "T19-138_5mC_vs_5mC5hmC_scatter.png")
ggsave(output_file2, p2, width = 8, height = 8, dpi = 300)
cat(sprintf("Saved scatter plot: %s\n", output_file2))

# Create density curve plot
p3 <- ggplot(long_data, aes(x = Beta, fill = Method, color = Method)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("5mC only" = "#E74C3C", "5mC + 5hmC" = "#3498DB")) +
  scale_color_manual(values = c("5mC only" = "#E74C3C", "5mC + 5hmC" = "#3498DB")) +
  labs(
    title = "T19-138: Beta Value Distribution",
    subtitle = sprintf("Common probes: %d", nrow(merged_data)),
    x = "Beta Value",
    y = "Density",
    fill = "Method",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "top"
  )

# Save density plot
output_file3 <- file.path(output_dir, "T19-138_5mC_vs_5mC5hmC_density.png")
ggsave(output_file3, p3, width = 10, height = 6, dpi = 300)
cat(sprintf("Saved density plot: %s\n", output_file3))

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("\n5mC only:\n")
print(summary(merged_data$beta_5mC))

cat("\n5mC + 5hmC:\n")
print(summary(merged_data$beta_5mC_5hmC))

cat("\nDifference (5mC+5hmC - 5mC):\n")
diff_values <- merged_data$beta_5mC_5hmC - merged_data$beta_5mC
print(summary(diff_values))

cat(sprintf("\nCorrelation: %.4f\n", cor(merged_data$beta_5mC, merged_data$beta_5mC_5hmC)))
cat(sprintf("Mean absolute difference: %.4f\n", mean(abs(diff_values))))

cat("\nDone!\n")
