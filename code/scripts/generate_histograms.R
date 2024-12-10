#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)

# Function to read sequence lengths from a file
read_lengths <- function(file) {
  lengths <- scan(file, what = numeric(), quiet = TRUE)
  return(lengths)
}

# Function to read GC percentages from a file
read_gc <- function(file) {
  gc_content <- scan(file, what = numeric(), quiet = TRUE)
  return(gc_content)
}

# Plotting function for histograms
plot_histogram <- function(data, title, xlab, file_name, log_scale = FALSE) {
  p <- ggplot(data, aes(x = value)) +
    geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
    labs(title = title, x = xlab, y = "Frequency") +
    theme_minimal()
  
  if (log_scale) {
    p <- p + scale_x_log10()
  }
  
  ggsave(file_name, plot = p, width = 8, height = 6)
}

# File paths for sequence lengths and GC content
le_100kb_lengths_file <- "data/processed/sequences_le_100kb_lengths.txt"
gt_100kb_lengths_file <- "data/processed/sequences_gt_100kb_lengths.txt"
le_100kb_gc_file <- "data/processed/sequences_le_100kb_gc.txt"
gt_100kb_gc_file <- "data/processed/sequences_gt_100kb_gc.txt"

# Read data
le_100kb_lengths <- read_lengths(le_100kb_lengths_file)
gt_100kb_lengths <- read_lengths(gt_100kb_lengths_file)
le_100kb_gc <- read_gc(le_100kb_gc_file)
gt_100kb_gc <- read_gc(gt_100kb_gc_file)

# Convert to data frames
le_100kb_lengths_df <- data.frame(value = le_100kb_lengths)
gt_100kb_lengths_df <- data.frame(value = gt_100kb_lengths)
le_100kb_gc_df <- data.frame(value = le_100kb_gc)
gt_100kb_gc_df <- data.frame(value = gt_100kb_gc)

# Ensure the output directory exists
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# Generate histograms
plot_histogram(le_100kb_lengths_df, "Sequence Length Distribution (≤ 100kb)", "Length (bp)", "output/figures/le_100kb_length_histogram.png", log_scale = TRUE)
plot_histogram(gt_100kb_lengths_df, "Sequence Length Distribution (> 100kb)", "Length (bp)", "output/figures/gt_100kb_length_histogram.png", log_scale = TRUE)
plot_histogram(le_100kb_gc_df, "GC% Distribution (≤ 100kb)", "GC%", "output/figures/le_100kb_gc_histogram.png")
plot_histogram(gt_100kb_gc_df, "GC% Distribution (> 100kb)", "GC%", "output/figures/gt_100kb_gc_histogram.png")