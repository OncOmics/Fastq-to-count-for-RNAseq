#!/usr/bin/env Rscript

# ==============================================================================
# Aggregate STAR Gene Counts
#
# This script reads all 'ReadsPerGene.out.tab' files from a specified directory,
# selects the appropriate count column based on library strandedness, and
# merges them into a single, analysis-ready count matrix.
#
# Author: [Seu Nome/Laboratório]
# Date: [Data Atual]
#
# Usage:
#   ./aggregate_counts.R <star_counts_dir> <count_column> <output_file.tsv>
#
# Example:
#   ./aggregate_counts.R ./results/05_aligned_reads 4 ./results/06_final_counts/final_counts.tsv
# ==============================================================================

# --- Load Libraries ---
# Suppress startup messages for a cleaner command-line experience
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# --- Parse Command-Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
  stop("Usage: ./aggregate_counts.R <star_counts_dir> <count_column> <output_file.tsv>", call. = FALSE)
}

input_dir <- args[1]
count_col <- as.integer(args[2])
output_file <- args[3]

# --- Validate Column Number ---
# Ensure the provided column makes scientific sense
if (!count_col %in% c(2, 3, 4)) {
  stop("Error: Count column must be 2 (unstranded), 3 (stranded-forward), or 4 (stranded-reverse).", call. = FALSE)
}

# --- Find All Count Files ---
count_files <- list.files(path = input_dir, pattern = "ReadsPerGene.out.tab$", full.names = TRUE)
if (length(count_files) == 0) {
  stop(paste("Error: No '*ReadsPerGene.out.tab' files found in the directory:", input_dir), call. = FALSE)
}

# --- Define Column Names Explicitly ---
# This makes the script robust to STAR versions that may or may not include a header row.
col_names_def <- c("GeneID", "Unstranded", "Stranded_Fwd", "Stranded_Rev")

# --- Read and Process All Files in a Loop ---
all_counts <- lapply(count_files, function(file) {
  # Extract a clean sample name from the filename
  sample_name <- sub("_ReadsPerGene.out.tab$", "", basename(file))
  
  # Read the STAR output file:
  # - comment = "N_": Skips the initial summary lines that start with "N_"
  # - col_names: Manually assigns our defined column names
  # - col_types: Specifies data types for faster reading (c=character, i=integer)
  read_tsv(
    file,
    comment = "N_",
    col_names = col_names_def,
    col_types = "ciii"
  ) %>%
    # Select only the GeneID and the count column specified by the user
    select(GeneID, all_of(count_col)) %>%
    # Rename the count column to the actual sample name for the final matrix
    rename(!!sample_name := colnames(.)[2])
})

# --- Join All Data Frames into a Single Matrix ---
# The Reduce function iteratively merges each data frame in the 'all_counts' list
# using a full_join on the 'GeneID' column.
final_matrix <- Reduce(function(x, y) full_join(x, y, by = "GeneID"), all_counts)

# --- Write the Final Matrix to a File ---
# Create the output directory if it doesn't already exist
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(final_matrix, output_file)

cat("Final count matrix created successfully at:", output_file, "\n")
