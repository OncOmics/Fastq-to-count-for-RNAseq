#!/usr/bin/env Rscript

# --- Load Libraries ---
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
if (!count_col %in% c(2, 3, 4)) {
  stop("Count column must be 2 (unstranded), 3 (stranded-forward), or 4 (stranded-reverse).", call. = FALSE)
}

# --- Find All Count Files ---
count_files <- list.files(path = input_dir, pattern = "ReadsPerGene.out.tab$", full.names = TRUE)
if (length(count_files) == 0) {
  stop("No '*ReadsPerGene.out.tab' files found in the specified directory.", call. = FALSE)
}

# --- Define Column Names Explicitly ---
col_names_def <- c("GeneID", "Unstranded", "Stranded_Fwd", "Stranded_Rev")

# --- Read and Process All Files ---
all_counts <- lapply(count_files, function(file) {
  sample_name <- sub("_ReadsPerGene.out.tab$", "", basename(file))
  
  read_tsv(
    file,
    comment = "N_", # Skip the summary lines at the top
    col_names = col_names_def, # Assign column names manually
    col_types = "ciii"
  ) %>%
    select(GeneID, all_of(count_col)) %>%
    rename(!!sample_name := colnames(.)[2]) # Rename count column to sample name
})

# --- Join All Data Frames into a Single Matrix ---
final_matrix <- Reduce(function(x, y) full_join(x, y, by = "GeneID"), all_counts)

# --- Write the Final Matrix ---
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(final_matrix, output_file)

cat("Final count matrix created successfully at:", output_file, "\n")
