# Fastq-to-count-for-RNAseq

## Overview

This automated pipeline generates gene count tables from raw paired-end RNA-Seq data (FASTQ) for downstream differential expression analysis (e.g., DESeq2, edgeR). It integrates Bash and R scripts, adhering to modern bioinformatics best practices for efficiency, robustness, and reproducibility.

Key Features:

*   **Efficiency**: Natively processes gzipped FASTQ files and uses symbolic links for single-lane samples to optimize disk space and I/O.
*   **Full Automation**: A single command executes the entire workflow, including automatic library strandedness detection using RSeQC.
*   **Robustness**: Includes checks for reference files and indexes, preventing redundant processing and ensuring clear error handling.
*   **Reproducibility**: Designed for execution within a self-contained Conda environment, guaranteeing consistent software versions and dependencies.


## Pipeline StagesSteps

The workflow follows these automated steps:

1.  **FASTQ Preparation**: Organizes input files, creating symbolic links or consolidating multi-lane samples.
2.  **Quality Control (Raw)**: Runs **FastQC** and aggregates reports with **MultiQC**.
3.  **Read Trimming**: Uses **Trimmomatic** for adapter and quality trimming.
4.  **Quality Report**: Re-runs **FastQC** and **MultiQC** on the cleaned data.
5.  **Reference Preparation**: Downloads necessary FASTA and GTF files.
6.  **Reference Indexing**: Creates STAR genome indexes.
7.  **Alignment & Quantification**: Aligns reads with STAR, generating per-sample count files.
8.  **Strandedness Inference**: Guesses library strandedness (unstranded, stranded-forward, or stranded-reverse) using `infer_experiment.py`.
9.  **Count Matrix Aggregation**: An R script consolidates per-sample counts into a final gene-count matrix, using the detected strandedness.

## Setup and Execution

Follow these steps chronologically to set up the environment and run the pipeline.

### Step 1: Install Miniforge (One-Time Setup)

To manage software dependencies in an isolated and reproducible way, we strongly recommend using Conda. **Miniforge** is a community-driven, open-source installer for Conda that defaults to the powerful `conda-forge` channel.

If you do not have Conda or Miniforge installed, follow these steps. If you do, proceed to Step 2.

1.  **Download the Miniforge installer for Linux:**
    ```bash
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    ```
    *For macOS or other architectures, please visit the [Miniforge GitHub page](https://github.com/conda-forge/miniforge) for the correct installer.*

2.  **Run the installer script:**
    ```bash
    bash Miniforge3-Linux-x86_64.sh
    ```
    Follow the on-screen prompts. It is recommended to accept the defaults, including allowing the installer to run `conda init`.

3.  **Close and reopen your terminal** for the changes to take effect.

### Step 2: Create and Activate the Conda Environment

This step creates a self-contained environment with all the necessary tools for the pipeline.

1.  **Create the environment.** This command will download and install all required software versions that are compatible with each other.
    ```bash
    conda create -n bio_env -c conda-forge -c bioconda fastqc multiqc trimmomatic star sra-tools pigz r-base r-dplyr r-readr rseqc bc
    ```
    *Note: `bc` (basic calculator) is included for performing floating-point comparisons in Bash.*

2.  **Activate the environment.** You must do this every time you open a new terminal session to work on this project.
    ```bash
    conda activate bio_env
    ```
    Your terminal prompt should now be prefixed with `(bio_env)`.

### Step 3: Set Up Project Directory and Acquire Data

1.  **Create a main project directory and navigate into it:**
    ```bash
    mkdir rnaseq_project
    cd rnaseq_project
    ```

2.  **Create subdirectories for your raw data and results:**
    ```bash
    mkdir raw_data
    mkdir results
    ```

3.  **Download your FASTQ files into the `raw_data` directory.** Choose one of the two methods below.

    #### Option A: Direct Download (via `wget`)
    Use this method for direct FTP/HTTP links (common for ENA).

    ```bash
    # Replace these placeholders with the actual links to your files
    wget -P raw_data/ [LINK_TO_YOUR_R1_FILE.fastq.gz]
    wget -P raw_data/ [LINK_TO_YOUR_R2_FILE.fastq.gz]
    ```

    #### Option B: Downloading from NCBI SRA (via SRA-Toolkit)
    Use this method for data from NCBI's SRA using an accession number (e.g., `SRR123456`). Repeat for each sample.

    1.  **Prefetch the data:**
        ```bash
        prefetch [SAMPLE_CODE]
        ```
    2.  **Convert to gzipped FASTQ:**
        ```bash
        fasterq-dump [SAMPLE_CODE] --split-files --gzip -e 8 -p -O raw_data
        ```
    3.  **(Recommended) Clean up the SRA file:**
        ```bash
        # The .sra file is typically in ~/ncbi/public/sra/. Adjust if needed.
        rm ~/ncbi/public/sra/[SAMPLE_CODE].sra
        ```

### Step 4: Create the Pipeline Scripts

You will need to create two script files in your `rnaseq_project` directory.

1.  **Create the main Bash script (`pipeline_rnaseq.sh`):**
    ```bash
    nano pipeline_rnaseq.sh
    ```
    Copy the code from the **"Pipeline Source Code (pipeline_rnaseq.sh)"** section below. Save and exit (`Ctrl+X`, `Y`, `Enter`).

2.  **Create the R aggregation script (`aggregate_counts.R`):**
    ```bash
    nano aggregate_counts.R
    ```
    Copy the code from the **"R Script Source Code (aggregate_counts.R)"** section below. Save and exit.

### Step 5: Make Scripts Executable

1.  **Use `chmod` to grant execute permissions:**
    ```bash
    chmod +x pipeline_rnaseq.sh
    chmod +x aggregate_counts.R
    ```

### Step 6: Run the Pipeline

With your `bio_env` Conda environment active, you can now run the **single** main script. It will handle everything automatically.

1.  **Use the format:** `./pipeline_rnaseq.sh [INPUT_DIR] [OUTPUT_DIR] [THREADS]`
2.  **Choosing threads:** Use about half the available cores on your server.
    *   On a **12-core server**, `6` threads is a good choice:
        ```bash
        ./pipeline_rnaseq.sh raw_data results 6
        ```
    *   On a **32-core server**, you could use `16`:
        ```bash
        ./pipeline_rnaseq.sh raw_data results 16
        ```

The script will now perform all steps, including the automated strandedness check and final aggregation, without requiring any further input.

## How the Automated Strandedness Check Works

The STAR aligner produces a count file (`ReadsPerGene.out.tab`) with three possible count values for each gene, corresponding to different library preparation strategies. Choosing the wrong one invalidates the analysis. Instead of requiring manual intervention, this pipeline automates the choice:

1.  After alignment, the script selects the first generated BAM file as a representative sample.
2.  It runs `infer_experiment.py` (from the RSeQC package) on this BAM file. This tool analyzes how reads align to gene annotations to infer the original library type.
3.  The script captures the output and extracts the fraction of reads corresponding to a **stranded-reverse** library and a **stranded-forward** library.
4.  Using a threshold (80%), it decides which library type is dominant:
    *   If >80% of reads are reverse, it selects **column 4**.
    *   If >80% of reads are forward, it selects **column 3**.
    *   Otherwise, it defaults to **unstranded (column 2)**.
5.  This determined column number is then automatically passed to the `aggregate_counts.R` script, ensuring a scientifically sound and fully automated workflow.

---

## Pipeline Source Code (`pipeline_rnaseq.sh`)

```bash
#!/bin/bash
set -euo pipefail

# ==============================================================================
# RNA-Seq Analysis Pipeline (Fully Automated Version)
#
# This script handles all steps from raw FASTQ to a final count matrix,
# including an automated check for library strandedness.
#
# Usage:
#   ./pipeline_rnaseq.sh [INPUT_DIR] [OUTPUT_DIR] [THREADS]
# ==============================================================================

# --- Parameters ---
INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-./processed_files}"
THREADS="${3:-4}"

# --- Output Directories ---
MERGED_DIR="${OUTPUT_DIR}/01_merged_files"
FASTQC_RAW_DIR="${OUTPUT_DIR}/02_fastqc_raw"
TRIMMED_DIR="${OUTPUT_DIR}/03_trimmed_files"
FASTQC_TRIM_DIR="${OUTPUT_DIR}/04_fastqc_trimmed"
REF_GENOME_DIR="${OUTPUT_DIR}/reference_genome"
STAR_INDEX_DIR="${OUTPUT_DIR}/STAR_index"
ALIGNED_DIR="${OUTPUT_DIR}/05_aligned_reads"
FINAL_COUNTS_DIR="${OUTPUT_DIR}/06_final_counts"

# --- Reference Files ---
GENOME_FA_URL="https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GENOME_GTF_URL="https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
GENOME_FA_NAME=$(basename "$GENOME_FA_URL" .gz)
GENOME_GTF_NAME=$(basename "$GENOME_GTF_URL" .gz)

mkdir -p "$MERGED_DIR" "$FASTQC_RAW_DIR" "$TRIMMED_DIR" "$FASTQC_TRIM_DIR" \
         "$REF_GENOME_DIR" "$STAR_INDEX_DIR" "$ALIGNED_DIR" "$FINAL_COUNTS_DIR"

# ------------------------------------------------------------------------------
# 1. Function: Merge or Link lanes, preserving compression [OPTIMIZED]
# ------------------------------------------------------------------------------
prepare_and_merge_lanes() {
    echo ">> STAGE 1: Preparing FASTQ files (Optimized with symbolic links)..."
    SAMPLES=$(find "$INPUT_DIR" -type f -name "*.fastq.gz" | xargs -n 1 basename | sed -E 's/(_L[0-9]+)?_([Rr])?[12](_001)?\.fastq\.gz$//' | sort | uniq)
    if [ -z "$SAMPLES" ]; then echo "ERROR: No .fastq.gz files found in '$INPUT_DIR'" >&2; exit 1; fi
    echo "Samples identified: $SAMPLES"
    for sample in $SAMPLES; do
        echo "Processing sample: $sample"
        files_r1_array=( $(find "$INPUT_DIR" -type f \( -name "${sample}*_[Rr]1*.fastq.gz" -o -name "${sample}*_1.fastq.gz" \) | sort) )
        if [ ${#files_r1_array[@]} -eq 0 ]; then echo "WARNING: No R1/1 files found for $sample." >&2; continue
        elif [ ${#files_r1_array[@]} -gt 1 ]; then echo "  Merging ${#files_r1_array[@]} R1 file(s)..."; cat "${files_r1_array[@]}" > "${MERGED_DIR}/${sample}_R1.fastq.gz"
        else echo "  Linking 1 R1 file..."; ln -sf "$(realpath "${files_r1_array[0]}")" "${MERGED_DIR}/${sample}_R1.fastq.gz"; fi
        files_r2_array=( $(find "$INPUT_DIR" -type f \( -name "${sample}*_[Rr]2*.fastq.gz" -o -name "${sample}*_2.fastq.gz" \) | sort) )
        if [ ${#files_r2_array[@]} -eq 0 ]; then echo "WARNING: No R2/2 files found for $sample." >&2; continue
        elif [ ${#files_r2_array[@]} -gt 1 ]; then echo "  Merging ${#files_r2_array[@]} R2 file(s)..."; cat "${files_r2_array[@]}" > "${MERGED_DIR}/${sample}_R2.fastq.gz"
        else echo "  Linking 1 R2 file..."; ln -sf "$(realpath "${files_r2_array[0]}")" "${MERGED_DIR}/${sample}_R2.fastq.gz"; fi
    done
    echo "File preparation complete."
}
# ------------------------------------------------------------------------------
# 2. Function: Quality Control with FastQC and MultiQC (Raw Data)
# ------------------------------------------------------------------------------
run_fastqc_raw() {
    echo ">> STAGE 2: Running FastQC on raw data..."
    fastqc "$MERGED_DIR"/*.fastq.gz -o "$FASTQC_RAW_DIR" --threads "$THREADS"
    multiqc "$FASTQC_RAW_DIR" -o "$FASTQC_RAW_DIR" --force -n "multiqc_report_raw"
    echo "Raw data QC report generated."
}
# ------------------------------------------------------------------------------
# 3. Function: Read Trimming with Trimmomatic
# ------------------------------------------------------------------------------
run_trimming() {
    echo ">> STAGE 3: Running Trimmomatic for adapter and quality trimming..."
    ADAPTERS_PATH="TruSeq3-PE-2.fa"
    if [ ! -f "$ADAPTERS_PATH" ]; then
        echo "WARNING: Adapter file '$ADAPTERS_PATH' not found. Downloading a standard one..."
        wget -q -O $ADAPTERS_PATH https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE-2.fa
    fi
    for r1_file in "$MERGED_DIR"/*_R1.fastq.gz; do
        sample=$(basename "$r1_file" _R1.fastq.gz)
        r2_file="${MERGED_DIR}/${sample}_R2.fastq.gz"
        echo "Trimming sample: $sample"
        trimmomatic PE -threads "$THREADS" -trimlog "${TRIMMED_DIR}/trimlog_${sample}.txt" \
            "$r1_file" "$r2_file" \
            "${TRIMMED_DIR}/${sample}_1P.fastq.gz" "${TRIMMED_DIR}/${sample}_1U.fastq.gz" \
            "${TRIMMED_DIR}/${sample}_2P.fastq.gz" "${TRIMMED_DIR}/${sample}_2U.fastq.gz" \
            "ILLUMINACLIP:${ADAPTERS_PATH}:2:30:10" SLIDINGWINDOW:4:20 MINLEN:36
    done
    echo "Trimming complete."
}
# ------------------------------------------------------------------------------
# 4. Function: Quality Control with FastQC and MultiQC (Trimmed Data)
# ------------------------------------------------------------------------------
run_fastqc_trimmed() {
    echo ">> STAGE 4: Running FastQC on trimmed data..."
    fastqc "$TRIMMED_DIR"/*_1P.fastq.gz "$TRIMMED_DIR"/*_2P.fastq.gz -o "$FASTQC_TRIM_DIR" --threads "$THREADS"
    multiqc "$FASTQC_TRIM_DIR" -o "$FASTQC_TRIM_DIR" --force -n "multiqc_report_trimmed"
    echo "Trimmed data QC report generated."
}
# ------------------------------------------------------------------------------
# 5. Function: Download and Prepare Reference Genome
# ------------------------------------------------------------------------------
prepare_reference_genome() {
    echo ">> STAGE 5: Preparing reference genome..."
    if [ ! -f "$REF_GENOME_DIR/$GENOME_FA_NAME" ] || [ ! -f "$REF_GENOME_DIR/$GENOME_GTF_NAME" ]; then
        echo "Downloading genome and annotation files..."
        wget -q -O "$REF_GENOME_DIR/${GENOME_FA_NAME}.gz" "$GENOME_FA_URL"
        wget -q -O "$REF_GENOME_DIR/${GENOME_GTF_NAME}.gz" "$GENOME_GTF_URL"
        gunzip "$REF_GENOME_DIR/${GENOME_FA_NAME}.gz"
        gunzip "$REF_GENOME_DIR/${GENOME_GTF_NAME}.gz"
    else
        echo "Reference files already exist."
    fi
}
# ------------------------------------------------------------------------------
# 6. Function: Index Reference Genome with STAR
# ------------------------------------------------------------------------------
index_reference_genome() {
    echo ">> STAGE 6: Indexing genome with STAR..."
    if [ -f "$STAR_INDEX_DIR/SA" ]; then
        echo "STAR index already exists. Skipping."
    else
        STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$STAR_INDEX_DIR" \
             --genomeFastaFiles "$REF_GENOME_DIR/$GENOME_FA_NAME" --sjdbGTFfile "$REF_GENOME_DIR/$GENOME_GTF_NAME" --sjdbOverhang 99
        echo "Indexing complete."
    fi
}
# ------------------------------------------------------------------------------
# 7. Function: Align Reads and Quantify with STAR
# ------------------------------------------------------------------------------
align_reads() {
    echo ">> STAGE 7: Aligning reads and quantifying with STAR..."
    for r1_paired_file in "$TRIMMED_DIR"/*_1P.fastq.gz; do
        sample=$(basename "$r1_paired_file" _1P.fastq.gz)
        r2_paired_file="${TRIMMED_DIR}/${sample}_2P.fastq.gz"
        echo "Aligning sample: $sample"
        STAR --runMode alignReads --runThreadN "$THREADS" --genomeDir "$STAR_INDEX_DIR" \
             --readFilesIn "$r1_paired_file" "$r2_paired_file" --readFilesCommand zcat \
             --limitBAMsortRAM 20000000000 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
             --outFileNamePrefix "${ALIGNED_DIR}/${sample}_"
    done
    echo "Alignment complete."
}

# ------------------------------------------------------------------------------
# 8. Function: Automatically Determine Strandedness and Aggregate Counts
# ------------------------------------------------------------------------------
determine_strandedness_and_aggregate() {
    echo ">> STAGE 8: Automatically determining library strandedness..."

    # Find a single BAM file to test
    local test_bam_file=$(find "$ALIGNED_DIR" -name "*_Aligned.sortedByCoord.out.bam" | head -n 1)
    if [ -z "$test_bam_file" ]; then
        echo "ERROR: No BAM files found to determine strandedness." >&2
        exit 1
    fi
    echo "Using sample BAM file for test: $test_bam_file"

    # Run infer_experiment.py and capture its output
    # Redirect stderr to stdout (2>&1) to capture all output
    local output=$(infer_experiment.py -i "$test_bam_file" -r "$REF_GENOME_DIR/$GENOME_GTF_NAME" 2>&1)

    # Parse the fractions from the output, setting to 0 if not found
    local frac_reverse=$(echo "$output" | grep 'Fraction of reads explained by "1+-,1-+"' | awk '{print $NF}' || echo "0")
    local frac_forward=$(echo "$output" | grep 'Fraction of reads explained by "1++,1--,2+-"' | awk '{print $NF}' || echo "0")

    local count_col=0
    local lib_type=""

    # Decide which column to use based on a >0.80 threshold
    if (( $(echo "$frac_reverse > 0.80" | bc -l) )); then
        count_col=4
        lib_type="Stranded (Reverse)"
    elif (( $(echo "$frac_forward > 0.80" | bc -l) )); then
        count_col=3
        lib_type="Stranded (Forward)"
    else
        count_col=2
        lib_type="Unstranded"
    fi

    echo "Library type determined as: $lib_type. Using column $count_col for counts."

    # Find the path to the R script (assuming it's in the same directory as this script)
    local script_dir=$(dirname "$(realpath "$0")")
    local r_script="$script_dir/aggregate_counts.R"

    if [ ! -f "$r_script" ]; then
        echo "ERROR: Aggregation script 'aggregate_counts.R' not found in the same directory as the main script." >&2
        exit 1
    fi

    echo ">> STAGE 9: Aggregating counts into final matrix..."
    # Execute the R script with the determined column number
    "$r_script" "$ALIGNED_DIR" "$count_col" "${FINAL_COUNTS_DIR}/final_counts.tsv"
}

# ------------------------------------------------------------------------------
# Main Pipeline Function
# ------------------------------------------------------------------------------
main() {
    echo "--- STARTING RNA-SEQ PIPELINE (FULLY AUTOMATED MODE) ---"
    
    prepare_and_merge_lanes
    run_fastqc_raw
    run_trimming
    run_fastqc_trimmed
    prepare_reference_genome
    index_reference_genome
    align_reads
    determine_strandedness_and_aggregate
    
    echo -e "\n--- PIPELINE COMPLETED SUCCESSFULLY! ---\n\nKey Results:\n  - BAM Alignments: ${ALIGNED_DIR}/\n  - Final Count Matrix: ${FINAL_COUNTS_DIR}/final_counts.tsv\n  - QC Reports: ${FASTQC_RAW_DIR}/ and ${FASTQC_TRIM_DIR}/\n-----------------------------------------"
}

# ==============================================================================
# Execution Trigger
# ==============================================================================
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main
fi
```

---

## R Script Source Code (`aggregate_counts.R`)

```R
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
```

---

## Understanding the Final Output Structure

Upon successful completion of all steps, the output directory (e.g., `results/`) will contain the following structure:

*   `01_merged_files/`: Symbolic links to or merged copies of the raw FASTQ files.
*   `02_fastqc_raw/`: Quality control reports for the raw data.
*   `03_trimmed_files/`: FASTQ files after adapter and quality trimming.
*   `04_fastqc_trimmed/`: Quality control reports for the trimmed data.
*   `05_aligned_reads/`: BAM alignment files and individual STAR count files (`*ReadsPerGene.out.tab`).
*   `06_final_counts/`: **The primary scientific result.** This directory is created by the `aggregate_counts.R` script and contains the final, aggregated count matrix.
*   `reference_genome/` & `STAR_index/`: Contain the reference genome and STAR index files.
