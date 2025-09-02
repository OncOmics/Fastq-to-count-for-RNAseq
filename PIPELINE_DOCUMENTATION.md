# Pipeline Documentation and User Manual

## Introduction

This document is the main user guide for the RNA-Seq Analysis Pipeline. It provides a comprehensive, step-by-step walkthrough of how to structure your project, create the necessary scripts, download data, and execute a full analysis.

This guide assumes you have already set up your Conda environment as detailed in the **[Environment Setup Guide](./ENVIRONMENT_SETUP.md)**. Please ensure your `bio_env` environment is active before proceeding:
```bash
conda activate bio_env
```

## The Two-Script Philosophy

This pipeline is built on a powerful two-script philosophy that separates orchestration from data manipulation, using the best tool for each job:

1.  **`pipeline.sh` (The Orchestrator)**: A Bash script that acts as the "engine" of the pipeline. Its job is to call external bioinformatics programs (`fastqc`, `star`, etc.) in the correct order, manage files and directories, and handle parallel processing. It is excellent at automation.

2.  **`aggregate_counts.R` (The Specialist)**: An R script that acts as a specialized tool for a single, critical task: accurately aggregating the gene counts. R is far superior to Bash for data manipulation, allowing us to robustly parse the aligner's output and build a clean, analysis-ready matrix.

You will create both of these scripts in the following steps.

## Step 1: Set Up Your Project Directory

A clean directory structure is essential. The pipeline is designed to work with the following layout.

1.  **Create a main project directory** for your analysis and navigate into it:
    ```bash
    mkdir my_rnaseq_analysis
    cd my_rnaseq_analysis
    ```

2.  **Create subdirectories for raw data and results.** The pipeline will read from `raw_data` and write all outputs into `results`.
    ```bash
    mkdir raw_data
    mkdir results
    ```
    Your final structure will eventually contain the scripts you create and the `environment.yml` file from this repository.

## Step 2: Create the Pipeline Scripts

You will now create the two core scripts for the pipeline. The source code for both is provided in this repository in the files `pipeline.sh` and `aggregate_counts.R`.

#### 2.1 The R Script (`aggregate_counts.R`)

This script's job is to intelligently parse the output from STAR, select the correct count data based on the library type, and generate the final matrix.

1.  **Create the file using a text editor:**
    ```bash
    nano aggregate_counts.R
    ```
2.  **Copy the entire contents** of the `aggregate_counts.R` file from this repository and paste it into the nano editor.
3.  **Save and exit** (`Ctrl+X`, then `Y`, then `Enter`).

#### 2.2 The Main Bash Script (`pipeline.sh`)

This is the main executable that runs the entire workflow.

1.  **Create the file:**
    ```bash
    nano pipeline.sh
    ```
2.  **Copy the entire contents** of the `pipeline.sh` file from this repository and paste it into the nano editor.
3.  **Save and exit.**

#### 2.3 Make the Scripts Executable

Finally, you must give the system permission to run these files as programs.
```bash
chmod +x pipeline.sh
chmod +x aggregate_counts.R
```

## Step 3: Acquire Your Raw Data

Place all your raw, paired-end FASTQ files (in `.fastq.gz` format) inside the `raw_data` directory.

#### Option A: Direct Download (via `wget`)
```bash
wget -P raw_data/ https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/000/ERR2208890/ERR2208890_1.fastq.gz
wget -P raw_data/ https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/000/ERR2208890/ERR2208890_2.fastq.gz
```

#### Option B: Downloading from NCBI SRA (via SRA-Toolkit)
```bash
prefetch SRR123456
fasterq-dump SRR123456 --split-files --gzip -e 8 -p -O raw_data
rm ~/ncbi/public/sra/SRR123456.sra
```

## Step 4: Run the Pipeline

The `pipeline.sh` script is designed with professional command-line features to be both user-friendly and highly flexible.

### Understanding the Command Structure and Design

The script was intentionally designed to move beyond simple positional arguments (like `$1`, `$2`) to a more robust "flag-based" system. Here’s why this is a superior approach:

*   **Professional Argument Parsing (Flags)**: Instead of relying on a strict order of arguments, the script uses named flags (e.g., `--input`, `-t`). This is the standard for professional command-line tools.
    *   **Why?** It makes the command self-documenting and less prone to user error. You don't have to remember if threads come before or after the output directory, and you can supply arguments in any order.

*   **Clear Default Values**: All optional arguments have sensible defaults (e.g., 4 threads, 20GB RAM).
    *   **Why?** This allows a user to run the pipeline with a minimal command, making it easy to get started. Advanced users can then override these defaults as needed without modifying the script's code.

*   **Customizable Tool Parameters**: Critical parameters for internal tools, which were previously "hardcoded," are now exposed as command-line flags.
    *   `--ram`: Controls the `--limitBAMsortRAM` parameter in STAR.
    *   `--overhang`: Controls the `--sjdbOverhang` parameter in STAR.
    *   **Why?** This provides essential flexibility. A user on a high-memory server can increase the RAM limit for better performance, while a user on a resource-constrained machine can lower it to prevent crashes. Similarly, the overhang can be tuned precisely to the read length of a specific experiment.

*   **Built-in Help Message**: The script includes a `--help` flag.
    *   **Why?** This makes the script self-contained. A user can quickly see all available options and how to use them without needing to read the entire script or documentation.

### Command Parameters

| Flag | Argument | Description | Required? | Default |
| :--- | :--- | :--- | :--- | :--- |
| `-i`, `--input` | Path | Directory containing the raw `.fastq.gz` files. | **Yes** | N/A |
| `-o`, `--output`| Path | Directory where all results will be saved. | **Yes** | N/A |
| `-t`, `--threads` | Integer | Number of CPU threads to use for parallel tasks. | No | `4` |
| `-r`, `--ram` | String | Memory limit for STAR BAM sorting (e.g., `30G`). | No | `20G` |
| `-s`, `--overhang`| Integer | Value for STAR's `sjdbOverhang` (Rule: `read_length - 1`). | No | `99` |
| `-h`, `--help` | N/A | Displays the help message and exits. | N/A | N/A |

### Example Executions

#### Basic Usage
This is the simplest way to run the pipeline, using the default settings for memory and overhang.

```bash
./pipeline.sh --input raw_data --output results --threads 8
```

#### Advanced Usage (Customizing Parameters)
This example is tailored for a more powerful server and longer reads, overriding the default RAM and overhang values.

```bash
# Example for a 32-core server with 64GB RAM, processing 150bp reads
./pipeline.sh \
  --input raw_data \
  --output results \
  --threads 16 \
  --ram 50G \
  --overhang 149
```

## Detailed Breakdown of Pipeline Functions

The `pipeline.sh` script is composed of several modular functions. Here is a description of what each one does.

*   `prepare_and_merge_lanes()`
    Prepares raw FASTQ files. It detects if a sample is split across multiple files (lanes) and concatenates them. If a sample has only one file per read, it creates efficient **symbolic links** instead of making slow file copies, saving time and disk space.

*   `run_fastqc_raw()`
    Performs initial quality control. It runs **FastQC** on all prepared FASTQ files and then uses **MultiQC** to aggregate the individual reports into a single HTML summary.

*   `run_trimming()`
    Cleans the raw reads using **Trimmomatic** to remove Illumina adapter sequences and trim low-quality bases from the ends, improving downstream alignment accuracy.

*   `run_fastqc_trimmed()`
    Re-runs **FastQC** and **MultiQC** on the cleaned FASTQ files, allowing for a direct "before and after" comparison to verify trimming effectiveness.

*   `prepare_reference_genome()`
    Ensures the necessary reference files (FASTA genome and GTF annotation) are available, downloading them from Ensembl if not found.

*   `index_reference_genome()`
    Creates a STAR index from the reference genome, a required pre-processing step for the aligner. It checks if an index already exists to avoid re-running this time-consuming step. The `sjdbOverhang` parameter (customizable via `--overhang`) is used here.

*   `align_reads()`
    The core alignment step. It maps reads to the reference genome using STAR. During this process, it also performs the initial gene-level quantification (`--quantMode GeneCounts`) and generates a sorted BAM file. The memory for this step is controlled by the `--ram` flag.

*   `determine_strandedness_and_aggregate()`
    The final, intelligent step. It automatically runs `infer_experiment.py` on a sample BAM file to determine the library's strandedness. It then calls the `aggregate_counts.R` script, passing it the correct column number to use for building the final, scientifically accurate count matrix.

## Understanding the Dual-Execution Mode

For advanced users and debugging, the pipeline supports an interactive execution mode. This is made possible by the "Execution Trigger" at the end of the `pipeline.sh` script.

*   **Automated Execution (`./pipeline.sh ...`)**: When run directly, the script executes the `main` function, running all steps in sequence.
*   **Interactive Execution (`source pipeline.sh ...`)**: When you `source` the script, it loads all functions and variables into your terminal but **does not** run `main`. This allows you to execute each function manually, one at a time, to inspect intermediate results or debug a specific step.

## Understanding the Final Output Structure

Upon successful completion, the output directory (e.g., `results/`) will contain the following structure:

*   `01_merged_files/`: Symbolic links to or merged copies of the raw FASTQ files.
*   `02_fastqc_raw/`: Quality control reports for the raw data.
*   `03_trimmed_files/`: FASTQ files after adapter and quality trimming.
*   `04_fastqc_trimmed/`: Quality control reports for the trimmed data.
*   `05_aligned_reads/`: BAM alignment files and individual STAR count files (`*ReadsPerGene.out.tab`).
*   `06_final_counts/`: **The primary scientific result.** This directory is created by the `aggregate_counts.R` script and contains the final, aggregated count matrix (`final_counts.tsv`).
*   `reference_genome/` & `STAR_index/`: Contain the reference genome and STAR index files.
