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

The `pipeline.sh` script is designed with two execution modes for maximum flexibility.

### Mode 1: Automated Full Execution (Recommended for Production)

This mode runs the entire pipeline from start to finish automatically. It is the standard way to process your data.

#### Command Structure and Parameters

| Flag | Argument | Description | Required? | Default |
| :--- | :--- | :--- | :--- | :--- |
| `-i`, `--input` | Path | Directory containing the raw `.fastq.gz` files. | **Yes** | N/A |
| `-o`, `--output`| Path | Directory where all results will be saved. | **Yes** | N/A |
| `-t`, `--threads` | Integer | Number of CPU threads to use for parallel tasks. | No | `4` |
| `-r`, `--ram` | String | Memory limit for STAR BAM sorting (e.g., `30G`). | No | `20G` |
| `-s`, `--overhang`| Integer | Value for STAR's `sjdbOverhang` (Rule: `read_length - 1`). | No | `99` |
| `-h`, `--help` | N/A | Displays the help message and exits. | N/A | N/A |

#### Example Execution

```bash
./pipeline.sh --input raw_data --output results --threads 8
```

### Mode 2: Interactive Step-by-Step Execution (for Debugging)

This mode is designed for debugging or running specific parts of the pipeline. By using the `source` command, you load all the script's functions into your terminal without executing them, allowing you to call them manually.

1.  **First, source the script, providing the arguments as you normally would:**
    ```bash
    source pipeline.sh --input raw_data --output results --threads 8
    ```
    This will set up all the necessary variables but will **not** start the pipeline.

2.  **Now, you can execute each function individually.** For example:
    ```bash
    # Run only the first two steps to check initial QC
    prepare_and_merge_lanes
    run_fastqc_raw

    # After inspecting the QC, you can continue
    run_trimming
    run_fastqc_trimmed

    # Or, if you want to run the entire sequence from this point
    main
    ```

### How Dual-Execution Works: The "Execution Trigger"

This flexibility is made possible by the final block in the `pipeline.sh` script:
```bash
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main
fi
```
*   When you run the script directly (`./pipeline.sh`), the script's name (`$0`) is the same as its source file path (`${BASH_SOURCE[0]}`). The condition is **true**, and the `main` function is automatically executed.
*   When you `source` the script, your terminal is the program executing (`$0`), which is different from the script's file path. The condition is **false**, and the `main` function is **not** executed, leaving the functions ready for you to call manually.

## Understanding the Pipeline Functions

The `pipeline.sh` script is composed of several modular functions. Here is a detailed description of what each one does.

*   `prepare_and_merge_lanes()`
    This function prepares the raw FASTQ files. It intelligently detects if a sample is split across multiple files (lanes). If so, it concatenates them. If a sample has only one R1 and one R2 file, it creates efficient **symbolic links** instead of making slow file copies, saving significant time and disk space.

*   `run_fastqc_raw()`
    This function performs the initial quality control step. It runs **FastQC** on all prepared FASTQ files and then uses **MultiQC** to aggregate all the individual reports into a single, easy-to-read HTML summary.

*   `run_trimming()`
    This function cleans the raw reads. It uses **Trimmomatic** to remove any remaining Illumina adapter sequences and to trim low-quality bases from the ends of the reads, improving the accuracy of the downstream alignment.

*   `run_fastqc_trimmed()`
    After trimming, this function re-runs **FastQC** and **MultiQC** on the cleaned FASTQ files. This allows you to directly compare the "before" and "after" quality reports to verify that the trimming step was effective.

*   `prepare_reference_genome()`
    This function ensures the necessary reference files are available. It checks for the reference genome (FASTA) and gene annotation (GTF) files and downloads them from Ensembl if they are not found.

*   `index_reference_genome()`
    Before alignment, the aligner (STAR) requires the reference genome to be in a special, pre-processed format called an index. This function runs the `STAR --runMode genomeGenerate` command to create this index. It also checks if an index already exists to avoid re-running this time-consuming step. The `sjdbOverhang` parameter, customizable via the `--overhang` flag, is used here.

*   `align_reads()`
    This is the core alignment step. The function iterates through each sample's trimmed FASTQ files and uses `STAR --runMode alignReads` to map them to the reference genome. During this process, it also performs the initial gene-level quantification (`--quantMode GeneCounts`) and generates a sorted BAM file for each sample. The memory for this step is controlled by the `--ram` flag.

*   `determine_strandedness_and_aggregate()`
    This is the final, intelligent step. The function automatically runs `infer_experiment.py` (from RSeQC) on a sample BAM file to determine the library's strandedness. It then parses the result to decide whether the data is unstranded, stranded-forward, or stranded-reverse. Finally, it calls the `aggregate_counts.R` script, passing it the correct column number to use for building the final, scientifically accurate count matrix.

## Understanding the Final Output Structure

Upon successful completion, the output directory (e.g., `results/`) will contain the following structure:

*   `01_merged_files/`: Symbolic links to or merged copies of the raw FASTQ files.
*   `02_fastqc_raw/`: Quality control reports for the raw data.
*   `03_trimmed_files/`: FASTQ files after adapter and quality trimming.
*   `04_fastqc_trimmed/`: Quality control reports for the trimmed data.
*   `05_aligned_reads/`: BAM alignment files and individual STAR count files (`*ReadsPerGene.out.tab`).
*   `06_final_counts/`: **The primary scientific result.** This directory is created by the `aggregate_counts.R` script and contains the final, aggregated count matrix (`final_counts.tsv`).
*   `reference_genome/` & `STAR_index/`: Contain the reference genome and STAR index files.
