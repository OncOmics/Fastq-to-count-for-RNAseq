# Fastq to Count Matrix for bulk RNAseq

[![DOI](https://zenodo.org/badge/DOI/YOUR_DOI_HERE.svg)](https://doi.org/YOUR_DOI_HERE)
<!--
NOTE: This DOI badge is a placeholder. After your first release,
you will get a real DOI from Zenodo to replace the link and text here.
-->

## Overview

This automated pipeline generates gene count matrices from raw, paired-end RNA-Seq data (FASTQ) for downstream differential expression analysis (e.g., with DESeq2 or edgeR). It integrates Bash and R scripts, adhering to modern bioinformatics best practices for efficiency, robustness, and reproducibility.

### Key Features

*   **Efficiency**: Natively processes gzipped FASTQ files and uses symbolic links for single-lane samples to optimize disk space and I/O.
*   **Full Automation**: A single command executes the entire workflow, including automatic library strandedness detection using RSeQC.
*   **Robustness**: Includes checks for reference files and indexes, preventing redundant processing and ensuring clear error handling.
*   **Reproducibility**: Uses a YAML file (`environment.yml`) to ensure a perfectly reproducible software environment.

---

## How to Use This Repository

This project is designed to cater to different user needs. Choose the path that best suits you.

#### Path 1: The Guided Tour (For a Complete Understanding)
If you are new to the pipeline or wish to understand its components in detail before execution, we recommend following this structured learning path.

1.  **Understand the Environment**: Begin by reading the setup guide to learn how to install the necessary software in a controlled and reproducible way.
    *   **>> Proceed to the [Environment Setup Guide](./ENVIRONMENT_SETUP.md)**

2.  **Understand the Pipeline**: Next, read the main user manual. It explains the project's structure, the logic behind each function, and the details of every customizable parameter.
    *   **>> Proceed to the [Pipeline Documentation](./PIPELINE_DOCUMENTATION.md)**

#### Path 2: The Quick Start (For a Rapid Test Run)
If you are an experienced user and want to quickly test the pipeline to ensure it works in your system, follow the steps in this section.

1.  **Clone the repository and navigate into it.**
    ```bash
    git clone https://github.com/OncOmics/Fastq-to-count-for-RNAseq.git
    cd Fastq-to-count-for-RNAseq
    ```

2.  **Create and activate the Conda environment.** This project recommends **Mamba**, a fast, parallel re-implementation of Conda.
    ```bash
    # Create the environment from the provided recipe file
    mamba env create -f environment.yml

    # Activate the environment
    conda activate bio_env
    ```

3.  **Download sample data.** This command creates the necessary directories and downloads a test dataset.
    ```bash
    mkdir -p raw_data results
    wget -P raw_data/ https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/000/ERR2208890/ERR2208890_1.fastq.gz
    wget -P raw_data/ https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/000/ERR2208890/ERR2208890_2.fastq.gz
    ```

4.  **Run the pipeline.** Execute the main script with standard parameters.
    ```bash
    ./pipeline.sh --input raw_data --output results --threads 8
    ```
    If the pipeline completes successfully, your setup is correct, and you are ready to use it on your own data.

---

## Repository Structure

This repository is organized into several key files. Here is a map of what each file does:

*   `README.md`: You are here. A high-level overview and navigational guide for the project.
*   `ENVIRONMENT_SETUP.md`: A detailed, step-by-step guide for the **one-time installation** of Miniforge and the creation of the specific Conda environment required to run this pipeline.
*   `PIPELINE_DOCUMENTATION.md`: The main **user manual** for the pipeline. It contains comprehensive instructions on how to structure your project, download data, and execute the `pipeline.sh` script with all its parameters.
*   `environment.yml`: The Conda environment **recipe file**. This YAML file lists all the software dependencies required by the pipeline.
*   `pipeline.sh`: The main **executable Bash script**. This script is the engine of the pipeline, orchestrating all the processing steps.
*   `aggregate_counts.R`: An **R script** that is called automatically by `pipeline.sh`. It parses the output from the aligner, selects the scientifically correct data based on the library type, and generates the final, clean count matrix.

## Citation

PANEPUCCI, E.M.; ANZOLINI-CASSIANO, M.H. (2025). Fastq to Count Matrix for bulk RNAseq. GitHub. Retrieved from https://github.com/OncOmics/Fastq-to-count-for-RNAseq.
