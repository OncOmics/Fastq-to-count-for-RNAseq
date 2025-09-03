# Fastq to Count Matrix for bulk RNAseq

## Overview

  This standard automated pipeline generates gene count tables from raw paired-end RNA-Seq data (FASTQ) for downstream differential expression analysis (e.g., DESeq2, edgeR). It integrates Bash and R scripts, adhering to modern bioinformatics best practices for efficiency, robustness, and reproducibility.

Key Features:

*  **Efficiency**: Natively processes gzipped FASTQ files and uses symbolic links for single-lane samples to optimize disk space and I/O.
*  **Full Automation**: A single command executes the entire workflow, including automatic library strandedness detection using RSeQC.
*  **Robustness**: Includes checks for reference files and indexes, preventing redundant processing and ensuring clear error handling.
*  **Reproducibility**: Designed for execution within a self-contained Conda environment, guaranteeing consistent software versions and dependencies.

## Pipeline Steps

In brief, the workflow follows these steps:

1.  **FASTQ Preparation**: Organizes input files, creating symbolic links or consolidating multi-lane samples.
2.  **Quality Control (Raw)**: Runs **FastQC** and aggregates reports with **MultiQC**.
3.  **Read Trimming**: Uses **Trimmomatic** for adapter and quality trimming.
4.  **Quality Report**: Re-runs **FastQC** and **MultiQC** on the cleaned data.
5.  **Reference Preparation**: Downloads necessary FASTA and GTF files.
6.  **Reference Indexing**: Creates STAR genome indexes.
7.  **Alignment & Quantification**: Aligns reads with STAR, generating per-sample count files.
8.  **Strandedness Inference**: Guesses library strandedness (unstranded, stranded-forward, or stranded-reverse) using `infer_experiment.py`.
9.  **Count Matrix Aggregation**: An R script consolidates per-sample counts into a final gene-count matrix, using the detected strandedness.

## Getting Started

To use this pipeline, follow this two-step process:

1.  **Clone the repository**: Once the environment is configured, download this repo, either by clicking the download button or cloning with `git`.

	```
	git clone https://github.com/OncOmics/Fastq-to-count-for-RNAseq.git
	cd Fastq-to-count-for-RNAseq.git
	```

2.  **Set up the environment**: This is a **one-time setup** to install all the necessary software in an isolated Conda environment. This step should ensure that the pipeline runs on any machine.

* **Create an enviroment using our yaml recipe**
Make sure you have any conda manager installation, create the environment and activate it.
	```
	mamba env create -f ./environment.yml
	conda activate count_env
	```
*  **Or proceed to the [environment setup guide](./ENVIRONMENT_SETUP.md) for detailed instructions.**

## Quick Start

1. **Place you fastq files in a folder**.

For example, you can download raw data from NCBI's SRA and dump fastqs.

```
mkdir raw_data
wget -P raw_data/ https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/000/ERR2208890/ERR2208890_1.fastq.gz
wget -P raw_data/ https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR220/000/ERR2208890/ERR2208890_2.fastq.gz
 ```

2.  **Run the Pipeline**: Once everithing is set:

```
./pipeline.sh --input raw_data --output results --threads 8
```

> **Consider proceeding to the [Pipeline Documentation](./PIPELINE_DOCUMENTATION.md) for detailled usage instructions.**

## Repository Structure

This repository is organized into several key files. Here is a map of what each file does:

*  `README.md`: You are here. This file provides a high-level overview of the project and directs you to the relevant documentation.
*  `ENVIRONMENT_SETUP.md`:  A detailed, step-by-step guide for the **one-time installation** of Miniforge and the creation of the specific environment required to run this pipeline.
*  `PIPELINE_DOCUMENTATION.md`:  The main **user manual** for the pipeline. It contains comprehensive instructions on how to structure your project, download data, and execute the `pipeline.sh` script with all its parameters.
*  `environment.yml`:  The conda environment **recipe file**. This yaml file lists all the software dependencies required by the pipeline.
*  `pipeline.sh`:  The main **executable bash script**. This script is the engine of the pipeline, orchestrating all the processing steps.
*  `aggregate_counts.R`: A **R script** that is called automatically by `pipeline.sh`. It parses the output from the aligner, select the scientifically correct data based on the library type, and generate the final, clean count matrix.

## Citation

PANEPUCCI, E.M.; ANZOLINI-CASSIANO, M.H. (2025). Fastq to Count Matrix for bulk RNAseq. GitHub. Retrieved from https://github.com/OncOmics/Fastq-to-count-for-RNAseq.
