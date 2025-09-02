# RNA-Seq Analysis Pipeline

This repository contains a fully automated, robust, and reproducible pipeline for processing raw, paired-end RNA-Seq data from FASTQ files to a final gene-count matrix.

The workflow combines the power of a Bash script for orchestrating core bioinformatics tools and an R script for accurate, final data aggregation. It is designed with best practices in mind, including full automation, reproducibility via a Conda environment, efficiency with data handling, and flexibility through customizable parameters.

## Getting Started

To use this pipeline, follow this two-step process:

1.  **Set Up the Environment**: This is a **one-time setup** to install all the necessary software in an isolated Conda environment. This step ensures that the pipeline runs reliably on any machine.
    *   **>> Proceed to the [Environment Setup Guide](./ENVIRONMENT_SETUP.md) for detailed instructions.**

2.  **Run the Pipeline**: Once the environment is configured, you can proceed to the main user manual, which details how to prepare your data and execute the pipeline.
    *   **>> Proceed to the [Pipeline Documentation](./PIPELINE_DOCUMENTATION.md) for usage instructions.**

## Repository Structure

This repository is organized into several key files. Here is a map of what each file does:

*   `README.md`
    *   You are here. This file provides a high-level overview of the project and directs you to the relevant documentation.

*   `ENVIRONMENT_SETUP.md`
    *   A detailed, step-by-step guide for the **one-time installation** of Miniforge and the creation of the specific Conda environment required to run this pipeline.

*   `PIPELINE_DOCUMENTATION.md`
    *   The main **user manual** for the pipeline. It contains comprehensive instructions on how to structure your project, download data, and execute the `pipeline.sh` script with its various parameters.

*   `environment.yml`
    *   The Conda environment **recipe file**. This text file lists all the software dependencies required by the pipeline. It is used to create a perfectly reproducible software environment, which is the cornerstone of reliable scientific analysis.

*   `pipeline.sh`
    *   The main **executable Bash script**. This script is the engine of the pipeline, orchestrating all the processing steps from quality control and trimming to alignment and quantification.

*   `aggregate_counts.R`
    *   A helper **R script** that is called automatically by `pipeline.sh`. Its specialized job is to intelligently parse the output from the aligner, select the scientifically correct data based on the library type, and generate the final, clean count matrix.
