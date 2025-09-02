# Environment Setup Guide

This guide provides detailed, one-time instructions for setting up the Conda environment required to run the RNA-Seq pipeline. Following these steps will ensure that you have all the necessary software, with the correct versions, to guarantee a reproducible and error-free analysis.

## The Importance of a Controlled Environment

Modern bioinformatics relies on many different software tools, each with its own specific dependencies. Managing these manually can lead to version conflicts and non-reproducible results (the dreaded "it works on my machine" problem).

We solve this using **Conda**, a package and environment management system. It allows us to create an isolated, self-contained environment with all the tools needed for this pipeline. The recipe for this environment is stored in the `environment.yml` file, making the entire software setup perfectly reproducible for anyone, on any machine.

## Step 1: Install Miniforge

If you already have a Conda distribution (like Anaconda or Miniconda) installed, you can skip to Step 2.

We recommend **Miniforge** because it is a lightweight, open-source installer for Conda that defaults to the community-driven `conda-forge` channel, a standard for scientific software.

1.  **Download the Miniforge installer for Linux.** Open your terminal and run the following command:
    ```bash
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    ```
    *Note: For macOS or other architectures, please visit the [Miniforge GitHub page](https://github.com/conda-forge/miniforge) to find the correct installer link.*

2.  **Run the installer script.** This will start the installation process.
    ```bash
    bash Miniforge3-Linux-x86_64.sh
    ```
    You will be guided through a series of prompts. We recommend accepting the default options by pressing **Enter** at each step. When asked if you want to initialize Miniforge3, it is important to answer **yes**, as this will configure your shell to use Conda.

3.  **Close and reopen your terminal.** This is a critical step to ensure that the changes to your system's configuration are loaded. After reopening, you should see `(base)` at the beginning of your terminal prompt, indicating that Conda is active.

## Step 2: Create the Conda Environment from the YAML File

Now that Conda is installed, you can create the specific environment for this pipeline using the provided `environment.yml` recipe file.

1.  **Navigate to the project directory** where you have cloned or downloaded the repository files (including `environment.yml`).

2.  **Create the environment.** This single command will read the `environment.yml` file, automatically resolve all dependencies, and install every tool listed within a new environment named `bio_env`.
    ```bash
    conda env create -f environment.yml
    ```
    This process will take several minutes as Conda downloads and installs all the necessary packages.

## Step 3: Activate the Environment

Once the environment is created, you must **activate** it to start using the installed tools.

1.  **Run the activation command:**
    ```bash
    conda activate bio_env
    ```
2.  **Verify activation.** Your terminal prompt should now change to be prefixed with `(bio_env)`, like this:
    ```
    (bio_env) your_username@your_machine:~$
    ```
    This prefix confirms that you are now inside the dedicated `bio_env` environment and that all the tools required for the pipeline (`fastqc`, `star`, `rseqc`, etc.) are ready to be used.

    *Note: You will need to run `conda activate bio_env` every time you open a new terminal session to work on this project.*

**Setup is now complete!** You can proceed to the main user guide to run your analysis.

**>> Next Step: [Pipeline Documentation](./PIPELINE_DOCUMENTATION.md)**
