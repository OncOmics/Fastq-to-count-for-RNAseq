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