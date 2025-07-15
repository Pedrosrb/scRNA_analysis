#!/bin/bash

# --- Configuration ---
# Path to the text file containing a list of SRA accession numbers (one per line)
# SRR20791199 -> AM
# SRR20791200 -> AM
# SRR20791203 -> AM
# SRR20791201 -> AM
# SRR20791205 -> AM
# SRR20791206 -> AM
# SRR20791198 -> CM
# SRR20791197 -> CM
# SRR20791204 -> CM
SRA_LIST="sra_list.txt"

# Path to the Cell Ranger reference transcriptome
TRANSCRIPTOME_PATH="/path/to/your/transcriptome"

# Directory to store the downloaded and processed files
OUTPUT_DIR="cellranger_analysis"

# Number of cores for Cell Ranger to use
LOCAL_CORES=8

# --- Script Start ---

echo "Starting SRA download and Cell Ranger analysis pipeline."

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if the SRA list file exists
if [ ! -f "$SRA_LIST" ]; then
    echo "Error: SRA list file not found at '$SRA_LIST'"
    exit 1
fi

# Read each SRA accession from the list
while IFS= read -r sra_accession || [[ -n "$sra_accession" ]]; do
    echo "--- Processing SRA: $sra_accession ---"

    # Define directories for this specific SRA
    sra_output_dir="$OUTPUT_DIR/$sra_accession"
    fastq_dir="$sra_output_dir/fastqs"
    cellranger_output_dir="$sra_output_dir/cellranger_count"

    # Create directories for the current SRA
    mkdir -p "$fastq_dir"
    mkdir -p "$cellranger_output_dir"

    # Download and dump to gzipped FASTQ
    echo "Downloading and converting $sra_accession to gzipped FASTQ..."
    fastq-dump --gzip --skip-technical --outdir "$fastq_dir" --split-files "$sra_accession"

    # Check if fastq-dump was successful
    if [ $? -ne 0 ]; then
        echo "Error: fastq-dump failed for $sra_accession. Skipping to next."
        continue
    fi

    echo "Download and conversion for $sra_accession complete."

    # Run Cell Ranger count
    echo "Running Cell Ranger for $sra_accession..."
    cellranger count --id="$sra_accession" \
                     --transcriptome="$TRANSCRIPTOME_PATH" \
                     --fastqs="$fastq_dir" \
                     --sample="$sra_accession" \
                     --localcores="$LOCAL_CORES" \
                     --output-dir="$cellranger_output_dir"

    # Check if Cell Ranger was successful
    if [ $? -ne 0 ]; then
        echo "Error: Cell Ranger failed for $sra_accession."
    else
        echo "Cell Ranger analysis complete for $sra_accession."
    fi

    echo "--- Finished processing $sra_accession ---"
    echo ""

done < "$SRA_LIST"

echo "All SRA accessions have been processed."