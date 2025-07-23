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

# --- Default Configuration ---
# These values are used if no parameters are provided via the command line.
SRA_LIST="sra_list.txt"
TRANSCRIPTOME_PATH="/path/to/your/transcriptome"
OUTPUT_DIR="cellranger_analysis"

# --- Functions ---
# Function to display usage information
usage() {
    echo "Usage: $0 -s <sra_list_file> -t <transcriptome_path> -o <output_dir>"
    echo ""
    echo "This script downloads SRA data and runs Cell Ranger count."
    echo ""
    echo "Options:"
    echo "  -s  Path to the text file containing a list of SRA accession numbers (one per line). Required."
    echo "  -t  Path to the Cell Ranger reference transcriptome. Required."
    echo "  -o  Directory to store the downloaded and processed files. Required."
    echo "  -h  Display this help message."
    exit 1
}

# --- Parse Command-Line Arguments ---
while getopts "s:t:o:h" opt; do
    case ${opt} in
        s)
            SRA_LIST=$OPTARG
            ;;
        t)
            TRANSCRIPTOME_PATH=$OPTARG
            ;;
        o)
            OUTPUT_DIR=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid Option: -$OPTARG" 1>&2
            usage
            ;;
        :)
            echo "Invalid Option: -$OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# --- Validate Required Parameters ---
if [ -z "$SRA_LIST" ] || [ -z "$TRANSCRIPTOME_PATH" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments."
    usage
fi


# --- Script Start ---
echo "--- Configuration ---"
echo "SRA List File: $SRA_LIST"
echo "Transcriptome Path: $TRANSCRIPTOME_PATH"
echo "Output Directory: $OUTPUT_DIR"
echo "---------------------"
echo ""
echo "Starting SRA download and Cell Ranger analysis pipeline."

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if the SRA list file exists
if [ ! -f "$SRA_LIST" ]; then
    echo "Error: SRA list file not found at '$SRA_LIST'"
    exit 1
fi

# Check if the transcriptome path exists
if [ ! -d "$TRANSCRIPTOME_PATH" ]; then
    echo "Error: Transcriptome directory not found at '$TRANSCRIPTOME_PATH'"
    exit 1
fi


# Read each SRA accession from the list
while IFS= read -r sra_accession || [[ -n "$sra_accession" ]]; do
    # Ignore empty lines or lines with only whitespace
    if [[ -z "${sra_accession// }" ]]; then
        continue
    fi

    echo "--- Processing SRA: $sra_accession ---"

    # Define directories for this specific SRA
    sra_output_dir="$OUTPUT_DIR/$sra_accession"
    fastq_dir="$sra_output_dir/fastqs"
    # Note: Cell Ranger creates its own output directory inside this path
    cellranger_base_dir="$sra_output_dir/cellranger_count" 

    # Create directories for the current SRA
    mkdir -p "$fastq_dir"
    mkdir -p "$cellranger_base_dir"

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
    # The actual output will be in a subfolder named after the --id
    cd "$cellranger_base_dir" && cellranger count --id="$sra_accession" \
                     --transcriptome="$TRANSCRIPTOME_PATH" \
                     --fastqs="$fastq_dir" \
                     --sample="$sra_accession"

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
