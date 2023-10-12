#!/bin/bash

# Set path to sra-tools binaries
export sra_tools_path=/Users/adrian/sratoolkit.3.0.7-mac64/bin

# Define the input CSV file
input_file=$1

# Read each line from the CSV file, skipping the header
tail -n +2 "$input_file" | while IFS=, read -r run etc
do
    # Extract the GSM sample name from the line
    # Step 1: Remove everything before GSM
    gsm_temp=$(echo "$etc" | sed 's/.*GSM/GSM/')
    # Step 2: Remove everything after the comma, including the comma
    gsm_id=$(echo "$gsm_temp" | sed 's/,.*//')

    echo "SRA id: ${run}"
    echo "GSM id: ${gsm_id}"

    # Run fasterq-dump for the current SRA ID
    ${sra_tools_path}/fasterq-dump -v -p -x -O "$run" "$run"
    
    # Check if the first fastq file exists
    if [[ -f "$run/${run}_1.fastq" ]]; then
        # Create an empty GSM fastq file for R1 if it doesn't exist
        [[ ! -f "${gsm_id}_R1.fastq" ]] && touch "${gsm_id}_R1.fastq"
        # Append the content of the SRA fastq file to the GSM fastq file for R1
        echo "[Dump Log] Dumping ${run}_1.fastq into ${gsm_id}_R1.fastq..."
        cat "$run/${run}_1.fastq" >> "${gsm_id}_R1.fastq"
    fi

    # Check if the second fastq file exists
    if [[ -f "$run/${run}_2.fastq" ]]; then
        # Create an empty GSM fastq file for R2 if it doesn't exist
        [[ ! -f "${gsm_id}_R2.fastq" ]] && touch "${gsm_id}_R2.fastq"
        # Append the content of the SRA fastq file to the GSM fastq file for R2
        echo "[Dump Log] Dumping ${run}_2.fastq into ${gsm_id}_R2.fastq..."
        cat "$run/${run}_2.fastq" >> "${gsm_id}_R2.fastq"
    fi
done

echo "Script finished for ${input_file}"
