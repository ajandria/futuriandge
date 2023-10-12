#!/bin/bash

# Define the input CSV file from SRA
input_file="metadata.txt"

# Read each line from the CSV file, skipping the header
tail -n +2 "$input_file" | while IFS=, read -r run etc
do
    # Extract the GSM sample name from the line
    gsm_sample_name=$(echo "$etc" | cut -d',' -f17)

    # Run fasterq-dump for the current SRA ID
    fasterq-dump -v -p -x -O "$run" "$run"
    
    # Check if the first fastq file exists
    if [[ -f "$run/${run}_1.fastq" ]]; then
        # Create an empty GSM fastq file for R1 if it doesn't exist
        [[ ! -f "${gsm_sample_name}_R1.fastq" ]] && touch "${gsm_sample_name}_R1.fastq"
        # Append the content of the SRA fastq file to the GSM fastq file for R1
        cat "$run/${run}_1.fastq" >> "${gsm_sample_name}_R1.fastq"
    fi

    # Check if the second fastq file exists
    if [[ -f "$run/${run}_2.fastq" ]]; then
        # Create an empty GSM fastq file for R2 if it doesn't exist
        [[ ! -f "${gsm_sample_name}_R2.fastq" ]] && touch "${gsm_sample_name}_R2.fastq"
        # Append the content of the SRA fastq file to the GSM fastq file for R2
        cat "$run/${run}_2.fastq" >> "${gsm_sample_name}_R2.fastq"
    fi
done
