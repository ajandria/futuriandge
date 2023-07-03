#!/bin/bash

# Directory containing FastQ files
fastq_dir="/archive/tests/test_FastQ/"

# Output CSV file
output_file="samplesheet.csv"

# Write header to CSV file
echo "sample,fastq_1,fastq_2,strandedness" > $output_file

# Find FastQ files in directory
for fastq_1 in $(ls ${fastq_dir}*R1_001.fastq.gz); do
  # Get the corresponding R2 FastQ file
  fastq_2="${fastq_1/R1_001.fastq.gz/R2_001.fastq.gz}"

  # Check if R2 FastQ file exists
  if [ -f "$fastq_2" ]; then
    # Extract sample id from R1 FastQ file name
    sample_id=$(basename "$fastq_1" | cut -d "_" -f1)

    # Write sample information to CSV file
    echo "$sample_id,$fastq_1,$fastq_2,reverse" >> $output_file
  else
    echo "No matching R2 FastQ file for $fastq_1"
  fi
done
