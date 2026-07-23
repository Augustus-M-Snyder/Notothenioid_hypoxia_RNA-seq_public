#!/usr/bin/bash

# Define the species array
species_array=("ACE" "GEO" "MAC")

# Define directories
chromosome_dir="/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/Genomes/Chromosome_Scaffolds"
annotation_dir="/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/Annotations"
output_dir="/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/Annotations/Subsets"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each species
for species in "${species_array[@]}"; do
  # Set filenames
  chromosome_list="${chromosome_dir}/${species}_chromosome_scaffolds.txt"
  annotation_file="${annotation_dir}/${species}_annotation.gtf"
  output_file="${output_dir}/${species}_subset.gtf"

  # Subset the GTF file using the chromosome list, including header lines
  awk 'BEGIN{while((getline < "'"$chromosome_list"'") > 0) chrs[$1]} /^#/ || ($1 in chrs)' "$annotation_file" > "$output_file"

  # Confirm output
  echo "Subset for $species saved to $output_file"
done

