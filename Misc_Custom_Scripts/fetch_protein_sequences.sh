#!/bin/bash

# Ensure efetch is available
if ! command -v efetch &> /dev/null; then
    echo "efetch command not found. Please ensure it is installed and in your PATH."
    exit 1
fi

# Array of species prefixes
species=("ACE" "COR" "GEO" "MAC" "STICK" "GIB")

# Loop through each species
for sp in "${species[@]}"; do
    input_file="${sp}_accessions.txt"
    output_file="${sp}_protein_sequences.fa"

    # Check if the input file exists
    if [[ -f "${input_file}" ]]; then
        echo "Fetching protein sequences for ${sp}..."
        efetch -input "${input_file}" -db protein -format fasta > "${output_file}"
        echo "Saved protein sequences to ${output_file}"
    else
        echo "Input file ${input_file} not found. Skipping ${sp}."
    fi
done

echo "Done fetching protein sequences."

