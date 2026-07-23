#!/bin/bash

# Input file containing GO IDs (one per line)
input_file="cleaned_go_ids_list.txt"
output_file="gene_symbols_from_go.tsv"

# Add header to output file
echo -e "GO_ID\tAccession\tProtein_Name\tGene_Symbol" > $output_file

# Loop through each GO ID and query UniProt
while read go_id; do
    echo "Querying GO ID: $go_id"
    curl -G "https://rest.uniprot.org/uniprotkb/search" \
        --data-urlencode "query=$go_id" \
        --data-urlencode "format=tsv" \
        --data-urlencode "fields=accession,protein_name,gene_primary" | \
        awk -v go_id="$go_id" 'NR > 1 {print go_id "\t" $0}' >> $output_file
done < $input_file

echo "Results saved to $output_file"

