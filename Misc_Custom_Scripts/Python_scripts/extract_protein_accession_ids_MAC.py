import csv
from collections import defaultdict

def extract_accession_ids(gtf_file, output_file):
    # Dictionary to map gene_id to a list of protein IDs
    gene_to_proteins = defaultdict(set)
    
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue  # Skip headers/comments
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue  # Only process CDS rows
            
            # Extract gene_id and protein_id from the 9th column
            attributes = fields[8]
            gene_id = None
            protein_id = None
            for attr in attributes.split(";"):
                if "gene_id" in attr:
                    gene_id = attr.split('"')[1]
                if "protein_id" in attr:
                    protein_id = attr.split('"')[1]
            
            if gene_id and protein_id:
                gene_to_proteins[gene_id].add(protein_id)
    
    # Write results to output file
    with open(output_file, 'w') as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["Gene_ID", "Protein_Accession_IDs"])
        for gene_id, protein_ids in gene_to_proteins.items():
            writer.writerow([gene_id, ",".join(protein_ids)])

# Species-specific file processing
species_files = {
    "MAC": "MAC_annotation.gtf"
}

output_dir = "./protein_mappings"

import os
os.makedirs(output_dir, exist_ok=True)

for species, gtf_file in species_files.items():
    output_file = f"{output_dir}/{species}_protein_mapping.tsv"
    print(f"Processing {species}...")
    extract_accession_ids(gtf_file, output_file)
    print(f"Results saved to {output_file}")

