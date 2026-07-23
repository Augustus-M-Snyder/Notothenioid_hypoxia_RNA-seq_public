import requests
import pandas as pd

# Input and output files
go_ids_file = "cleaned_go_ids_list.txt"
output_file = "gene_symbols_from_go_fixed.csv"

# Load GO IDs
with open(go_ids_file, "r") as f:
    go_ids = [line.strip() for line in f if line.strip()]

# Base URL for UniProt API
base_url = "https://rest.uniprot.org/uniprotkb/search"

# Initialize results list
results = []

# Query UniProt for each GO ID
for go_id in go_ids:
    print(f"Querying UniProt for GO ID: {go_id}")
    params = {
        "query": f"annotation:(type:go {go_id})",
        "format": "tsv",
        "fields": "accession,protein_name,gene_primary",
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        lines = response.text.strip().split("\n")[1:]  # Skip header row
        for line in lines:
            accession, protein_name, gene_primary = line.split("\t")
            results.append({
                "GO_ID": go_id,
                "Accession": accession,
                "Protein_Name": protein_name,
                "Gene_Symbol": gene_primary,
            })
    else:
        print(f"Failed to fetch data for GO ID {go_id}: {response.status_code}")

# Convert results to DataFrame and save to CSV
df = pd.DataFrame(results)
df.to_csv(output_file, index=False)

print(f"Gene symbols saved to {output_file}")

