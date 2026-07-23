import requests
import pandas as pd
import time

# Load the input file
input_file = "./ACE_gene_id_map_file.csv"
output_file = "./ACE_gene_symbol_mapping.csv"
data = pd.read_csv(input_file)

# Base URLs for NCBI Gene and UniProt
ncbi_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
uniprot_base_url = "https://rest.uniprot.org/uniprotkb/search"

# Function to query NCBI Gene
def query_ncbi_gene(query):
    params = {
        "db": "gene",
        "term": query,
        "retmode": "json",
    }
    response = requests.get(ncbi_base_url, params=params)
    if response.status_code == 200:
        results = response.json()
        if "esearchresult" in results and "idlist" in results["esearchresult"]:
            return results["esearchresult"]["idlist"]
    return []

# Function to query UniProt
def query_uniprot(query):
    params = {
        "query": query,
        "format": "tsv",
        "fields": "gene_primary",
    }
    response = requests.get(uniprot_base_url, params=params)
    if response.status_code == 200 and response.text.strip():
        lines = response.text.strip().split("\n")[1:]  # Skip header
        symbols = {line.split("\t")[0] for line in lines}
        return list(symbols)
    return []

# Initialize results list
results = []

# Query both databases for each row
for _, row in data.iterrows():
    original_id = row["Original_Gene_ID"]
    description = row["Description"]
    gene_symbol_ncbi = []
    gene_symbol_uniprot = []

    # Query NCBI Gene
    if pd.notna(original_id):
        gene_symbol_ncbi = query_ncbi_gene(original_id)
    if not gene_symbol_ncbi and pd.notna(description):
        gene_symbol_ncbi = query_ncbi_gene(description)

    # Query UniProt
    if pd.notna(original_id):
        gene_symbol_uniprot = query_uniprot(original_id)
    if not gene_symbol_uniprot and pd.notna(description):
        gene_symbol_uniprot = query_uniprot(description)

    # Merge results into the dataframe
    results.append({
        "Gene_ID_x": row["Gene_ID_x"],
        "Original_Gene_ID": original_id,
        "Description": description,
        "Gene_Symbol_NCBI": ";".join(gene_symbol_ncbi) if gene_symbol_ncbi else "NA",
        "Gene_Symbol_UniProt": ";".join(gene_symbol_uniprot) if gene_symbol_uniprot else "NA",
    })

    # Add a delay to avoid overwhelming the APIs
    time.sleep(0.5)

# Save the results to a new CSV file
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, index=False)

print(f"Gene symbol mapping saved to {output_file}")

