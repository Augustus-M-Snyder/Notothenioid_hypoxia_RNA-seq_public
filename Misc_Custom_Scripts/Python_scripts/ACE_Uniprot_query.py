import requests
import csv
import pandas as pd
import time
import re

# Input and output file paths
input_file = "./ACE_gene_id_map_file.csv"  # Input file with gene IDs and descriptions
output_file = "./ACE_gene_symbol_results.csv"  # Output file to save results incrementally
completed_log_file = "./completed_queries.txt"  # File to log completed queries

# Base URL for UniProt API
base_url = "https://rest.uniprot.org/uniprotkb/search"

# Function to query UniProt for a gene description
def query_uniprot(description):
    query = f'protein_name:"{description}"'
    params = {
        "query": query,
        "fields": "accession,gene_names,protein_name",
        "format": "json",
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.json()  # Return the parsed JSON response
    else:
        print(f"Failed to fetch data for: {description} (Status Code: {response.status_code})")
        return None

# Load the log of completed queries
try:
    with open(completed_log_file, "r") as log:
        completed_queries = set(line.strip() for line in log.readlines())
except FileNotFoundError:
    completed_queries = set()

# Open the output file for incremental writing
with open(output_file, "a", newline="") as csvfile:
    fieldnames = ["Gene_ID", "Description", "Accession", "Protein_Name", "Gene_Names", "Synonyms"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # Write the header if the file is empty
    if csvfile.tell() == 0:
        writer.writeheader()
        csvfile.flush()

    # Load the input data
    data = pd.read_csv(input_file)

    # Query UniProt for each gene description
    for _, row in data.iterrows():
        gene_id = row["Gene_ID_x"]
        description = row["Description"]

        # Normalize description for case-insensitive checks
        normalized_description = description.lower() if pd.notna(description) else ""

        # Skip completed or invalid descriptions
        if (
            gene_id in completed_queries
            or not description
            or "predicted protein" in normalized_description
            or "predicted:" in normalized_description
            or "uncharacterized" in normalized_description
            or "---na---" in normalized_description
        ):
            print(f"Skipping invalid or completed entry: {description} (Gene_ID: {gene_id})")
            continue

        print(f"Querying UniProt for description: {description} (Gene_ID: {gene_id})")
        response = query_uniprot(description)

        if response:
            for entry in response.get("results", []):
                # Extract relevant fields
                accession = entry.get("primaryAccession", "NA")
                protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "NA")
                gene_names = ", ".join([gene.get("geneName", {}).get("value", "NA") for gene in entry.get("genes", [])])
                synonyms = ", ".join([syn.get("value", "NA") for gene in entry.get("genes", []) for syn in gene.get("synonyms", [])])

                # Exclude entries with LOC patterns
                if (
                    "predicted protein" in protein_name.lower()
                    or "predicted:" in protein_name.lower()
                    or "uncharacterized" in protein_name.lower()
                    or re.search(r"loc", gene_names.lower())
                    or re.search(r"loc", synonyms.lower())
                ):
                    continue

                # Write each result to the output file and flush
                writer.writerow({
                    "Gene_ID": gene_id,
                    "Description": description,
                    "Accession": accession,
                    "Protein_Name": protein_name,
                    "Gene_Names": gene_names,
                    "Synonyms": synonyms,
                })
                csvfile.flush()
        else:
            # Handle cases where no results are returned and flush
            writer.writerow({
                "Gene_ID": gene_id,
                "Description": description,
                "Accession": "NA",
                "Protein_Name": "NA",
                "Gene_Names": "NA",
                "Synonyms": "NA",
            })
            csvfile.flush()

        # Log the completed query
        with open(completed_log_file, "a") as log:
            log.write(f"{gene_id}\n")

        time.sleep(0.5)  # Rate-limit requests to avoid overwhelming the server

print(f"Results are being written incrementally to {output_file}")

