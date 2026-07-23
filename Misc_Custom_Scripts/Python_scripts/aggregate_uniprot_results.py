import pandas as pd

# Input and output file paths
input_file = "./ACE_gene_symbol_results.csv"  # UniProt results
output_file = "./ACE_gene_symbol_aggregated_cleaned.csv"  # Cleaned and aggregated output

# Load the UniProt results
data = pd.read_csv(input_file)

# Helper function to normalize and deduplicate names
def normalize_and_deduplicate(names):
    # Split by `, `, normalize to lowercase, remove duplicates, and join back
    return "; ".join(sorted(set(name.lower() for name in names.split(", "))))

# Group by Gene_ID and aggregate unique values for relevant fields
aggregated_data = (
    data.groupby("Gene_ID")
    .agg({
        "Description": "first",  # Keep the first description (assuming it's the same for all)
        "Accession": lambda x: "; ".join(sorted(set(x.dropna()))),  # Combine unique accessions
        "Protein_Name": lambda x: "; ".join(sorted(set(x.dropna()))),  # Combine unique protein names
        "Gene_Names": lambda x: "; ".join(
            sorted(set(normalize_and_deduplicate(names) for names in x.dropna()))
        ),  # Normalize and deduplicate gene names
        "Synonyms": lambda x: "; ".join(
            sorted(set(normalize_and_deduplicate(names) for names in x.dropna()))
        ),  # Normalize and deduplicate synonyms
    })
    .reset_index()
)

# Save the cleaned and aggregated results to a new CSV file
aggregated_data.to_csv(output_file, index=False)

print(f"Cleaned and aggregated results saved to {output_file}")

