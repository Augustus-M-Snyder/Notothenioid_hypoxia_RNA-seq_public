import pandas as pd
from collections import Counter

# Input files
ace_file = "./ACE_gene_symbol_aggregated_cleaned.csv"
cor_file = "./COR_protein_mapping.tsv"
geo_file = "./GEO_protein_mapping.tsv"
mac_file = "./MAC_gene_info.csv"

# Output file
output_file = "./species_gene_mapping_split_rows.csv"

# Load data
ace_data = pd.read_csv(ace_file)
cor_data = pd.read_csv(cor_file, sep="\t")
geo_data = pd.read_csv(geo_file, sep="\t")
mac_data = pd.read_csv(mac_file)

# Normalize to lowercase
ace_data["Gene_Names"] = ace_data["Gene_Names"].fillna("").str.lower()
ace_data["Synonyms"] = ace_data["Synonyms"].fillna("").str.lower()
cor_data["Gene_ID"] = cor_data["Gene_ID"].str.lower()
geo_data["Gene_ID"] = geo_data["Gene_ID"].str.lower()
mac_data["Gene_Name"] = mac_data["Gene_Name"].str.lower().fillna("NA")
mac_data["Gene_ID"] = mac_data["Gene_ID"].str.lower()

# Combine ACE gene symbols and synonyms into sets
ace_data["Synonyms"] = ace_data.apply(
    lambda row: set(row["Gene_Names"].split("; ") + row["Synonyms"].split("; ")),
    axis=1,
)

# Initialize the mapping output
final_mapping = []

# Process ACE data to generate rows
for _, row in ace_data.iterrows():
    ace_id = row["Gene_ID"]
    synonyms = row["Synonyms"]

    # Find matches in other species
    cor_matches = cor_data[cor_data["Gene_ID"].isin(synonyms)]["Gene_ID"].unique()
    geo_matches = geo_data[geo_data["Gene_ID"].isin(synonyms)]["Gene_ID"].unique()
    mac_matches = mac_data[mac_data["Gene_Name"].isin(synonyms)]["Gene_Name"].unique()
    mac_gene_ids = mac_data[mac_data["Gene_Name"].isin(synonyms)]["Gene_ID"].unique()

    # Split synonyms if they appear independently in other species
    matched_synonyms = set(cor_matches) | set(geo_matches) | set(mac_matches)
    independent_synonyms = synonyms & matched_synonyms
    remaining_synonyms = synonyms - independent_synonyms

    # Create rows for independent synonyms
    for synonym in independent_synonyms:
        cor_id = synonym if synonym in cor_matches else "NA"
        geo_id = synonym if synonym in geo_matches else "NA"
        mac_id = synonym if synonym in mac_matches else "NA"
        mac_gene_id = "; ".join(
            mac_data[mac_data["Gene_Name"] == synonym]["Gene_ID"].unique()
        ) if synonym in mac_matches else "NA"

        final_mapping.append({
            "Common_ID": synonym,
            "Synonyms": synonym,
            "ACE_ID": ace_id,
            "COR_ID": cor_id,
            "GEO_ID": geo_id,
            "MAC_ID": mac_id,
            "MAC_gene_ID": mac_gene_id,
        })

    # Create a row for remaining synonyms
    if remaining_synonyms:
        common_id = max(Counter(remaining_synonyms), key=lambda x: x)  # Most frequent or first alphabetically
        final_mapping.append({
            "Common_ID": common_id,
            "Synonyms": "; ".join(remaining_synonyms),
            "ACE_ID": ace_id,
            "COR_ID": "; ".join(cor_matches) if any(s in cor_matches for s in remaining_synonyms) else "NA",
            "GEO_ID": "; ".join(geo_matches) if any(s in geo_matches for s in remaining_synonyms) else "NA",
            "MAC_ID": "; ".join(mac_matches) if any(s in mac_matches for s in remaining_synonyms) else "NA",
            "MAC_gene_ID": "; ".join(mac_data[mac_data["Gene_Name"].isin(remaining_synonyms)]["Gene_ID"].unique()) if any(s in mac_matches for s in remaining_synonyms) else "NA",
        })

# Convert to DataFrame and save
final_mapping_df = pd.DataFrame(final_mapping)
final_mapping_df.to_csv(output_file, index=False)

print(f"Mapping completed and saved to {output_file}")

