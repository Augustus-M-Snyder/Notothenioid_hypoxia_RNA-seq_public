import pandas as pd

# Load the mapping file
df = pd.read_csv("gene_id_mapping.csv")

# Prioritize rows with non-"NA" Original_Gene_ID
df_filtered = df.sort_values(by="Original_Gene_ID", na_position="last").drop_duplicates(subset="Gene_ID", keep="first")

# Remove the ".1" from the original gene id column
df["Original_Gene_ID"] = df["Original_Gene_ID"].str.replace(r"\.1$", "", regex=True)

# Save the filtered data
df_filtered.to_csv("gene_id_mapping_filtered.csv", index=False)

print("Filtered mapping file saved as gene_id_mapping_filtered.csv")

