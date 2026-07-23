import pandas as pd

# Input file paths
mapping_file = "./gene_id_mapping_final.csv"
description_file = "./CA.gene.desc.xlsx"

# Load the gene mapping file
gene_mapping_df = pd.read_csv(mapping_file)

# Load the CA.gene.desc file
ca_desc_df = pd.read_excel(description_file)

# Extract only the base Gene ID from the SeqName column by removing '-RA'
ca_desc_df["Gene_ID_Base"] = ca_desc_df["SeqName"].str.split('-').str[0]

# Select relevant columns from CA.gene.desc
ca_desc_df_subset = ca_desc_df[["Gene_ID_Base", "Description", "GO IDs", "GO Names"]]
ca_desc_df_subset.columns = ["Gene_ID", "Description", "GO_IDs", "GO_Names"]

# Merge the two files on Original_Gene_ID from gene mapping and Gene_ID_Base from CA.gene.desc
merged_df = pd.merge(
    gene_mapping_df,
    ca_desc_df_subset,
    left_on="Original_Gene_ID",
    right_on="Gene_ID",
    how="left"
)

# Save the final output
output_file = "./gene_mapping_with_descriptions_fixed.csv"
merged_df.to_csv(output_file, index=False)

print(f"Mapping file saved as {output_file}")

