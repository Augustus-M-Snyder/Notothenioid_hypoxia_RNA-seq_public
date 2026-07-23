import pandas as pd


# Load the output file
final_mapping_df = pd.read_csv("./species_gene_mapping.csv")

# Find rows where COR_ID, GEO_ID, or MAC_ID differ in the same row
diff_species_rows = final_mapping_df[
    (final_mapping_df["COR_ID"] != final_mapping_df["GEO_ID"]) |
    (final_mapping_df["COR_ID"] != final_mapping_df["MAC_ID"]) |
    (final_mapping_df["GEO_ID"] != final_mapping_df["MAC_ID"])
]

# Save to a new CSV for inspection
diff_species_rows.to_csv("./species_id_differences.csv", index=False)

print(f"Differences across species saved to species_id_differences.csv")
