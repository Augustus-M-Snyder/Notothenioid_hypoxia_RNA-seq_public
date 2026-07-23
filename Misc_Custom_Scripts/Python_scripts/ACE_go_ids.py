import pandas as pd

# Load the file
input_file = "./ACE_gene_id_map_file.csv"
data = pd.read_csv(input_file)

# Extract GO IDs and save them, splitting if necessary
go_ids = data["GO_IDs"].dropna().str.split(";").explode().drop_duplicates()

# Save to a file for searching
go_ids_file = "./ACE_go_ids_list.txt"
go_ids.to_csv(go_ids_file, index=False, header=False)

print(f"GO IDs saved to {go_ids_file}")

