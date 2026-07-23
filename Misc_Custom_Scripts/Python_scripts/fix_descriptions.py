import pandas as pd

# Load the original file with descriptions
input_file = "./gene_mapping_with_descriptions_fixed.csv"

# Read the file into a DataFrame
data = pd.read_csv(input_file)

# Extract the 4th column (Description) as a list
description_list = data["Description"].tolist()

# Create a cleaned file with placeholders for invalid descriptions
output_file = "./description_list.txt"
with open(output_file, "w") as f:
    for description in description_list:
        if pd.isna(description) or description.strip() == "---NA---":
            f.write("NO_DESCRIPTION\n")
        else:
            f.write(description + "\n")

print(f"Cleaned description list saved to {output_file}")

