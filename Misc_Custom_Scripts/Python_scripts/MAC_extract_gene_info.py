import re
import pandas as pd

# Input and output file paths
gtf_file = "MAC_annotation.gtf"  # Replace with the path to your GTF file
output_file = "MAC_gene_info.csv"

# Initialize storage for gene information
gene_info = []

# Parse the GTF file
with open(gtf_file, "r") as gtf:
    for line in gtf:
        if line.startswith("#"):  # Skip header lines
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        attributes = fields[8]
        
        # Extract gene_id, gene_name, and gene_desc using regex
        gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
        gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
        gene_desc_match = re.search(r'gene_desc "([^"]+)"', attributes)
        
        if gene_id_match and gene_name_match and gene_desc_match:
            gene_id = gene_id_match.group(1)
            gene_name = gene_name_match.group(1)
            gene_desc = gene_desc_match.group(1)
            
            # Append the extracted information
            gene_info.append({"Gene_ID": gene_id, "Gene_Name": gene_name, "Gene_Desc": gene_desc})

# Convert to DataFrame and remove duplicates
df = pd.DataFrame(gene_info).drop_duplicates()

# Save to CSV
df.to_csv(output_file, index=False)

print(f"Gene information extracted and saved to {output_file}")

