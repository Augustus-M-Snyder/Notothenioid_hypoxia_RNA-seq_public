import re
import pandas as pd

# Input and output file paths
gtf_file = "ACE_annotation.gtf"
output_file = "gene_id_mapping.csv"

# Initialize storage for mapping
mapping = []

# Open and parse the GTF file
with open(gtf_file, "r") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue  # Skip header lines
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue  # Skip malformed lines
        
        # Extract the attributes column
        attributes = fields[8]
        
        # Extract gene_id
        gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
        if not gene_id_match:
            continue  # Skip if no gene_id is found
        gene_id = gene_id_match.group(1)
        
        # Extract orig_transcript_id
        orig_id_match = re.search(r'orig_transcript_id "gnl\|WGS:[^|]+\|([^"]+)"', attributes)
        if orig_id_match:
            orig_gene_id = orig_id_match.group(1)
        else:
            orig_gene_id = "NA"  # If no orig_transcript_id is found
        
        # Store the mapping
        mapping.append({"gene_id": gene_id, "orig_gene_id": orig_gene_id})

# Convert mapping to a DataFrame
df = pd.DataFrame(mapping)

# Remove duplicates (if multiple entries exist for the same gene_id)
df = df.drop_duplicates()

# Save to CSV
df.to_csv(output_file, index=False, header=["Gene_ID", "Original_Gene_ID"])

print(f"Mapping file created: {output_file}")

