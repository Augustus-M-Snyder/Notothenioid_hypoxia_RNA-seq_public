import pandas as pd

def extract_chromosome_locations(gtf_path, protein_mapping_path, output_path):
    # Load the protein mapping file
    protein_mapping = pd.read_csv(protein_mapping_path, sep='\t', header=0)
    gene_ids = protein_mapping.iloc[:, 0].tolist()

    # Function to extract details for all IDs from the GTF file
    def verify_gtf_entry(gtf_path, gene_id):
        positions = []
        with open(gtf_path, 'r') as gtf_file:
            for line in gtf_file:
                if gene_id in line and "gene_id" in line:
                    fields = line.strip().split('\t')
                    if len(fields) >= 9:
                        positions.append({
                            "chromosome": fields[0],
                            "start": int(fields[3]),
                            "end": int(fields[4]),
                            "strand": fields[6]
                        })
        return positions

    rows = []
    for gene_id in gene_ids:
        gene_positions = verify_gtf_entry(gtf_path, gene_id)
        if gene_positions:
            # Use the first entry (gene level) for chromosome details
            entry = gene_positions[0]
            protein_accession = protein_mapping[protein_mapping.iloc[:, 0] == gene_id].iloc[0, 1] if not protein_mapping.empty else "N/A"
            rows.append({
                "Gene ID": gene_id,
                "Protein Accession": protein_accession,
                "Chromosome": entry["chromosome"],
                "Start": entry["start"],
                "End": entry["end"],
                "Strand": entry["strand"]
            })

    # Create a DataFrame and save it as a CSV file
    output_df = pd.DataFrame(rows)
    output_df.to_csv(output_path, index=False)
    print(f"Chromosome information saved to {output_path}")

# Usage example
extract_chromosome_locations(
    gtf_path='./ACE_annotation.gtf',
    protein_mapping_path='./protein_mappings/ACE_protein_mapping.tsv',
    output_path='./ACE_gene_locations.csv'
)

