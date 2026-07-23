import csv

def extract_cds_with_gene_id(gtf_file, output_file):
    with open(gtf_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        # Write header
        writer.writerow(["chromosome", "start", "end", "strand", "gene_id", "transcript_id"])
        
        for line in infile:
            if line.startswith("#"):  # Skip comment lines
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != "CDS":  # Skip non-CDS rows
                continue
            
            chromosome = fields[0]
            start = fields[3]
            end = fields[4]
            strand = fields[6]
            
            # Extract gene_id and transcript_id from the attributes field
            attributes = fields[8]
            gene_id = None
            transcript_id = None
            for attr in attributes.split(';'):
                if "gene_id" in attr:
                    gene_id = attr.split('"')[1]
                if "transcript_id" in attr:
                    transcript_id = attr.split('"')[1]
            
            if gene_id and transcript_id:
                writer.writerow([chromosome, start, end, strand, gene_id, transcript_id])

# Input and output file paths
gtf_file = "GEO_annotation.gtf"  # Replace with your actual GTF file path
output_file = "GEO_CDS_mapping.tsv"  # Output file for CDS mappings

# Run the function
extract_cds_with_gene_id(gtf_file, output_file)
print(f"CDS entries with gene_id and transcript_id extracted to {output_file}")

