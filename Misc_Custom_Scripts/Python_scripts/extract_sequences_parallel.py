import os
import re
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor

# Define input GTF files and corresponding genome files
gtf_files = {
    "ACE": ("ACE_annotation.gtf", "../Genomes/ACE_genome.fa"),
    "COR": ("COR_annotation.gtf", "../Genomes/COR_genome.fa"),
    "GEO": ("GEO_annotation.gtf", "../Genomes/GEO_genome.fa"),
    "MAC": ("MAC_annotation.gtf", "../Genomes/MAC_genome.fa"),
}

# Output directory for extracted sequences
output_dir = "sequence_blast_mapping"
os.makedirs(output_dir, exist_ok=True)

# Function to parse GTF and extract gene metadata
def parse_gtf(gtf_file):
    genes = []
    with open(gtf_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                if fields[2] == "gene":
                    attributes = fields[8]
                    gene_id = re.search(r'gene_id "([^"]+)"', attributes)
                    gene_name = re.search(r'gene_name "([^"]+)"', attributes)
                    genes.append({
                        "chromosome": fields[0],
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                        "strand": fields[6],
                        "gene_id": gene_id.group(1) if gene_id else None,
                        "gene_name": gene_name.group(1) if gene_name else None,
                    })
    return genes

# Function to extract gene sequences from genome in batches
def extract_sequences_batch(genes, genome_file, output_file, batch_size=500):
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    with open(output_file, "w") as out_f:
        for i in range(0, len(genes), batch_size):
            batch = genes[i:i + batch_size]
            for gene in batch:
                seq_id = gene["chromosome"]
                if seq_id in genome:
                    sequence = genome[seq_id].seq[gene["start"] - 1:gene["end"]]
                    if gene["strand"] == "-":
                        sequence = sequence.reverse_complement()
                    out_f.write(
                        f">{gene['gene_id']} {gene['gene_name']} {gene['chromosome']}:{gene['start']}-{gene['end']}\n"
                    )
                    out_f.write(str(sequence) + "\n")
                else:
                    print(f"Chromosome {seq_id} not found in genome file {genome_file}")
            print(f"Processed batch {i // batch_size + 1} for {output_file}")

# Process a single species
def process_species(species_abbr, gtf_file, genome_file):
    print(f"Processing {species_abbr}...")
    genes = parse_gtf(gtf_file)
    output_fasta = os.path.join(output_dir, f"{species_abbr}_genes.fasta")
    extract_sequences_batch(genes, genome_file, output_fasta)
    print(f"Extracted sequences for {species_abbr} saved to {output_fasta}")

# Run the processing in parallel
with ThreadPoolExecutor() as executor:
    futures = [
        executor.submit(process_species, species_abbr, gtf_file, genome_file)
        for species_abbr, (gtf_file, genome_file) in gtf_files.items()
    ]

