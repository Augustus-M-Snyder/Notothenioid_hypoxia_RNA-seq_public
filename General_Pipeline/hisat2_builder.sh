#!/bin/bash
#SBATCH --partition=t1standard
#SBATCH --job-name=build_hgfm_indexes     # Job name
#SBATCH --nodes=3		 	  # Number of nodes
#SBATCH --ntasks=48                       # Number of threads
#SBATCH --tasks-per-node=16               # Number of tasks/threads per node  (MPI processes)
#SBATCH --cpus-per-task=1                 # Number of CPUs per task (cores)
#SBATCH --mem=128GB                       # Memory per node
#SBATCH --output=build_hgfm_indexes_%j.out  # Output log file
#SBATCH --error=build_hgfm_indexes_%j.err   # Error log file

# Load necessary modules (adjust based on your HPC environment)
module purge
module load hisat2   # Load HISAT2 module if necessary

# Load conda environment
conda init
source ~/.bashrc
conda activate NOTO_hyp

# Directories containing genome and annotation files
genome_dir="Genomes"
annotation_dir="Annotations"

# Define arrays for genome and annotation files
genome_files=("ACE_genome.fa" "COR_genome.fa" "GEO_genome.fa" "GIB_genome.fa" "MAC_genome.fa")
annotation_files=("ACE_annotation.gtf" "COR_annotation.gtf" "GEO_annotation.gtf" "MAC_annotation.gtf")

# Function to build HGFM index with transcripts
build_hgfm_index() {
    local species="$1"
    local genome_file="$2"
    local exon_file="$3"
    local ss_file="$4"

    echo "Building HGFM index with transcripts for $species"

    hisat2-build -p 16 --exon "$exon_file" --ss "$ss_file" "$genome_file" "${species}_genome_tran"

    echo "HGFM index built for $species"
    echo
}

# Iterate over each genome file
for genome_file in "${genome_files[@]}"; do
    species=$(basename "$genome_file" | cut -d '_' -f 1)

    # Check if there is a corresponding annotation file for the species
    annotation_file=""
    for anno_file in "${annotation_files[@]}"; do
        if [[ "$anno_file" == "${species}_annotation.gtf" ]]; then
            annotation_file="$anno_file"
            break
        fi
    done

    if [[ -n "$annotation_file" ]]; then
        echo "Processing species: $species"

        # Extract splice sites
        hisat2_extract_splice_sites.py "$annotation_dir/$annotation_file" > "$species.genome.ss"

        # Extract exons
        hisat2_extract_exons.py "$annotation_dir/$annotation_file" > "$species.genome.exon"

        # Build HGFM index
        build_hgfm_index "$species" "$genome_dir/$genome_file" "$species.genome.exon" "$species.genome.ss"
    else
        echo "No annotation file found for $species"
    fi
done

# End of script

