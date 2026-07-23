#!/bin/bash
#SBATCH --job-name=Fastqc_untrimmed_Emac
#SBATCH --partition=t1small
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=amsnyder7@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output_%j_%x.txt

# Function to calculate percentage completion
calculate_percentage() {
    local current=$1
    local total=$2
    echo $(awk "BEGIN { pc=100*${current}/${total}; i=int(pc); print (pc-i<0.5)?i:i+1 }")
}

# Source ~/.bashrc file and Activate conda environment
conda init
source ~/.bashrc
conda activate NOTO_hyp

# Check if there are any *.fastq.gz files
files=$(ls *.fastq.gz 2>/dev/null | wc -l)
if [ "$files" -eq 0 ]; then
    echo "No *.fastq.gz files found in the directory."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p /center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/untrimmed_reads/Emac_fastq_raw_reads/untrimmed_fastqc_output

# Count total number of files
total_files=$(ls *.fastq.gz | wc -l)
processed_files=0

# Run FastQC on each file
for file in *.fastq.gz; do
    echo "Running FastQC on $file"
    fastqc -o /center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/untrimmed_reads/Emac_fastq_raw_reads/untrimmed_fastqc_output "$file"

    # Update progress
    processed_files=$((processed_files + 1))
    percentage=$(calculate_percentage $processed_files $total_files)
    echo "Progress: $percentage%"
done
